#Simulating abundance commmunities to model
# Ana Miller-ter Kuile
# June 12, 2024

#this script simulates some abundance communities that will
#run faster than real communities for the MSAM model

# Load packages -----------------------------------------------------------


package.list <- c("tidyverse", 'here')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

set.seed(1)

# Attributes of the dataset -----------------------------------------------

#This model is structured with a few aspects of the data:
#1. Data are collected on species counts for multiple species
## in survey units (e.g., transects) that are clustered in larger
## sampling schemes (e.g., sites) across a series of years, 
## with 1+ replicate surveys for each site in each year (at 
## least a subset need to have 2+ surveys for this model structure)

#2. For this dataset, detection probabilities are dependent on a 
## survey covariate, such that in each site, year, replicate combo,
## a species detectability is dependent on some aspect of the environment
## at that time (e.g, weather, canopy cover, etc.) We have provided a 
## commented-out addition to the code that would enable addition of 
## a species covariate that would impact detectability, but don't include
## it here for ease of application

#3. As coded, not all transects have to have the same sampling years,
## such that you could incorporate a common pattern at long-term monitoring
## sites where some sites are surveyed for different periods of time. This
## model structure does assume you have yearly sampling of your communities.
## For a model where this isn't the case, please refer to the code for 
## the PFNP plant dataset in our main paper (LINK)

# Data needed for model ---------------------------------------------------

#here is a list of data objects we will need to run this model:

#Site.ID[t] - vector of larger site ids of length of number of transects
#Year.ID[y] - vector of year ids of length of number of years
#survey.covariate[t,y,r] - array of dimensions transect, year, replicate within year
#count[s,t,y,r] (response data) - array of dimensions species, transect, year, replicate
#n.species - number of species
#n.transects - number of transects
#n.sites - number of larger sites in survey design
#n.years - number of total years in the datset
#n.start[t] - vector of end years length number of transects 
#n.end[t] - vector of end years of length number of transects
#n.rep[t,y] - matrix of dimensions transects by years with number of surveys in that combo

# Specifics of our dataset ------------------------------------------------

#n.sites- 1
#n.transects - 3 (all in same site)
#n.years - 5
#n.species - 10
#n.rep - always 2 in each transect/year combo

# Generate data list ------------------------------------------------------

#numbers for indexing
n.sites <- 1
n.transects <- 3
n.years <- 5
n.species <- 5
n.rep <- matrix(data = 2, nrow = n.transects, ncol = n.years)
n.start <- rep(1, n.transects)
n.end <- rep(n.years, n.transects)

#random effect ID vectors
Site.ID <- c(rep(1, 3), rep(2, 3), rep(3, 3))
Site.ID <- c(rep(1,3))
Year.ID <- 1:n.years


# Create DF for all data --------------------------------------------------

species <- 1:n.species
transects <- 1:n.transects
years <- 1:n.years
reps <- 1:2

df <- expand.grid(species = species, 
                  transect = transects, 
                  year = years, 
                  rep = reps)

#get a survey covariate value for each of these transects, years, reps
df <- df %>%
  group_by(transect, year, rep) %>%
  mutate(survey.covariate = rnorm(1, 0, 1)) %>%
  ungroup()

#get a species intercept for each species for detection model
#and abundance model

mu.b0species <- rnorm(1, 0, 1)
sig.b0species <- runif(1, 0,5)
mu.a0 <- rnorm(1, 0, 1)
sig.a0 <- runif(1, 0,5)

df <- df %>%
  group_by(species) %>%
  mutate(a0 = rnorm(1, mu.a0, sig.a0),
         b0.species = rnorm(1, mu.b0species,sig.b0species)) %>%
  ungroup() 

#get a random effect of site and year within each species
df <- df %>%
  mutate(site = case_when(transect %in% c(1:3) ~ 1,
                          transect %in% c(4:6) ~ 2,
                          transect %in% c(7:9) ~ 3,
                          transect %in% c(10:12) ~ 4,
                          transect %in% c(13:15) ~ 5,
                          transect %in% c(16:18) ~ 6,
                          transect %in% c(19:21) ~ 7,
                          transect %in% c(22:24) ~ 8,
                          transect %in% c(25:27) ~ 9,
                          transect %in% c(28:30) ~ 10,
                          TRUE ~ NA_real_)) %>%
  group_by(species, site) %>%
  mutate(eps.site = rnorm(1, 0,1)) %>%
  ungroup() %>%
  group_by(species, year) %>%
  mutate(eps.year = rnorm(1, 0,1)) %>%
  ungroup()

# Get p, N, and y in df ---------------------------------------------------

#y is dependent on p and N
#p is based on a regression:
#logit(p[s,t,y,r]) <- a0[s] + a1*survey.covariate[t,y,r]

#N is dependent on a lambda:
# N[s,t,y] ~ dpois(lambda[s,t,y])

#which is dependent on a regression:
#log(lambda[s,t,y]) <- b0.species[s] + 
#  eps.site[Site.ID[t], s] +
#  eps.year[Year.ID[y], s]

#y[s,t,y,r] (response data) - array of dimensions species, site, year, replicate
#y[s,t,y,r] ~ dbin(p[s,t,y,r], N[s,t,y])

a1 <- 0.2 #set to a value

df <- df %>%
  rowwise() %>%
  mutate(logit.p = a0 + a1*survey.covariate,
         p = plogis(logit.p),
         log.lambda = b0.species + eps.site + eps.year,
         lambda = exp(log.lambda),
         N = rpois(1, lambda),
         count = rbinom(1, N, p))

saveRDS(df, here("MSAM",
                 "data",
                 "simulated_data",
                 "MSAM_simulated_community_data.RDS"))

# Pull out covariate and response data --------------------------------------------------

#y[s,t,y,r] (response data) - array of dimensions species, transect, year, replicate
#survey.covariate[t,y,r] - array of dimensions transect, year, replicate within year

spec <- df$species
trans <- df$transect
yr <- df$year
rp <- df$rep

#make a blank array with dims of species x years
count <- array(NA, dim = c(n.species, n.transects, n.years, 2))

#fill taht array based on the values in those columns
for(i in 1:dim(df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the species of row i,
  #the transect of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the observed number
  # for that speciesxtransectxyearxreplicate combo
  count[spec[i], trans[i], yr[i], rp[i]] <- as.numeric(df[i,16])
}

#get the survey covariate array
survey.covariate <- array(dim = c(n.transects, n.years, 2))
for(i in 1:dim(df)[1]){ #dim[1] = n.rows
  #using info from the dataframe on the transect of row i,
  # the year of row i and the replicate of row i,
  # populate that space in the array with the column in
  # the dataframe that corresponds to the covariate value
  # for that transectxyearxreplicate combo
  survey.covariate[trans[i], yr[i], rp[i]] <- as.numeric(df[i,5])
}


# Initial values ----------------------------------------------------------

#Need to add max count to set initial value for N!!! TODOOOOO
Ndf <- df %>%
  group_by(transect, year, species) %>%
  summarise(maxcount = max(count, na.rm = T))

#now, generate IDs for the for loop where
# we will populate the matrix
yr2 <- Ndf$year #get a yearID for each iteration of the loop
trans2 <-Ndf$transect #site ID for each iteration fo the loop
spec2 <- Ndf$species #get a species ID for each iteration of the loop

#make a blank array with dims of species x years
countmax <- array(NA, dim = c(n.species, n.transects, n.years))

#fill taht array based on the values in those columns

for(i in 1:dim(Ndf)[1]){ #dim[1] = n.rows
  countmax[spec2[i], trans2[i], yr2[i]] <- as.numeric(Ndf[i,4])
}

#set all zeros (which could be true or false) to NA
countmax[countmax == 0] <- NA

a0 <- df %>%
  distinct(a0, species) %>%
  dplyr::select(a0) %>%
  as_vector()

b0.species <- df %>%
  distinct(species, b0.species) %>%
  dplyr::select(b0.species) %>%
  as_vector()

eps.site <- df %>%
  distinct(species, site, eps.site) %>%
  pivot_wider(names_from = 'species',
              values_from = eps.site) %>%
  column_to_rownames(var = 'site') %>%
  as.matrix()
  
eps.year <- df %>%
  distinct(species, year, eps.year) %>%
  pivot_wider(names_from = 'species',
              values_from = eps.year) %>%
  column_to_rownames(var = 'year') %>%
  as.matrix()

# Create data list --------------------------------------------------------

data_list <- list(Site.ID = Site.ID, 
                  Year.ID = Year.ID,
                  survey.covariate = survey.covariate,
                  count = count,
                  n.species = n.species,
                  n.transects = n.transects,
                  n.sites = n.sites,
                  n.years = n.years,
                  n.start = n.start,
                  n.end = n.end,
                  n.rep = n.rep)

saveRDS(data_list, here("MSAM",
                        'data',
                        'model_inputs',
                        'MSAM_data_list.RDS'))

inits_list <- list(a0 = a0,
                   a1 = a1,
                   b0.species = b0.species,
                   eps.site = eps.site,
                   eps.year = eps.year,
                   countmax = countmax,
                   mu.b0species = mu.b0species,
                   mu.a0 = mu.a0)

saveRDS(inits_list, here("MSAM", 
                         'data',
                         'model_inputs',
                         'MSAM_inits_list.RDS'))
