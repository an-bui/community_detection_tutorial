#MSAM wrapper script

# Load packages -----------------------------------------------------------


package.list <- c("tidyverse", 'here', #general packages
                  'jagsUI', #jags wrapper
                  'coda', #gelman.diag() function
                  'mcmcplots') #trace plot function

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data -----------------------------------------------------------

data_list <- readRDS(here("MSAM",
                          'data',
                          'model_inputs',
                          'MSAM_data_list.RDS'))

# Define model path -----------------------------------------------------------

#The model file is here:
model_file <- here("MSAM",
                   'code',
                   'MSAM_model.R')

# Specify parameters to save -----------------------------------------------------------

#These are the parameters we will want to track 
#to assess convergence
parameters <- c('b0.star',
                'a0',
                'a1',
                'eps.site.star',
                'eps.year.star',
                'mu.a0',
                'sig.a0',
                'mu.b0species',
                'sig.b0species',
                'sig.eps.site',
                'sig.eps.year')

# Set initial values for model -----------------------------------------------------------

#this was also generated when we simulated data - it aids in 
#convergence when we have known values
inits_list <- readRDS(here("MSAM",
                           'data',
                           "model_inputs",
                           'MSAM_inits_list.RDS'))

#IMPORTANT: You will need to set initial values for the number of 
#individuals (N), which will keep the model from running into errors
#when a sampled number in the distribution exceeds the maximum
#number observed at that site

#We can also set initial values for various variables 
# with known values to get convergence faster, but these are not
#necessary to get the model to run
inits <- list(list(a1 = inits_list$a1,
                   mu.a0 = inits_list$mu.a0,
                   mu.b0species = inits_list$mu.b0species,
                   a0 = inits_list$a0,
                   b0.species = inits_list$b0.species,
                   eps.site = inits_list$eps.site,
                   eps.year = inits_list$eps.year,
                   N = inits_list$countmax),
              list(a1 = inits_list$a1 + 0.05,
                   mu.a0 = inits_list$mu.a0 + 0.05,
                   mu.b0species = inits_list$mu.b0species + 0.05,
                   a0 = inits_list$a0 + 0.05,
                   b0.species = inits_list$b0.species + 0.05,
                   eps.site = inits_list$eps.site + 0.05,
                   eps.year = inits_list$eps.year + 0.05,
                   N = inits_list$countmax),
              list(a1 = inits_list$a1 - 0.05,
                   mu.a0 = inits_list$mu.a0 - 0.05,
                   mu.b0species = inits_list$mu.b0species - 0.05,
                   a0 = inits_list$a0 - 0.05,
                   b0.species = inits_list$b0.species - 0.05,
                   eps.site = inits_list$eps.site - 0.05,
                   eps.year = inits_list$eps.year - 0.05,
                   N = inits_list$countmax))

# Run the model -----------------------------------------------------------

#run the model 
model <- jagsUI::jags(data = data_list,
                      inits = inits,
                      parameters.to.save = parameters,
                      model.file = model_file,
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 1335,
                      DIC = TRUE)

# Assess convergence -----------------------------------------------------------

#assess convergence with Gelman statistic
gelman.diag(model$samples, multivariate = F)
#generate trace plots
mcmcplot(model$samples)