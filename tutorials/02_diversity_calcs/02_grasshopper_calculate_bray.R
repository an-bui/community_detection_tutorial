#Get Bray-Curtis out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises Bray Curtis for visualizations
#purposes for An's project

# Load packages ---------------------------------------------------------------
Sys.time()

# Load packages,
package.list <- c('here', 'tidyverse') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

meta <- read.csv(here('examples',
                      'data_output',
                      'grasshopper_stability',
                      'tidy_data',
                      'grasshopper_msam_tidy_data.csv'))
  
#load the formatted data for the JAGS model
data <- readRDS(here('examples',
                     'data_output',
                     'grasshopper_stability',
                     'model_inputs',
                     'grasshopper_msam_input_data_list.RDS'))

samples <- readRDS(here('examples',
                        'data_output',
                        'grasshopper_stability',
                        'computing_cluster_outputs',
                        'sev_N_samples.RDS'))

# Get dataframe out of samples --------------------------------------------

ndf <- as.data.frame.table(samples$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         speciesID = as.numeric(as.factor(Var2)),
         siteID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, speciesID, siteID, yrID, N)


hopper_sum <- meta %>%
  filter(CNT > 0) %>%
  group_by(speciesID) %>%
  summarise(max = max(CNT, na.rm = T),
            sd = sd(CNT, na.rm = T),
            var = var(CNT, na.rm = T)) %>%
  mutate(sd = case_when(is.na(sd) ~ max,
                        TRUE ~ sd),
         var = case_when(is.na(var) ~ max,
                         TRUE ~ var))

ndf2 <- ndf %>%
  left_join(hopper_sum, by = "speciesID")%>%
  mutate(N = case_when(N > (max+var) ~ rnorm(1, (max+sd), sd),
                       TRUE ~ N)) 

n.transects <- data$n.transects
n.start <- data$n.start
n.end <- data$n.end
n.species <- data$n.species
n.years <- data$n.years

bray_fun <- function(iter){
  
  ndf3 <- ndf2 %>%
    filter(Iteration == iter)
  
  #make a blank array with dims of species x years x reps
  N <- array(NA, dim = c(n.species, n.transects, n.years))
  
  spec <- ndf3$speciesID
  site <- ndf3$siteID
  yr <- ndf3$yrID
  
  #fill taht array based on the values in those columns
  # for occupancy
  for(i in 1:dim(ndf3)[1]){ #dim[1] = n.rows
    #using info from the dataframe on the species of row i,
    #the site of row i, 
    # the year of row i and the replicate of row i,
    # populate that space in the array with the column in
    # the dataframe that corresponds to the 1-0 occupancy
    # for that speciesxyearxreplicate combo
    N[spec[i], site[i], yr[i]] <- as.numeric(ndf3[i,5])
  }
  
#spec,site,yr
a <- array(NA, dim = c(n.species, n.transects, n.years))
b <- array(NA, dim = c(n.species, n.transects, n.years))
c <- array(NA, dim = c(n.species, n.transects, n.years))
#site, yr
A <- matrix(NA, nrow = n.transects, ncol = n.years)
B<- matrix(NA, nrow = n.transects, ncol = n.years)
C <- matrix(NA, nrow = n.transects, ncol = n.years)
num <- matrix(NA, nrow = n.transects, ncol = n.years)
denom1 <- matrix(NA, nrow = n.transects, ncol = n.years)
denom <- matrix(NA, nrow = n.transects, ncol = n.years)
bray <- matrix(NA, nrow = n.transects, ncol = n.years)
num.bal <- matrix(NA, nrow = n.transects, ncol = n.years)
denom.bal1 <- matrix(NA, nrow = n.transects, ncol = n.years)
denom.bal <- matrix(NA, nrow = n.transects, ncol = n.years)
bray_balanced <- matrix(NA, nrow = n.transects, ncol = n.years)
bray_gradient <- matrix(NA, nrow = n.transects, ncol = n.years)


#BRAY CURTIS DERIVED QUANTIIES
#lots of ways to calculate this, but I did this way
#IF WE WANT TO partition components of Bray:
for(i in 1:n.transects){
  for(t in (n.start[i]+1):n.end[i]){
    for(k in 1:n.species){
      # num individuals in both time periods per species
      a[k,i,t] <- min(N[k,i,t-1], N[k,i,t])
      # num individuals only in first time point
      b[k,i,t] <- N[k,i,t-1] - a[k,i,t]
      # num individuals only in second time point
      c[k,i,t] <- N[k,i,t] - a[k,i,t]
    }
    #for all years 2 onward:
    #total number of shared individuals across time periods
    A[i,t] <- sum(a[,i,t])
    #total number of individuals in only first time period
    B[i,t] <- sum(b[,i,t])
    #total number of individuals in only second time period
    C[i,t] <- sum(c[,i,t])
    
    #total bray-curtis (B+C)/(2A+B+C)
    num[i,t] <- B[i,t] + C[i,t]
    denom1[i,t] <- 2*A[i,t]+B[i,t]+C[i,t]
    #if all values are zero - this just keeps the eqn. from
    #dividing by zero
    denom[i,t] <- ifelse(denom1[i,t]==0,1, denom1[i,t])
    bray[i,t] <- num[i,t]/denom[i,t]
    
    #how much is dissimilarity shaped by
    # individuals of one species being replaced by individuals
    #of another species?
    num.bal[i,t] <- min(B[i,t], C[i,t])
    denom.bal1[i,t] <- A[i,t] + num.bal[i,t]
    denom.bal[i,t] <- ifelse(denom.bal1[i,t] == 0,1, denom.bal1[i,t])
    bray_balanced[i,t] <- num.bal[i,t]/denom.bal[i,t]
    
    #how much is dissimilarity shaped by
    # individuals that are lost without substitution?
    bray_gradient[i,t] <- bray[i,t] - bray_balanced[i,t]
  }
}

bray1 <- as.data.frame(bray) %>%
  rownames_to_column(var = "siteID") %>%
  pivot_longer(cols = 2:last_col(),
               names_to = 'yrID',
               values_to = 'bray')

balanced1 <- as.data.frame(bray_balanced) %>%
  rownames_to_column(var = "siteID") %>%
  pivot_longer(cols = 2:last_col(),
               names_to = 'yrID',
               values_to = 'bray_balanced')

gradient1 <- as.data.frame(bray_gradient) %>%
  rownames_to_column(var = "siteID") %>%
  pivot_longer(cols = 2:last_col(),
               names_to = 'yrID',
               values_to = 'bray_gradient')

all_df <- bray1 %>%
  left_join(balanced1, by = c("siteID", "yrID")) %>%
  left_join(gradient1, by = c("siteID", "yrID")) %>%
  mutate(Iteration = iter) %>%
  mutate(yrID = str_sub(yrID, 2, nchar(yrID))) %>%
  mutate(yrID = as.numeric(yrID),
         siteID = as.numeric(siteID))

return(all_df)
}

iter_set <- 1:max(ndf$Iteration)

bray_list <- lapply(iter_set, bray_fun)

bray_df <- bind_rows(bray_list)

saveRDS(bray_df, here('examples',
                      'data_output',
                      'grasshopper_stability',
                      'tidy_data',
                      'grasshopper_modeled_bray_samples.RDS'))

bray_df_sum <- bray_df %>%
  group_by(siteID, yrID) %>%
  summarise(mean_bray = mean(bray, na.rm = T),
            sd_bray = sd(bray, na.rm = T),
            mean_balanced = mean(bray_balanced, na.rm = T),
            sd_balanced = sd(bray_balanced, na.rm = T),
            mean_gradient = mean(bray_gradient, na.rm = T),
            sd_gradient = sd(bray_gradient, na.rm = T))

saveRDS(bray_df_sum, here('examples',
                      'data_output',
                      'grasshopper_stability',
                      'tidy_data',
                      'grasshopper_stability_meanSD_corrected.RDS'))
