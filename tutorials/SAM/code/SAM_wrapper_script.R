
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", #general packages for data input/manipulation
                  "jagsUI", #to run JAGS models
                  'mcmcplots', #to look at trace plots
                  "coda") #to evaluate convergence
                   

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('SAM',
                          'data',
                          "model_inputs",
                          "SAM_input_data.RDS"))


# Define model path -------------------------------------------------------

model <- here('SAM',
              'code',
              'SAM_model.R')

# Parameters to save ------------------------------------------------------

#these parameters we can track to assess convergence
params <- c('b0.site', #site-level intercepts
            'b0', #overall intercept
            'b', #covariate effects
            'wA', #plant biomass weights vector
            'wB', #temperature weights
            'sig.site', #sd of site effects
            'var.process') #unknown variance



# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 1500,
                    n.iter = 5000,
                    DIC = TRUE)


# Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)


# Get the parameter summaries ---------------------------------------------

sum <- summary(mod$samples)
# 
saveRDS(sum, here('SAM',
                  'data',
                  'model_outputs',
                  'SAM_model_summary.RDS'))

