#Prepping data object for the stability SAM models
#Ana Miller-ter Kuile
#July 5, 2023

#this is a script that preps data list for the turnover SAM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

all_data <- read.csv(here("SAM",
                          'data',
                          "tidy_data",
                          "SBC_fish_stability_metrics_with_covariates.csv"))

# Prep data for jags ------------------------------------------------------

#indexing for model main likelihood loop
n.data <- nrow(all_data)

#indexing for site-level intercept random effects 
n.sites <- length(unique(all_data$SITE))

#number of plant lags
n.plantlag <- all_data %>%
  dplyr::select(DRY_GM2:DRY_GM2_l5) %>%
  ncol()

#number of temperatuere lags
n.templag <- all_data %>%
  dplyr::select(TEMP_C:TEMP_C_l10) %>%
  ncol()


# Response data -----------------------------------------------------------

#get d and the known variance
d <- as.vector(all_data$bray)
var.estimate <- as.vector(all_data$SD^2)


# Random effect and covariaes ---------------------------------------------

#get these site IDs 
Site.ID <- all_data %>%
    mutate(SITE = as.numeric(as.factor(SITE))) %>%
    dplyr::select(SITE) %>%
    as_vector()
  
#plant covariate data (needs to be a matrix of sites x lags)
Plant <- all_data %>%
  dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>%
  pivot_longer(DRY_GM2:DRY_GM2_l5,
               names_to = 'lag',
               values_to = 'kelp') %>%
  mutate(kelp = scale(kelp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "kelp") %>%
  dplyr::select(DRY_GM2:DRY_GM2_l5) %>%
  as.matrix()

#temperature covariate data (needs to be a matrix of sites x lags)
Temp <- all_data %>%
  dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l10) %>%
  pivot_longer(TEMP_C:TEMP_C_l10,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(TEMP_C:TEMP_C_l10) %>%
  as.matrix()

# Compile data list -------------------------------------------------------

data <- list(n.data = n.data,
             n.sites = n.sites,
             Site.ID = Site.ID,
             d = d,
             var.estimate = var.estimate,
             n.plantlag = n.plantlag,
             Plant = Plant,
             n.templag = n.templag,
             Temp = Temp)

saveRDS(data, here('SAM',
                   'data',
                   "model_inputs",
                   "SAM_input_data.RDS"))


