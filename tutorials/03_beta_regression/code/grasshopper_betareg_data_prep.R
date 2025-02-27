#Sevilleta grasshoppers environmental variable prep
# October 12, 2023

#this script preps the environmental variables with lags
#for the Sevilleta Grasshopper dataset


#Sites for Sevilleta grasshoppers, based on temporal coverage:

#NPP is from 1999 for two sites: 
#Core-Blak gramma and Core-Creosote
#these sites are also called "Five Points" in the datasets
#I think for meteorological station this includes station:
#49 (Five points) 

#NPP are seasonal - spring and fall of each year
#looks like grouping them by "web" and then taking an average of 
#the plots on that web for all species biomass in each season
#would be the way to combine with grasshopper data

#Climate are hourly - only one site near both the creosote
#and black gramma sites
#I think for meteorological station this includes station:
#49 (Five points) 

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

grasshoppers <- read.csv(here('tutorials',
                              '03_beta_regression',
                              'data',
                              'tidy_data',
                              'grasshopper_betareg_tidydata.csv'))


# Prep data for jags ------------------------------------------------------


# prep modeled data for jags ----------------------------------------------

# Filter years with data on climate/npp -----------------------------------

all_data3 <- grasshoppers %>%
  filter(YEAR > 1998) %>%
  #make site-web ID
  unite(c(site, web),
        col = "site_web",
        remove = F,
        sep = "_") %>%
  mutate(Web.ID = as.numeric(as.factor(site_web)))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(all_data3)

n.webs <- length(unique(all_data3$site_web))

n.transects <- length(unique(all_data3$site_web_trans))

n.templag <- all_data3 %>%
  dplyr::select(Temp:Temp_l5) %>%
  ncol()

n.pptlag <- all_data3 %>%
  dplyr::select(PPT:PPT_l5) %>%
  ncol()

n.npplag <- all_data3 %>%
  dplyr::select(NPP:NPP_l10) %>%
  ncol()

# Response data -----------------------------------------------------------

bray <- as.vector(all_data3$mean_bray)
var.estimate <- as.vector(all_data3$sd_bray^2)

# Random effects ----------------------------------------------------------

# Web.ID <- all_data3 %>%
#   distinct(siteID, Web.ID) %>%
#   dplyr::select(Web.ID) %>%
#   as_vector()

Web.ID <- as.vector(all_data3$Web.ID)

Transect.ID <- as.vector(all_data3$siteID)

# Covariates --------------------------------------------------------------

Temp <- all_data3  %>%
  dplyr::select(site_web_trans, YEAR, Temp:Temp_l5) %>%
  pivot_longer(Temp:Temp_l5,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(Temp:Temp_l5) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))
#7% missing data
#with more lags - 0.4%

PPT <- all_data3  %>%
  dplyr::select(site_web_trans, YEAR, PPT:PPT_l5) %>%
  pivot_longer(PPT:PPT_l5,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l5) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#<5% missing data

NPP <- all_data3 %>%
  dplyr::select(site_web_trans, YEAR, NPP:NPP_l10) %>%
  pivot_longer(NPP:NPP_l10,
               names_to = 'lag',
               values_to = 'npp') %>%
  mutate(npp = scale(npp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "npp") %>%
  dplyr::select(NPP:NPP_l10) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))
#~10% missing data


# Combine data into a data list -------------------------------------------


data1 <- list(n.data = n.data,
             n.webs = n.webs,
             n.templag = n.templag,
             n.npplag = n.npplag,
             n.pptlag = n.pptlag,
             n.transects = n.transects,
             bray = bray,
             var.estimate = var.estimate,
             Web.ID = Web.ID,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT, 
             NPP = NPP)

saveRDS(data1, here('tutorials',
                    '03_beta_regression',
                    'data',
                    'model_inputs',
                    'grasshopper_betareg_input_data_list.RDS'))

