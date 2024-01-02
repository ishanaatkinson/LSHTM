
######################## LOAD PACKAGES #############################

library(pacman)
p_load(malariasimulation, foreSIGHT, malariaEquilibrium, tidyr, dplyr, ggplot2,
       reshape2, ggpubr, gridExtra, readxl, stringi, scene, XML, maps, readr,
       here, sf)

#  Set colour palette:
cols  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

################### IDENTIFY COUNTRY AND AREA CLASSIFICATION ##################

country_code <- "CMR" 

# Set WD 
# All model output dataframes will get saved in a folder called "simulation_results" within this work directory
setwd(paste0("C:/Users/IshanaAtkinson/OneDrive - London School of Hygiene and Tropical Medicine/Documents/R files/", as.character(country_code)))


set.seed(222)
################# READ FORESITE DATA AND GET PARAMETERS #######################

full_data <- readRDS(paste0(country_code, ".rds")) 

# convert all "-" into "_" to ease in analysis later on and remove any accents from admin-1 area names 
full_data$sites$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$sites$name_1), id = "Latin-ASCII")
full_data$prevalence$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$prevalence$name_1), id = "Latin-ASCII") 
full_data$interventions$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$interventions$name_1), id = "Latin-ASCII")
full_data$population$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$population$name_1), id = "Latin-ASCII") 
full_data$vectors$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$vectors$name_1), id = "Latin-ASCII")
full_data$pyrethroid_resistance$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$pyrethroid_resistance$name_1), id = "Latin-ASCII")
full_data$seasonality$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$seasonality$name_1), id = "Latin-ASCII")
full_data$eir$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$eir$name_1), id = "Latin-ASCII") 

full_data$sites$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$sites$name_2), id = "Latin-ASCII")
full_data$prevalence$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$prevalence$name_2), id = "Latin-ASCII") 
full_data$interventions$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$interventions$name_2), id = "Latin-ASCII")
full_data$population$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$population$name_2), id = "Latin-ASCII") 
full_data$vectors$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$vectors$name_2), id = "Latin-ASCII")
full_data$pyrethroid_resistance$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$pyrethroid_resistance$name_2), id = "Latin-ASCII")
full_data$seasonality$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$seasonality$name_2), id = "Latin-ASCII")
full_data$eir$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$eir$name_2), id = "Latin-ASCII") 


################## SET UP DATAFRAMES TO FILL #############################

if (rur_or_urb == "rural") {
  sim_output_rural <- c() # simulation output dataframe
  incidence_ppy_df_rural <- c() # incidence per person per year dataframe 
  incidence_df_rural <- c() # incidence dataframe
  population_df_rural <- c() # population size dataframe
  store_total_M_rural <- c() # mosquito population dataframe
}

if (rur_or_urb == "urban") {
  sim_output_urban <- c() # simulation output dataframe
  incidence_ppy_df_urban <- c() # incidence per person per year dataframe
  incidence_df_urban <- c() # incidence dataframe
  population_df_urban <- c() # population size dataframe
  store_total_M_urban <- c() # mosquito population dataframe
}
