#source("CMR_model.R")
#source("running_model_CMR.R")

###################### READ IN RESULTS FROM IMPERIAL MODEL ####################

#sim_output_rural <- read.csv(paste0(getwd(), "/simulation_results/sim_output_rural.csv"))
incidence_ppy_df_rural <- read.csv(paste0(getwd(), "/simulation_results/incidence_ppy_df_rural.csv"))
incidence_df_rural <- read.csv(paste0(getwd(), "/simulation_results/incidence_df_rural.csv"))
population_df_rural <- read.csv(paste0(getwd(), "/simulation_results/population_df_rural.csv"))

#sim_output_urban <- read.csv(paste0(getwd(), "/simulation_results/sim_output_urban.csv"))
incidence_ppy_df_urban <- read.csv(paste0(getwd(), "/simulation_results/incidence_ppy_df_urban.csv"))
incidence_df_urban <- read.csv(paste0(getwd(), "/simulation_results/incidence_df_urban.csv"))
population_df_urban <- read.csv(paste0(getwd(), "/simulation_results/population_df_urban.csv"))


############# MERGE DATAFRAMES TOGETHER BY TAKING WEIGHTED AVERAGE ############# 

########################## SIMULATION OUTPUT ########################## 
# 
# merge_sim_output <- bind_rows(sim_output_rural, sim_output_urban) %>% 
#   group_by(area, timestep)
# 
# # USING WEIGHTED MEAN 
# merged_sim_output_df <- data.frame()
# 
# for (i in 1:length(area_names)) {
#   current_merge_sim_output <- merge_sim_output %>% filter(area==area_names[i])
#   for (j in 1:sim_length) {
#     current_merge_sim_output_by_time <- current_merge_sim_output %>% filter(timestep == j)
#     
#     if (dim(current_merge_sim_output_by_time)[1] == 2){
#       rural_df <- (current_merge_sim_output_by_time %>% filter(rur_or_urb == "rural"))[,5:dim(current_merge_sim_output_by_time)[2]]
#       urban_df <- (current_merge_sim_output_by_time %>% filter(rur_or_urb == "urban"))[,5:dim(current_merge_sim_output_by_time)[2]]
#       
#       rural_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "rural"))$pop
#       urban_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "urban"))$pop
#       
#       rural_weight <- rural_pop / (sum(rural_pop, urban_pop))
#       urban_weight <- urban_pop / (sum(rural_pop, urban_pop))
#       
#       weighted_means <- (rural_weight * rural_df) + (urban_weight * urban_df)
#       
#     } else {
#       weighted_means <- current_merge_sim_output_by_time[,5:dim(current_merge_sim_output_by_time)[2]]
#     }
#     
#     merged_sim_output_df <- rbind(merged_sim_output_df, weighted_means)
#     
#     
#   }
# }
# 
# names(merged_sim_output_df) <- colnames(merge_sim_output[,5:dim(merge_sim_output)[2]])
# merged_sim_output_df$timestep <- rep(1:sim_length, times=length(area_names))
# merged_sim_output_df$iso_code <- rep(country_code, dim(merged_sim_output_df)[1])
# merged_sim_output_df$rur_or_urb <- rep("rural/urban merged", dim(merged_sim_output_df)[1])
# merged_sim_output_df$area <- rep(area_names, each = sim_length)
# merged_sim_output_df <- merged_sim_output_df %>%
#   relocate(timestep, iso_code, rur_or_urb, area)


########################## INCIDENCE PPY ########################## 

merge_incidence_ppy <- bind_rows(incidence_ppy_df_rural, incidence_ppy_df_urban) %>% 
  group_by(area, infection_class, age_in_days_midpoint)

infection_classes <- c("clinical", "severe", "total", "asymptomatic")
age_midpoint <- age_in_days_midpoint

# USING WEIGHTED MEAN 
merged_incidence_ppy_df <- data.frame()

for (i in 1:length(area_names)) {
  current_merge_incidence_ppy <- merge_incidence_ppy %>% filter(area==area_names[i])
  for (j in 1:length(infection_classes)) {
    current_merge_incidence_ppy_by_infection <- current_merge_incidence_ppy %>% filter(infection_class == infection_classes[j])
    for (k in 1:length(age_in_days_midpoint)) {
      current_merge_incidence_ppy_by_infection_and_age <- current_merge_incidence_ppy_by_infection %>% filter(age_in_days_midpoint == age_midpoint[k])
      
      if (dim(current_merge_incidence_ppy_by_infection_and_age)[1] == 2){
        rural_df <- (current_merge_incidence_ppy_by_infection_and_age %>% filter(rur_or_urb == "rural"))$value
        urban_df <- (current_merge_incidence_ppy_by_infection_and_age %>% filter(rur_or_urb == "urban"))$value
        
        rural_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "rural"))$pop
        urban_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "urban"))$pop
        
        rural_weight <- rural_pop / (sum(rural_pop, urban_pop))
        urban_weight <- urban_pop / (sum(rural_pop, urban_pop))
        
        weighted_means <- (rural_weight * rural_df) + (urban_weight * urban_df)
        
      } else {
        weighted_means <- current_merge_incidence_ppy_by_infection_and_age$value
      }
      merged_incidence_ppy_df <- rbind(merged_incidence_ppy_df, weighted_means)
    }
    
  }
}

colnames(merged_incidence_ppy_df) = "value"
merged_incidence_ppy_df$area <- rep(area_names, each = length(infection_classes) * length(age_in_days_midpoint))
merged_incidence_ppy_df$infection_class <- rep(rep(infection_classes, each = length(age_in_days_midpoint)), times=length(area_names))
merged_incidence_ppy_df$iso_code <- rep(country_code, times=dim(merged_incidence_ppy_df)[1])
merged_incidence_ppy_df$rur_or_urb <- rep("rural/urban merged", dim(merged_incidence_ppy_df)[1])
merged_incidence_ppy_df$age_group <- rep(age_group_names, times = length(infection_classes) * length(area_names))
merged_incidence_ppy_df$age_in_days_midpoint <- rep(age_in_days_midpoint, times = length(infection_classes) * length(area_names) )
merged_incidence_ppy_df$units <- rep("ppy", times=dim(merged_incidence_ppy_df)[1])
merged_incidence_ppy_df <- merged_incidence_ppy_df %>%
  relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units, value)


########################## INCIDENCE ########################## 

merge_incidence <- bind_rows(incidence_df_rural, incidence_df_urban) %>% 
  group_by(area, timestep, infection_class)

infection_classes <- c("clinical", "severe", "total", "asymptomatic")

# USING WEIGHTED MEAN 
merged_incidence_df <- data.frame()

for (i in 1:length(area_names)) {
  current_merge_incidence <- merge_incidence %>% filter(area==area_names[i])
  for (j in 1:length(infection_classes)) {
    current_merge_incidence_by_infection <- current_merge_incidence %>% filter(infection_class == infection_classes[j])
    for (k in 1:sim_length) {
      current_merge_incidence_by_infection_and_time <- current_merge_incidence_by_infection %>% filter(timestep == k)
      
      if (dim(current_merge_incidence_by_infection_and_time)[1] == 2){
        rural_df <- (current_merge_incidence_by_infection_and_time %>% filter(rur_or_urb == "rural"))[,6:dim(current_merge_incidence_by_infection_and_time)[2]]
        urban_df <- (current_merge_incidence_by_infection_and_time %>% filter(rur_or_urb == "urban"))[,6:dim(current_merge_incidence_by_infection_and_time)[2]]
        
        rural_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "rural"))$pop
        urban_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "urban"))$pop
        
        rural_weight <- rural_pop / (sum(rural_pop, urban_pop))
        urban_weight <- urban_pop / (sum(rural_pop, urban_pop))
        
        weighted_means <- (rural_weight * rural_df) + (urban_weight * urban_df)
        
      } else {
        weighted_means <- current_merge_incidence_by_infection_and_time[,6:dim(current_merge_incidence_by_infection_and_time)[2]]
      }
      merged_incidence_df <- rbind(merged_incidence_df, weighted_means)
    }
    
  }
}

merged_incidence_df$timestep <- rep(1:sim_length, times = length(infection_classes) * length(area_names))
merged_incidence_df$area <- rep(area_names, each = length(infection_classes) * sim_length)
merged_incidence_df$infection_class <- rep(rep(infection_classes, each = sim_length), times=length(area_names))
merged_incidence_df$iso_code <- rep(country_code, times=dim(merged_incidence_df)[1])
merged_incidence_df$rur_or_urb <- rep("rural/urban merged", dim(merged_incidence_df)[1])
merged_incidence_df$units <- rep("ppy", times=dim(merged_incidence_df)[1])
merged_incidence_df <- merged_incidence_df %>%
  relocate(timestep, iso_code, area, rur_or_urb, infection_class, units)


########################## POPULATION ########################## 

merge_population <- bind_rows(population_df_rural, population_df_urban) %>% 
  group_by(area, timestep)

# USING WEIGHTED MEAN 
merged_population_df <- data.frame()

for (i in 1:length(area_names)) {
  current_merge_population <- merge_population %>% filter(area==area_names[i])
  print(i)
  for (k in 1:sim_length) {
      current_merge_population_by_time <- current_merge_population %>% filter(timestep == k)
      
      if (dim(current_merge_population_by_time)[1] == 2){
        rural_df <- (current_merge_population_by_time %>% filter(rur_or_urb == "rural"))[,5:dim(current_merge_population_by_time)[2]]
        urban_df <- (current_merge_population_by_time %>% filter(rur_or_urb == "urban"))[,5:dim(current_merge_population_by_time)[2]]
        
        rural_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "rural"))$pop
        urban_pop <- (full_data$population %>% filter(year == 2023, name_2 == area_names[i], urban_rural == "urban"))$pop
        
        rural_weight <- rural_pop / (sum(rural_pop, urban_pop))
        urban_weight <- urban_pop / (sum(rural_pop, urban_pop))
        
        weighted_means <- (rural_weight * rural_df) + (urban_weight * urban_df)
        
      } else {
        weighted_means <- current_merge_population_by_time[,5:dim(current_merge_population_by_time)[2]]
      }
      merged_population_df <- rbind(merged_population_df, weighted_means)
    }
    
  }


merged_population_df$timestep <- rep(1:sim_length, times = length(area_names))
merged_population_df$area <- rep(area_names, each = sim_length)
merged_population_df$iso_code <- rep(country_code, times=dim(merged_population_df)[1])
merged_population_df$rur_or_urb <- rep("rural/urban merged", dim(merged_population_df)[1])
merged_population_df$units <- rep("ppy", times=dim(merged_population_df)[1])
merged_population_df <- merged_population_df %>%
  relocate(timestep, iso_code, area, rur_or_urb, units)


######################### STORE DATAFRAMES IN CSV FORMAT ########################

# all simulation data
write.csv(merged_sim_output_df, paste0(getwd(), "/simulation_results/sim_output_merged.csv"), row.names=FALSE)

# average incidence per person per year for each site
write.csv(merged_incidence_ppy_df, paste0(getwd(), "/simulation_results/incidence_ppy_df_merged.csv"), row.names=FALSE)

# number of infections in each age group at each timestep, for each site
write.csv(merged_incidence_df, paste0(getwd(), "/simulation_results/incidence_df_merged.csv"), row.names=FALSE)

# population sizes in each age group at each timestep, for each site 
write.csv(merged_population_df, paste0(getwd(), "/simulation_results/population_df_merged.csv"), row.names=FALSE)




