########################## COMMON PARAMETERS/VARIABLES ########################

# areas to loop over within the country 
area_names <- c("Sud", "Nord")#, "Sud_Ouest")#, "Ouest", "Nord_Ouest") 
#area_names <- unique(full_data$sites$name_2)


human_population <- 1000 # population size in model 
min_age_to_model <- 0 # minimum age to model (in days)
max_age_to_model <- 0.05*365 # enter in multiples of 0.5, maximum age to model (in days)

# Set the time span over which to simulate

# NOTE: currently 23 years is the max as we only have ITN/IRS data for this length
year <- 365; years <- 1; sim_length_for_data <- year * years

years_proj_forward <- 1 # number of years to project forward from known data
years_of_simulation <- years + years_proj_forward # simulation length (years) 

sim_length <- (years + years_proj_forward) * year # simulation length (days)

# interval between modeled age groups (1 = days, 7 = weeks etc)
step_length <- 1

# vector of min and max ages to be used in each age bracket
age_min <- seq(min_age_to_model, max_age_to_model, step_length) 
age_max <- seq(min_age_to_model, max_age_to_model, step_length) + step_length

sixmonth_intervals <- c(0, 183, 365, 548, 730, 913, 1095, 1278, 1460, 1643,
                        1825, 2008, 2190, 2373, 2555, 2738, 2920, 3285, 3468,
                        3650, 4015)

# time steps used in model (number of rows in the simulations output data frame)
timesteps<-seq(0,round(sim_length),1)

# vector for the midpoint age for each age bracket
age_in_days <- seq(age_min[1], age_max[length(age_max)], length.out=length(age_min) + 1)
age_in_days_midpoint <- age_in_days[-length(age_in_days)] + diff(age_in_days)/2

no_sixmonth_intervals <- round(max(age_max)/182.5)

# Write vectors for the column names of interest 
age_group_names <- c() # age group columns
age_group_names_sixmonth <- c()
clin_inc_cols <- c() # clinical incidence column names
sev_inc_cols <- c() # severe incidence column names
tot_inc_cols <- c() # total incidence column names
asym_inc_cols <- c() # asymptomatic incidence column names

# fill vectors
for (i in 1:(length(age_min))) {
  age_group_names <- append(age_group_names, paste0("n_age_", as.character(age_min[i]), "_", as.character(age_max[i])))
}

sixmonth_intervals_midpoint <- c()

for (i in 1:(length(sixmonth_intervals) - 1)) {
  
  if (i == (length(sixmonth_intervals) - 1)) {
    age_group_names_sixmonth <- append(age_group_names_sixmonth, paste0("n_age_", as.character(sixmonth_intervals[i]), "_", as.character(sixmonth_intervals[i+1])))
    sixmonth_intervals_midpoint <- append(sixmonth_intervals_midpoint, (sixmonth_intervals[i]+sixmonth_intervals[i+1])/2)
    }
  
  if (i != (length(sixmonth_intervals) - 1)) {
    age_group_names_sixmonth <- append(age_group_names_sixmonth, paste0("n_age_", as.character(sixmonth_intervals[i]), "_", as.character(sixmonth_intervals[i+1] - 1)))
    sixmonth_intervals_midpoint <- append(sixmonth_intervals_midpoint, (sixmonth_intervals[i]+(sixmonth_intervals[i+1] - 1))/2)
    }
  
}

# write column names 
for (i in 1:length(age_min)) {
  clin_inc_cols <- append(clin_inc_cols, paste0("n_inc_clinical_", as.character(age_min[i]), "_", as.character(age_max[i])))
  sev_inc_cols <- append(sev_inc_cols, paste0("n_inc_severe_", as.character(age_min[i]), "_", as.character(age_max[i])))
  tot_inc_cols <- append(tot_inc_cols, paste0("n_inc_", as.character(age_min[i]), "_", as.character(age_max[i])))
  asym_inc_cols <- append(asym_inc_cols, paste0("n_inc_asym_", as.character(age_min[i]), "_", as.character(age_max[i])))
}

############################## RURAL MODEL SIMULATION ########################

rur_or_urb <- "rural"

source("CMR_model.R")

if (rur_or_urb == "rural") {
  
  for (i in 1:length(area_names)){
    
    # Time simulation
    start_time <- Sys.time()
    
    ################### SET UP PARAMETERS AND INTERVENTIONS ##########################
    
    # Current area 
    current_area_name <-area_names[i]
    
    # Starting EIR
    starting_EIR <- (full_data$eir %>%
                       filter(name_2==current_area_name) %>%
                       filter(spp=="pf") %>%
                       filter(urban_rural == rur_or_urb))$eir
    
    
    # If theres no rural data for the current admin-1 area
    if (length(starting_EIR) == 0) {
      print(paste0("No ", rur_or_urb, " data for the area: ", current_area_name))
      next
    }
    
    # If the EIR is 0 
    if (starting_EIR == 0) {
      print(paste0("Starting EIR == 0 for ", rur_or_urb, " area in ", current_area_name, ". Therefore model cannot run."))
      next
    }
    
    
    # SEASONALITY
    seas_data <- full_data$seasonality %>% filter(name_2==current_area_name) 
    
    g0 <- seas_data$g0
    g <- as.double(c(seas_data$g1, seas_data$g2, seas_data$g3))
    h <- as.double(c(seas_data$h1, seas_data$h2, seas_data$h3))
    
    
    # Initialise the age structures to use and seasonality in parameters
    
    simparams <- get_parameters(
      list(
        human_population = human_population,
        age_group_rendering_min_ages = age_min,
        age_group_rendering_max_ages = age_max,
        incidence_rendering_min_ages = age_min,
        incidence_rendering_max_ages = age_max,
        clinical_incidence_rendering_min_ages = age_min,
        clinical_incidence_rendering_max_ages = age_max,
        severe_incidence_rendering_min_ages = age_min,
        severe_incidence_rendering_max_ages = age_max#,
        # model_seasonality = TRUE,
        # g0 = g0,
        # g = g,
        # h = h
      )
    )
    
    # SPECIES
    
    # proportions of An. arabiensis, An. funestus, An. gambiae (in that order)
    
    species_data <- full_data$vectors %>% 
      filter(name_2 == current_area_name)
    
    # for when only two species are present
    # arabiensis
    if (any(species_data$species == "arabiensis") == TRUE) {
      arab_prop <- (species_data %>% filter(species == "arabiensis"))$prop
    } else {
      arab_prop <- 0
    }
    
    # funestus
    if (any(species_data$species == "funestus") == TRUE) {
      fun_prop <- (species_data %>% filter(species == "funestus"))$prop
    } else {
      fun_prop <- 0
    }
    
    # gambiae
    if (any(species_data$species == "gambiae") == TRUE) {
      gamb_prop <- (species_data %>% filter(species == "gambiae"))$prop
    } else {
      gamb_prop <- 0
    }
    
    species_prop <- round(c(arab_prop, fun_prop, gamb_prop),digits=10)
    
    simparams <- set_species(
      simparams,
      species = list(arab_params, fun_params, gamb_params), # must be in same order with prop above
      proportions = species_prop
    )
    
    
    # INTERVENTIONS
    
    # Subset intervention data by area/classification
    intervention_data <- full_data$interventions %>%
      filter(name_2 == current_area_name) %>%
      filter(urban_rural == rur_or_urb)
    
    # subset intervention data by last n years of known data
    intervention_data <- tail(intervention_data, n=years)
     
    # Assume interventions continue at the same coverage and schedule from last known data
    interv_data_proj_forwards <- intervention_data |> expand_interventions(max_year=(tail(intervention_data$year, n=1)+years_proj_forward), group_var="name_2")
    interv_data_proj_forwards <- interv_data_proj_forwards |> fill_extrapolate(group_var = "name_2")
    
    
    # BED NETS
    
    # bed nets assumed to be distributed at the end of each year
    itn_timesteps <- c(1:(sim_length/year)) * year
    itn_timesteps <- itn_timesteps[itn_timesteps <= tail(timesteps, n=1)] # remove values higher than simulation length


    # simparams <- set_bednets(
    #   simparams,
    #   timesteps = itn_timesteps,
    #   coverages = interv_data_proj_forwards$itn_use,  # Each round is distributed to x% of the population randomly
    #   retention = 5 * year, # Nets are kept on average 5 years
    #   dn0 = matrix(interv_data_proj_forwards$dn0, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of death probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   rn = matrix(interv_data_proj_forwards$rn0, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of repelling probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   rnm = matrix(interv_data_proj_forwards$rnm, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of minimum repelling probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   gamman = interv_data_proj_forwards$gamman # Vector of bed net half-lives for each distribution timestep
    # )
    
    
    # IRS
    
    # irs is assumed to be implemented at the beginning of each year
    irs_timesteps <- c(1,(1:(sim_length/year) * year + 1))
    irs_timesteps <- irs_timesteps[irs_timesteps <= tail(timesteps, n=1)] #  remove values higher than simulation length


    # simparams <- set_spraying(
    #   simparams,
    #   timesteps = irs_timesteps,
    #   coverages = interv_data_proj_forwards$irs_cov, # NOTE: coverage is very low so has very little impact on model. tail() to access the last X number of years of data
    #   ls_theta = matrix(interv_data_proj_forwards$ls_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of mortality parameters; nrows=length(timesteps), ncols=length(species)
    #   ls_gamma = matrix(interv_data_proj_forwards$ls_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of mortality parameters per round of IRS and per species
    #   ks_theta = matrix(interv_data_proj_forwards$ks_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of feeding success parameters per round of IRS and per species
    #   ks_gamma = matrix(interv_data_proj_forwards$ks_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of feeding success parameters per round of IRS and per species
    #   ms_theta = matrix(interv_data_proj_forwards$ms_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of deterrence parameters per round of IRS and per species
    #   ms_gamma = matrix(interv_data_proj_forwards$ms_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)) # GUESS: Matrix of deterrence parameters per round of IRS and per species
    # )
    
    
    
    # PMC 
    
    # # PMC timesteps 
    
    # schedule_CMR_8doses<- c(10*7, 14*7, 6*30, 9*30, 12*30, 15*30, 18*30, 24*30)
    # 
    # cov_CMR_8doses<- rep(0.5, sim_length)
    # 
    # testing_timesteps <- (1:sim_length)
    # 
    # #SP_AQ_params<- c(1,0.32,4.3,38.1)
    # simparams <- set_drugs(simparams, list(SP_AQ_params))
    # 
    # 
    # simparams <- set_pmc(
    #   simparams,
    #   drug = 1,
    #   timesteps = testing_timesteps,
    #   coverages = cov_CMR_8doses,
    #   ages = schedule_CMR_8doses
    # )
    
    
    # TREATMENT
    
    # treatment_policy <- read_xlsx("C:/Users/IshanaAtkinson/OneDrive - London School of Hygiene and Tropical Medicine/Documents/malariasimulation work/drug_policies.xlsx", sheet=1)
    # 
    # # most recent value for the country
    # country_policy <- (treatment_policy %>% filter(ISO_CODE == country_code))["2015"]
    # 
    # overall_treatment_cov <- as.numeric(interv_data_proj_forwards$tx_cov)
    # ACT_treatment_cov <- as.numeric(interv_data_proj_forwards$prop_act) * overall_treatment_cov
    # non_ACT_treatment_cov <- overall_treatment_cov - ACT_treatment_cov
    # 
    # treatment_timestep <- seq(from=1, to=years_of_simulation*year, by=year) # change on the new year
    # 
    # sim_params <- set_drugs(simparams, list(AL_params, SP_AQ_params))
    # 
    # if (country_policy == "AL") {
    # sim_params <- set_clinical_treatment(
    #   parameters = sim_params,
    #   drug = 1, # AL (ACT)
    #   timesteps =  treatment_timestep, # Treatment coverage changes
    #   coverages =  ACT_treatment_cov)
    # } else if (country_policy == "AS+AQ") {
    #   sim_params <- set_clinical_treatment(
    #     parameters = sim_params,
    #     drug = 2, # AS-AQ (ACT) ####### INSERT AS-AQ PARAMETERS !!!!!
    #     timesteps =  treatment_timestep, # Treatment coverage changes
    #     coverages =  ACT_treatment_cov)
    # } else {
    #   sim_params <- set_clinical_treatment(
    #     parameters = sim_params,
    #     drug = 1, # assumed AL (ACT) in absence of data 
    #     timesteps =  treatment_timestep, # Treatment coverage changes
    #     coverages =  ACT_treatment_cov)
    # }
    # 
    # sim_params <- set_clinical_treatment(
    #   parameters = sim_params,
    #   drug = 2, # SP-AQ (non-ACT)
    #   timesteps =  treatment_timestep, # Treatment coverage changes
    #   coverages =  non_ACT_treatment_cov)
    
    
    ######################## RUN MODEL ########################################
    
    # Adjusts IBM parameters to match equilibrium parameters
    simparams <- set_equilibrium(simparams, starting_EIR)
    
    # Run simulation 
    test_sim <- run_simulation(sim_length, simparams)
    
    # Assign location data to simulation results 
    test_sim$iso_code <- rep(country_code, dim(test_sim)[1])
    test_sim$rur_or_urb <- rep(rur_or_urb, dim(test_sim)[1])
    test_sim$area <- rep(current_area_name, dim(test_sim)[1])
    test_sim <- test_sim %>%
      relocate(timestep, iso_code, rur_or_urb, area)
    
    
    # Store all simulation results
    sim_output_rural <- rbind(sim_output_rural, test_sim)
    
    
    ##################### CREATE VECTORS OF INTEREST ###########################

    clin_inc_data_mats <- c() # clinical incidence
    sev_inc_data_mats <- c() # severe incidence
    tot_inc_data_mats <- c() # total incidence
    asym_inc_data_mats <- c() # asymptomatic incidence
    pop_size_data_mats <- c() # population size 
    new_clin_inc_per_person <- c() # clinical incidence (per capita)
    new_sev_inc_per_person <- c() # severe incidence (per capita)
    new_tot_inc_per_person <- c() # total incidence (per capita)
    new_asym_inc_per_person <- c() # asymptomatic incidence (per capita)
    
    for (j in 1:length(test_sim$timestep)) {
      
      # Create matrix with data for population size for each age group
      current_pop_size_data <- test_sim[j, age_group_names]
      
      # create matrix with data for new infections for each age group
      current_clin_inc_data <- test_sim[j, clin_inc_cols]
      current_sev_inc_data <- test_sim[j, sev_inc_cols]
      current_tot_inc_data <- test_sim[j, tot_inc_cols]
      current_asym_inc_data <- current_tot_inc_data - current_clin_inc_data
      
      # Populate data frame with data at each timestep 
      clin_inc_data_mats<-rbind(clin_inc_data_mats, current_clin_inc_data)
      sev_inc_data_mats<-rbind(sev_inc_data_mats, current_sev_inc_data)
      tot_inc_data_mats<-rbind(tot_inc_data_mats, current_tot_inc_data)
      asym_inc_data_mats<-rbind(asym_inc_data_mats, current_asym_inc_data)
      pop_size_data_mats<-rbind(pop_size_data_mats, current_pop_size_data)
      
      # New infections per capita at each timestep (NaN appear is pop size = 0)
      new_clin_inc_per_person <- rbind(new_clin_inc_per_person, current_clin_inc_data/current_pop_size_data)
      new_sev_inc_per_person <- rbind(new_sev_inc_per_person, current_sev_inc_data/current_pop_size_data)
      new_tot_inc_per_person <- rbind(new_tot_inc_per_person, current_tot_inc_data/current_pop_size_data)
      new_asym_inc_per_person <- rbind(new_asym_inc_per_person, current_asym_inc_data/current_pop_size_data)
      
    }
    
    # Assign location data to simulation results 
    colnames(clin_inc_data_mats) <- age_group_names
    clin_inc_data_mats$timestep <- 1:sim_length
    clin_inc_data_mats$infection_class <- "clinical"
    clin_inc_data_mats$iso_code <- rep(country_code, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats$area <- rep(current_area_name, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats <- clin_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(sev_inc_data_mats) <- age_group_names
    sev_inc_data_mats$timestep <- 1:sim_length
    sev_inc_data_mats$infection_class <- "severe"
    sev_inc_data_mats$iso_code <- rep(country_code, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats$area <- rep(current_area_name, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats <- sev_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(tot_inc_data_mats) <- age_group_names
    tot_inc_data_mats$timestep <- 1:sim_length
    tot_inc_data_mats$infection_class <- "total"
    tot_inc_data_mats$iso_code <- rep(country_code, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats$area <- rep(current_area_name, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats <- tot_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(asym_inc_data_mats) <- age_group_names
    asym_inc_data_mats$timestep <- 1:sim_length
    asym_inc_data_mats$infection_class <- "asymptomatic"
    asym_inc_data_mats$iso_code <- rep(country_code, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats$area <- rep(current_area_name, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats <- asym_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(pop_size_data_mats) <- age_group_names
    pop_size_data_mats$timestep <- 1:sim_length
    pop_size_data_mats$iso_code <- rep(country_code, dim(pop_size_data_mats)[1])
    pop_size_data_mats$rur_or_urb <- rep(rur_or_urb, dim(pop_size_data_mats)[1])
    pop_size_data_mats$area <- rep(current_area_name, dim(pop_size_data_mats)[1])
    pop_size_data_mats <- pop_size_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area)
    
    
    # store clinical incidence results from each simulation run 
    
    store_clin_inc_data_mats_rural <- data.frame()
    store_sev_inc_data_mats_rural <- data.frame()
    store_tot_inc_data_mats_rural <- data.frame()
    store_asym_inc_data_mats_rural <- data.frame()
    
    store_clin_inc_data_mats_rural <- rbind(store_clin_inc_data_mats_rural, clin_inc_data_mats)
    store_sev_inc_data_mats_rural <- rbind(store_sev_inc_data_mats_rural, sev_inc_data_mats)
    store_tot_inc_data_mats_rural <- rbind(store_tot_inc_data_mats_rural, tot_inc_data_mats)
    store_asym_inc_data_mats_rural <- rbind(store_asym_inc_data_mats_rural, asym_inc_data_mats)
    
    population_df_rural <- rbind(population_df_rural, pop_size_data_mats)
    
    # Create average new clinical infections per person per year for each age group
    
    clin_inc_ppy <- store_clin_inc_ppy_rural <- data.frame()
    sev_inc_ppy <- store_sev_inc_ppy_rural <- data.frame()
    tot_inc_ppy <- store_tot_inc_ppy_rural <- data.frame()
    asym_inc_ppy <- store_asym_inc_ppy_rural <- data.frame()

    for (j in 1:dim(new_clin_inc_per_person)[2]) {
      clin_inc_ppy <- rbind(clin_inc_ppy, sum(new_clin_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_sev_inc_per_person)[2]) {
      sev_inc_ppy <- rbind(sev_inc_ppy, sum(new_sev_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_tot_inc_per_person)[2]) {
      tot_inc_ppy <- rbind(tot_inc_ppy, sum(new_tot_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_asym_inc_per_person)[2]) {
      asym_inc_ppy <- rbind(asym_inc_ppy, sum(new_asym_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    
    # Assign location and age group data to simulation results
    
    colnames(clin_inc_ppy) = "value"
    clin_inc_ppy$area <- current_area_name
    clin_inc_ppy$infection_class <- "clinical"
    clin_inc_ppy$iso_code <- country_code
    clin_inc_ppy$rur_or_urb <- rur_or_urb
    clin_inc_ppy$age_group <- age_group_names 
    clin_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    clin_inc_ppy$units <- "ppy" 
    clin_inc_ppy <- clin_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(sev_inc_ppy) = "value"
    sev_inc_ppy$area <- current_area_name
    sev_inc_ppy$infection_class <- "severe"
    sev_inc_ppy$iso_code <- country_code
    sev_inc_ppy$rur_or_urb <- rur_or_urb
    sev_inc_ppy$age_group <- age_group_names
    sev_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    sev_inc_ppy$units <- "ppy" 
    sev_inc_ppy <- sev_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(tot_inc_ppy) = "value"
    tot_inc_ppy$area <- current_area_name
    tot_inc_ppy$infection_class <- "total"
    tot_inc_ppy$iso_code <- country_code
    tot_inc_ppy$rur_or_urb <- rur_or_urb
    tot_inc_ppy$age_group <- age_group_names
    tot_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    tot_inc_ppy$units <- "ppy" 
    tot_inc_ppy <- tot_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(asym_inc_ppy) = "value"
    asym_inc_ppy$area <- current_area_name    
    asym_inc_ppy$infection_class <- "asymptomatic"
    asym_inc_ppy$iso_code <- country_code
    asym_inc_ppy$rur_or_urb <- rur_or_urb
    asym_inc_ppy$age_group <- age_group_names
    asym_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    asym_inc_ppy$units <- "ppy" 
    asym_inc_ppy <- asym_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    
    store_clin_inc_ppy_rural <- rbind(store_clin_inc_ppy_rural, clin_inc_ppy)
    store_sev_inc_ppy_rural <- rbind(store_sev_inc_ppy_rural, sev_inc_ppy)
    store_tot_inc_ppy_rural <- rbind(store_tot_inc_ppy_rural, tot_inc_ppy)
    store_asym_inc_ppy_rural <- rbind(store_asym_inc_ppy_rural, asym_inc_ppy)
    
    # store total mosquito population at each timestep
    total_M <- c()
    
    for (i in 1:length(test_sim$timestep)) {
      total_M <- c(total_M, sum(test_sim[i, c("total_M_gamb", "total_M_arab", "total_M_fun")]))
    }
    
    store_total_M_rural <- rbind(store_total_M_rural, total_M)
    
    
    # merge incidence dataframes together 
    
    current_incidence_ppy_df_rural <- rbind(store_clin_inc_ppy_rural, store_sev_inc_ppy_rural,
                                    store_tot_inc_ppy_rural, store_asym_inc_ppy_rural)
    
    current_incidence_df_rural <- rbind(store_clin_inc_data_mats_rural, store_sev_inc_data_mats_rural,
                                store_tot_inc_data_mats_rural, store_asym_inc_data_mats_rural)
    
    incidence_ppy_df_rural <- rbind(incidence_ppy_df_rural, current_incidence_ppy_df_rural)
    incidence_df_rural <- rbind(incidence_df_rural, current_incidence_df_rural)
    
    # End timer
    end_time <- Sys.time()
    
    # Time taken for simulation
    time_elapsed <- end_time - start_time
    
    # PRINT UPDATES
    print(paste0("Successful run for the ", rur_or_urb, " region in: ", current_area_name, ". ",
                 "Time taken to run model loop: ", time_elapsed))
    
    
  }

  
  
######################### store simulation data in CSV format ########################

  # all simulation data
  write.csv(sim_output_rural, paste0(getwd(), "/simulation_results/sim_output_rural.csv"), row.names=FALSE)
  
  # average clinical incidence per person per year for each site
  write.csv(incidence_ppy_df_rural, paste0(getwd(), "/simulation_results/incidence_ppy_df_rural.csv"), row.names=FALSE)
  
  # number of clinical infections in each age group at each timestep, for each site
  write.csv(incidence_df_rural, paste0(getwd(), "/simulation_results/incidence_df_rural.csv"), row.names=FALSE)
  
  # population sizes in each age group at each timestep, for each site 
  write.csv(population_df_rural, paste0(getwd(), "/simulation_results/population_df_rural.csv"), row.names=FALSE)
  

  
}


############################# URBAN SIMULATION ################################

rur_or_urb <- "urban"

source("CMR_model.R")

if (rur_or_urb == "urban") {
  
  for (i in 1:length(area_names)){
    
    # Time simulation
    start_time <- Sys.time()
    
    ################### SET UP PARAMETERS AND INTERVENTIONS ##########################
    
    # Current area 
    current_area_name <-area_names[i]
    
    # Starting EIR
    starting_EIR <- (full_data$eir %>%
                       filter(name_2==current_area_name) %>%
                       filter(spp=="pf") %>%
                       filter(urban_rural == rur_or_urb))$eir
    
    
    # If theres no urban data for the current admin-1 area
    if (length(starting_EIR) == 0) {
      print(paste0("No ", rur_or_urb, " data for the area: ", current_area_name))
      next
    }
    
    # If the EIR is 0 
    if (starting_EIR == 0) {
      print(paste0("Starting EIR == 0 for ", rur_or_urb, " area in ", current_area_name, ". Therefore model cannot run."))
      next
    }
    
    
    # SEASONALITY
    seas_data <- full_data$seasonality %>% filter(name_2==current_area_name) 
    
    g0 <- seas_data$g0
    g <- as.double(c(seas_data$g1, seas_data$g2, seas_data$g3))
    h <- as.double(c(seas_data$h1, seas_data$h2, seas_data$h3))
    
    
    # Initialise the age structures to use and seasonality in parameters
    
    simparams <- get_parameters(
      list(
        human_population = human_population,
        age_group_rendering_min_ages = age_min,
        age_group_rendering_max_ages = age_max,
        incidence_rendering_min_ages = age_min,
        incidence_rendering_max_ages = age_max,
        clinical_incidence_rendering_min_ages = age_min,
        clinical_incidence_rendering_max_ages = age_max,
        severe_incidence_rendering_min_ages = age_min,
        severe_incidence_rendering_max_ages = age_max#,
        # model_seasonality = TRUE,
        # g0 = g0,
        # g = g,
        # h = h
      )
    )
    
    # SPECIES
    
    # proportions of An. arabiensis, An. funestus, An. gambiae (in that order)
    
    species_data <- full_data$vectors %>% 
      filter(name_2 == current_area_name)
    
    # for when only two species are present
    # arabiensis
    if (any(species_data$species == "arabiensis") == TRUE) {
      arab_prop <- (species_data %>% filter(species == "arabiensis"))$prop
    } else {
      arab_prop <- 0
    }
    
    # funestus
    if (any(species_data$species == "funestus") == TRUE) {
      fun_prop <- (species_data %>% filter(species == "funestus"))$prop
    } else {
      fun_prop <- 0
    }
    
    # gambiae
    if (any(species_data$species == "gambiae") == TRUE) {
      gamb_prop <- (species_data %>% filter(species == "gambiae"))$prop
    } else {
      gamb_prop <- 0
    }
    
    species_prop <- round(c(arab_prop, fun_prop, gamb_prop),digits=10)
    
    simparams <- set_species(
      simparams,
      species = list(arab_params, fun_params, gamb_params), # must be in same order with prop above
      proportions = species_prop
    )
    
    
    # INTERVENTIONS
    
    # Subset intervention data by area/classification
    intervention_data <- full_data$interventions %>%
      filter(name_2 == current_area_name) %>%
      filter(urban_rural == rur_or_urb)
    
    # subset intervention data by last n years of known data
    intervention_data <- tail(intervention_data, n=years)
    
    # Assume interventions continue at the same coverage and schedule from last known data
    interv_data_proj_forwards <- intervention_data |> expand_interventions(max_year=(tail(intervention_data$year, n=1)+years_proj_forward), group_var="name_2")
    interv_data_proj_forwards <- interv_data_proj_forwards |> fill_extrapolate(group_var = "name_2")
    
    
    # BED NETS
    
    # bed nets assumed to be distributed at the end of each year
    itn_timesteps <- c(1:(sim_length/year)) * year
    itn_timesteps <- itn_timesteps[itn_timesteps <= tail(timesteps, n=1)] # remove values higher than simulation length
    
    
    # simparams <- set_bednets(
    #   simparams,
    #   timesteps = itn_timesteps,
    #   coverages = interv_data_proj_forwards$itn_use,  # Each round is distributed to x% of the population randomly
    #   retention = 5 * year, # Nets are kept on average 5 years
    #   dn0 = matrix(interv_data_proj_forwards$dn0, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of death probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   rn = matrix(interv_data_proj_forwards$rn0, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of repelling probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   rnm = matrix(interv_data_proj_forwards$rnm, nrow = length(itn_timesteps), ncol = length(simparams$species)), # Matrix of minimum repelling probabilities for each mosquito species (here 3 species only so 3 columns) for last X years
    #   gamman = interv_data_proj_forwards$gamman # Vector of bed net half-lives for each distribution timestep
    # )
    
    
    # IRS
    
    # irs is assumed to be implemented at the beginning of each year
    irs_timesteps <- c(1,(1:(sim_length/year) * year + 1))
    irs_timesteps <- irs_timesteps[irs_timesteps <= tail(timesteps, n=1)] #  remove values higher than simulation length
    
    
    # simparams <- set_spraying(
    #   simparams,
    #   timesteps = irs_timesteps,
    #   coverages = interv_data_proj_forwards$irs_cov, # NOTE: coverage is very low so has very little impact on model. tail() to access the last X number of years of data
    #   ls_theta = matrix(interv_data_proj_forwards$ls_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of mortality parameters; nrows=length(timesteps), ncols=length(species)
    #   ls_gamma = matrix(interv_data_proj_forwards$ls_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of mortality parameters per round of IRS and per species
    #   ks_theta = matrix(interv_data_proj_forwards$ks_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of feeding success parameters per round of IRS and per species
    #   ks_gamma = matrix(interv_data_proj_forwards$ks_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of feeding success parameters per round of IRS and per species
    #   ms_theta = matrix(interv_data_proj_forwards$ms_theta, nrow = length(irs_timesteps), ncol = length(simparams$species)), # GUESS: Matrix of deterrence parameters per round of IRS and per species
    #   ms_gamma = matrix(interv_data_proj_forwards$ms_gamma, nrow = length(irs_timesteps), ncol = length(simparams$species)) # GUESS: Matrix of deterrence parameters per round of IRS and per species
    # )
    
    
    
    # PMC 
    
    # # PMC timesteps 
    # schedule_CMR_8doses<- c(10*7, 14*7, 6*30, 9*30, 12*30, 15*30, 18*30, 24*30)
    # schedule_CMR_5doses<- c(10*7, 14*7, 6*30, 9*30, 15*30)
    # 
    # # PMC coverage
    # cov_CMR_8doses<- c(0.82, 0.73, 0.26, 0.52, 0.15, 0.24 ,0.09, 0.03)
    # cov_CMR_5doses<- c(0.82, 0.73, 0.26, 0.52, 0.24)
    # 
    # testing_timesteps <- c(500,500,500,500,500,500,500,500)
    # 
    # 
    # simparams <- set_drugs(simparams, list(SP_AQ_params))
    # 
    # simparams <- set_pmc(
    #   simparams,
    #   drug = 1,
    #   timesteps = testing_timesteps,
    #   coverages = cov_CMR_8doses,
    #   ages = schedule_CMR_8doses
    # )
    
    
    # TREATMENT
    
    treatment_cov <- as.numeric(interv_data_proj_forwards$tx_cov)
    treatment_timestep <- seq(from=1, to=years_of_simulation*year, by=year) # change on the new year
    
    # sim_params <- set_drugs(simparams, list(SP_AQ_params))
    # 
    # sim_params <- set_clinical_treatment(
    #   parameters = sim_params,
    #   drug = 1, # SP-AQ
    #   timesteps =  treatment_timestep, # Treatment coverage changes
    #   coverages =  treatment_cov)
    
    
    
    ######################## RUN MODEL ########################################
    
    # Adjusts IBM parameters to match equilibrium parameters
    simparams <- set_equilibrium(simparams, starting_EIR)
    
    # Run simulation 
    test_sim <- run_simulation(sim_length, simparams)
    
    # Assign location data to simulation results 
    test_sim$iso_code <- rep(country_code, dim(test_sim)[1])
    test_sim$rur_or_urb <- rep(rur_or_urb, dim(test_sim)[1])
    test_sim$area <- rep(current_area_name, dim(test_sim)[1])
    test_sim <- test_sim %>%
      relocate(timestep, iso_code, rur_or_urb, area)
    
    
    # Store all simulation results
    sim_output_urban <- rbind(sim_output_urban, test_sim)
    
    
    ##################### CREATE VECTORS OF INTEREST ###########################
    
    clin_inc_data_mats <- c() # clinical incidence
    sev_inc_data_mats <- c() # severe incidence
    tot_inc_data_mats <- c() # total incidence
    asym_inc_data_mats <- c() # asymptomatic incidence
    pop_size_data_mats <- c() # population size 
    new_clin_inc_per_person <- c() # clinical incidence (per capita)
    new_sev_inc_per_person <- c() # severe incidence (per capita)
    new_tot_inc_per_person <- c() # total incidence (per capita)
    new_asym_inc_per_person <- c() # asymptomatic incidence (per capita)
    
    for (j in 1:length(test_sim$timestep)) {
      
      # Create matrix with data for population size for each age group
      current_pop_size_data <- test_sim[j, age_group_names]
      
      # create matrix with data for new infections for each age group
      current_clin_inc_data <- test_sim[j, clin_inc_cols]
      current_sev_inc_data <- test_sim[j, sev_inc_cols]
      current_tot_inc_data <- test_sim[j, tot_inc_cols]
      current_asym_inc_data <- current_tot_inc_data - current_clin_inc_data
      
      # Populate data frame with data at each timestep 
      clin_inc_data_mats<-rbind(clin_inc_data_mats, current_clin_inc_data)
      sev_inc_data_mats<-rbind(sev_inc_data_mats, current_sev_inc_data)
      tot_inc_data_mats<-rbind(tot_inc_data_mats, current_tot_inc_data)
      asym_inc_data_mats<-rbind(asym_inc_data_mats, current_asym_inc_data)
      pop_size_data_mats<-rbind(pop_size_data_mats, current_pop_size_data)
      
      # New infections per capita at each timestep (NaN appear is pop size = 0)
      new_clin_inc_per_person <- rbind(new_clin_inc_per_person, current_clin_inc_data/current_pop_size_data)
      new_sev_inc_per_person <- rbind(new_sev_inc_per_person, current_sev_inc_data/current_pop_size_data)
      new_tot_inc_per_person <- rbind(new_tot_inc_per_person, current_tot_inc_data/current_pop_size_data)
      new_asym_inc_per_person <- rbind(new_asym_inc_per_person, current_asym_inc_data/current_pop_size_data)
      
    }
    
    # Assign location data to simulation results 
    colnames(clin_inc_data_mats) <- age_group_names
    clin_inc_data_mats$timestep <- 1:sim_length
    clin_inc_data_mats$infection_class <- "clinical"
    clin_inc_data_mats$iso_code <- rep(country_code, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats$area <- rep(current_area_name, dim(clin_inc_data_mats)[1])
    clin_inc_data_mats <- clin_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(sev_inc_data_mats) <- age_group_names
    sev_inc_data_mats$timestep <- 1:sim_length
    sev_inc_data_mats$infection_class <- "severe"
    sev_inc_data_mats$iso_code <- rep(country_code, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats$area <- rep(current_area_name, dim(sev_inc_data_mats)[1])
    sev_inc_data_mats <- sev_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(tot_inc_data_mats) <- age_group_names
    tot_inc_data_mats$timestep <- 1:sim_length
    tot_inc_data_mats$infection_class <- "total"
    tot_inc_data_mats$iso_code <- rep(country_code, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats$area <- rep(current_area_name, dim(tot_inc_data_mats)[1])
    tot_inc_data_mats <- tot_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(asym_inc_data_mats) <- age_group_names
    asym_inc_data_mats$timestep <- 1:sim_length
    asym_inc_data_mats$infection_class <- "asymptomatic"
    asym_inc_data_mats$iso_code <- rep(country_code, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats$rur_or_urb <- rep(rur_or_urb, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats$area <- rep(current_area_name, dim(asym_inc_data_mats)[1])
    asym_inc_data_mats <- asym_inc_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area, infection_class)
    
    colnames(pop_size_data_mats) <- age_group_names
    pop_size_data_mats$timestep <- 1:sim_length
    pop_size_data_mats$iso_code <- rep(country_code, dim(pop_size_data_mats)[1])
    pop_size_data_mats$rur_or_urb <- rep(rur_or_urb, dim(pop_size_data_mats)[1])
    pop_size_data_mats$area <- rep(current_area_name, dim(pop_size_data_mats)[1])
    pop_size_data_mats <- pop_size_data_mats %>%
      relocate(timestep, iso_code, rur_or_urb, area)
    
    
    # store clinical incidence results from each simulation run 
    
    store_clin_inc_data_mats_urban <- data.frame()
    store_sev_inc_data_mats_urban <- data.frame()
    store_tot_inc_data_mats_urban <- data.frame()
    store_asym_inc_data_mats_urban <- data.frame()
    
    store_clin_inc_data_mats_urban <- rbind(store_clin_inc_data_mats_urban, clin_inc_data_mats)
    store_sev_inc_data_mats_urban <- rbind(store_sev_inc_data_mats_urban, sev_inc_data_mats)
    store_tot_inc_data_mats_urban <- rbind(store_tot_inc_data_mats_urban, tot_inc_data_mats)
    store_asym_inc_data_mats_urban <- rbind(store_asym_inc_data_mats_urban, asym_inc_data_mats)
    
    population_df_urban <- rbind(population_df_urban, pop_size_data_mats)
    
    # Create average new clinical infections per person per year for each age group
    
    clin_inc_ppy <- store_clin_inc_ppy_urban <- data.frame()
    sev_inc_ppy <- store_sev_inc_ppy_urban <- data.frame()
    tot_inc_ppy <- store_tot_inc_ppy_urban <- data.frame()
    asym_inc_ppy <- store_asym_inc_ppy_urban <- data.frame()
    
    for (j in 1:dim(new_clin_inc_per_person)[2]) {
      clin_inc_ppy <- rbind(clin_inc_ppy, sum(new_clin_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_sev_inc_per_person)[2]) {
      sev_inc_ppy <- rbind(sev_inc_ppy, sum(new_sev_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_tot_inc_per_person)[2]) {
      tot_inc_ppy <- rbind(tot_inc_ppy, sum(new_tot_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    for (j in 1:dim(new_asym_inc_per_person)[2]) {
      asym_inc_ppy <- rbind(asym_inc_ppy, sum(new_asym_inc_per_person[,j], na.rm=TRUE)/years_of_simulation)
    }
    
    
    # Assign location and age group data to simulation results
    
    colnames(clin_inc_ppy) = "value"
    clin_inc_ppy$area <- current_area_name
    clin_inc_ppy$infection_class <- "clinical"
    clin_inc_ppy$iso_code <- country_code
    clin_inc_ppy$rur_or_urb <- rur_or_urb
    clin_inc_ppy$age_group <- age_group_names 
    clin_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    clin_inc_ppy$units <- "ppy" 
    clin_inc_ppy <- clin_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(sev_inc_ppy) = "value"
    sev_inc_ppy$area <- current_area_name
    sev_inc_ppy$infection_class <- "severe"
    sev_inc_ppy$iso_code <- country_code
    sev_inc_ppy$rur_or_urb <- rur_or_urb
    sev_inc_ppy$age_group <- age_group_names
    sev_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    sev_inc_ppy$units <- "ppy" 
    sev_inc_ppy <- sev_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(tot_inc_ppy) = "value"
    tot_inc_ppy$area <- current_area_name
    tot_inc_ppy$infection_class <- "total"
    tot_inc_ppy$iso_code <- country_code
    tot_inc_ppy$rur_or_urb <- rur_or_urb
    tot_inc_ppy$age_group <- age_group_names
    tot_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    tot_inc_ppy$units <- "ppy" 
    tot_inc_ppy <- tot_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    colnames(asym_inc_ppy) = "value"
    asym_inc_ppy$area <- current_area_name    
    asym_inc_ppy$infection_class <- "asymptomatic"
    asym_inc_ppy$iso_code <- country_code
    asym_inc_ppy$rur_or_urb <- rur_or_urb
    asym_inc_ppy$age_group <- age_group_names
    asym_inc_ppy$age_in_days_midpoint <- age_in_days_midpoint
    asym_inc_ppy$units <- "ppy" 
    asym_inc_ppy <- asym_inc_ppy %>%
      relocate(iso_code, area, rur_or_urb, age_group, age_in_days_midpoint, infection_class, units)
    
    
    store_clin_inc_ppy_urban <- rbind(store_clin_inc_ppy_urban, clin_inc_ppy)
    store_sev_inc_ppy_urban <- rbind(store_sev_inc_ppy_urban, sev_inc_ppy)
    store_tot_inc_ppy_urban <- rbind(store_tot_inc_ppy_urban, tot_inc_ppy)
    store_asym_inc_ppy_urban <- rbind(store_asym_inc_ppy_urban, asym_inc_ppy)
    
    # store total mosquito population at each timestep
    total_M <- c()
    
    for (i in 1:length(test_sim$timestep)) {
      total_M <- c(total_M, sum(test_sim[i, c("total_M_gamb", "total_M_arab", "total_M_fun")]))
    }
    
    store_total_M_urban <- rbind(store_total_M_urban, total_M)
    
    
    # merge incidence dataframes together 
    
    current_incidence_ppy_df_urban <- rbind(store_clin_inc_ppy_urban, store_sev_inc_ppy_urban,
                                            store_tot_inc_ppy_urban, store_asym_inc_ppy_urban)
    
    current_incidence_df_urban <- rbind(store_clin_inc_data_mats_urban, store_sev_inc_data_mats_urban,
                                        store_tot_inc_data_mats_urban, store_asym_inc_data_mats_urban)
    
    incidence_ppy_df_urban <- rbind(incidence_ppy_df_urban, current_incidence_ppy_df_urban)
    incidence_df_urban <- rbind(incidence_df_urban, current_incidence_df_urban)
    
    # End timer
    end_time <- Sys.time()
    
    # Time taken for simulation
    time_elapsed <- end_time - start_time
    
    # PRINT UPDATES
    print(paste0("Successful run for the ", rur_or_urb, " region in: ", current_area_name, ". ",
                 "Time taken to run model loop: ", time_elapsed))
    
    
  }
  
  
  
  ######################### store simulation data in CSV format ########################
  
  # all simulation data
  write.csv(sim_output_urban, paste0(getwd(), "/simulation_results/sim_output_urban.csv"), row.names=FALSE)
  
  # average incidence per person per year for each site
  write.csv(incidence_ppy_df_urban, paste0(getwd(), "/simulation_results/incidence_ppy_df_urban.csv"), row.names=FALSE)
  
  # number of infections in each age group at each timestep, for each site
  write.csv(incidence_df_urban, paste0(getwd(), "/simulation_results/incidence_df_urban.csv"), row.names=FALSE)
  
  # population sizes in each age group at each timestep, for each site 
  write.csv(population_df_urban, paste0(getwd(), "/simulation_results/population_df_urban.csv"), row.names=FALSE)
  
}


