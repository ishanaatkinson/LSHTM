# source("CMR_model.R")
#source("running_model_CMR.R")
#source("CMR_merge_results.R")
#source("CMR_PMC_post_processing.R")

############################### READ IN PMC IMPACT RESULTS ####################

#sim_output <- read.csv(paste0(getwd(), "/simulation_results/sim_output_merged.csv"))
incidence_ppy_df <- read.csv(paste0(getwd(), "/simulation_results/incidence_ppy_df_merged.csv"))
#incidence_df <- read.csv(paste0(getwd(), "/simulation_results/incidence_df_merged.csv"))
population_df <- read.csv(paste0(getwd(), "/simulation_results/population_df_merged.csv"))

PMC_impact_ppy_8dose<-read.csv(paste0(getwd(), "/simulation_results/PMC_impact_ppy_8dose.csv"))
PMC_impact_ppy_5dose<-read.csv(paste0(getwd(), "/simulation_results/PMC_impact_ppy_5dose.csv"))

################################ CALCULATIONS #################################

# Cases averted 
# Using model results on clinical incidence, foresite population
# estimates and age distributions from the malariasimulation model 

cases_no_PMC<- data.frame(age_group = rep(age_group_names,length(area_names)), age_in_days_midpoint = rep(age_in_days_midpoint,length(area_names)),
                          NAME_2 = rep(area_names, each=length(age_in_days_midpoint)),
                          units=rep("raw cases", length(area_names)*length(age_in_days_midpoint)))

annual_cases_no_PMC<- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                NAME_2 = area_names,
                                units=rep("raw cases", length(area_names)))

cases_8dose<- data.frame(age_group = rep(age_group_names,length(area_names)), age_in_days_midpoint = rep(age_in_days_midpoint,length(area_names)),
                         NAME_2 = rep(area_names, each=length(age_in_days_midpoint)),
                         units=rep("raw cases", length(area_names)*length(age_in_days_midpoint)))

annual_cases_8dose<- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                        NAME_2 = area_names,
                                        units=rep("raw cases", length(area_names)))

cases_5dose<- data.frame(age_group = rep(age_group_names,length(area_names)), age_in_days_midpoint = rep(age_in_days_midpoint,length(area_names)),
                         NAME_2 = rep(area_names, each=length(age_in_days_midpoint)),
                         units=rep("raw cases", length(area_names)*length(age_in_days_midpoint)))

annual_cases_5dose<- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                        NAME_2 = area_names,
                                        units=rep("raw cases", length(area_names)))


cases_averted_8dose<- data.frame(age_group = rep(age_group_names,length(area_names)), age_in_days_midpoint = rep(age_in_days_midpoint,length(area_names)),
                                 NAME_2 = rep(area_names, each=length(age_in_days_midpoint)),
                                 units=rep("raw cases", length(area_names)*length(age_in_days_midpoint)))

annual_cases_averted_8dose<- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                 NAME_2 = area_names,
                                 units=rep("raw cases", length(area_names)))

cases_averted_5dose<- data.frame(age_group = rep(age_group_names,length(area_names)), age_in_days_midpoint = rep(age_in_days_midpoint,length(area_names)),
                                 NAME_2 = rep(area_names, each=length(age_in_days_midpoint)),
                                 units=rep("raw cases", length(area_names)*length(age_in_days_midpoint)))

annual_cases_averted_5dose<- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                        NAME_2 = area_names,
                                        units=rep("raw cases", length(area_names)))

cases_reduction_8dose <- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                    units=rep("average reduction (%)", length(area_names)), NAME_2 = area_names)

cases_reduction_5dose <- data.frame(age_group = paste0("n_age_0_", max(age_max)),
                                    units=rep("average reduction (%)", length(area_names)), NAME_2 = area_names)


current_cases_no_PMC_clin <- current_cases_no_PMC_sev <- current_cases_no_PMC_tot <- current_cases_no_PMC_asym <- c()
current_cases_averted_8dose_clin <- current_cases_averted_8dose_sev <- current_cases_averted_8dose_tot <- current_cases_averted_8dose_asym <- c()
current_cases_averted_5dose_clin <- current_cases_averted_5dose_sev <- current_cases_averted_5dose_tot <- current_cases_averted_5dose_asym <- c()
current_cases_8dose_clin <- current_cases_8dose_sev <- current_cases_8dose_tot <- current_cases_8dose_asym <- c()
current_cases_5dose_clin <- current_cases_5dose_sev <- current_cases_5dose_tot <- current_cases_5dose_asym <- c()
current_cases_reduction_8dose_clin <- current_cases_reduction_8dose_sev <- current_cases_reduction_8dose_tot <- current_cases_reduction_8dose_asym <- c()
current_cases_reduction_5dose_clin <- current_cases_reduction_5dose_sev <- current_cases_reduction_5dose_tot <- current_cases_reduction_5dose_asym <- c()

annual_current_cases_no_PMC_clin <- annual_current_cases_no_PMC_sev <- annual_current_cases_no_PMC_tot <- annual_current_cases_no_PMC_asym <- c()
annual_current_cases_8dose_clin <- annual_current_cases_8dose_sev <- annual_current_cases_8dose_tot <- annual_current_cases_8dose_asym <- c()
annual_current_cases_5dose_clin <- annual_current_cases_5dose_sev <- annual_current_cases_5dose_tot <- annual_current_cases_5dose_asym <- c()
annual_current_cases_averted_8dose_clin <- annual_current_cases_averted_8dose_sev <- annual_current_cases_averted_8dose_tot <- annual_current_cases_averted_8dose_asym <- c()
annual_current_cases_averted_5dose_clin <- annual_current_cases_averted_5dose_sev <- annual_current_cases_averted_5dose_tot <- annual_current_cases_averted_5dose_asym <- c()


for (i in 1:length(area_names)){
  
  # extracting average age proportions from the model age distributions 
  age_distribution <- rep(colMeans((population_df %>% 
                                      filter(area == area_names[i]))[, 6:(dim(population_df)[2])] / human_population), 4)
  
  # sum of rural and urban areas as weighted mean is used 
  area_population <- sum((full_data$population %>% 
                        filter(name_2 == area_names[i], year==2023))$pop)
  
  # average number of cases in each age group per year (no PMC)
  current_cases_no_PMC <- (incidence_ppy_df %>% filter(area==area_names[i]))$value * age_distribution * area_population
  
  current_cases_no_PMC_clin <- c(current_cases_no_PMC_clin, current_cases_no_PMC[1:max(age_max)])
  current_cases_no_PMC_sev <- c(current_cases_no_PMC_sev, current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))])
  current_cases_no_PMC_tot <- c(current_cases_no_PMC_tot, current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))])
  current_cases_no_PMC_asym <- c(current_cases_no_PMC_asym, current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))])
  
  annual_current_cases_no_PMC_clin <- c(annual_current_cases_no_PMC_clin, sum(current_cases_no_PMC[1:max(age_max)]))
  annual_current_cases_no_PMC_sev <- c(annual_current_cases_no_PMC_sev, sum(current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))]))
  annual_current_cases_no_PMC_tot <- c(annual_current_cases_no_PMC_tot, sum(current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))]))
  annual_current_cases_no_PMC_asym <- c(annual_current_cases_no_PMC_asym, sum(current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))]))
  
  # average number of cases in each age group per person per year (with PMC)
  cases_averted_8dose_ppy <- (incidence_ppy_df %>% filter(area==area_names[i]))$value - (PMC_impact_ppy_8dose %>% filter(area==area_names[i]))$value
  cases_averted_5dose_ppy <- (incidence_ppy_df %>% filter(area==area_names[i]))$value - (PMC_impact_ppy_5dose %>% filter(area==area_names[i]))$value
  
  # difference in cases per person per year in each age group X proportion of population in each age group * total population in the area
  current_cases_averted_8dose <- cases_averted_8dose_ppy * age_distribution * area_population
  current_cases_averted_5dose <- cases_averted_5dose_ppy * age_distribution * area_population
  
  current_cases_averted_8dose_clin <- c(current_cases_averted_8dose_clin, current_cases_averted_8dose[1:max(age_max)])
  current_cases_averted_8dose_sev <- c(current_cases_averted_8dose_sev, current_cases_averted_8dose[(1+max(age_max)):(2*max(age_max))])
  current_cases_averted_8dose_tot <- c(current_cases_averted_8dose_tot, current_cases_averted_8dose[(1+2*max(age_max)):(3*max(age_max))])
  current_cases_averted_8dose_asym <- c(current_cases_averted_8dose_asym, current_cases_averted_8dose[(1+3*max(age_max)):(4*max(age_max))])
  
  current_cases_averted_5dose_clin <- c(current_cases_averted_5dose_clin, current_cases_averted_5dose[1:max(age_max)])
  current_cases_averted_5dose_sev <- c(current_cases_averted_5dose_sev, current_cases_averted_5dose[(1+max(age_max)):(2*max(age_max))])
  current_cases_averted_5dose_tot <- c(current_cases_averted_5dose_tot, current_cases_averted_5dose[(1+2*max(age_max)):(3*max(age_max))])
  current_cases_averted_5dose_asym <- c(current_cases_averted_5dose_asym, current_cases_averted_5dose[(1+3*max(age_max)):(4*max(age_max))])
  
  annual_current_cases_averted_8dose_clin <- c(annual_current_cases_averted_8dose_clin, sum(current_cases_averted_8dose[1:max(age_max)]))
  annual_current_cases_averted_8dose_sev <- c(annual_current_cases_averted_8dose_sev, sum(current_cases_averted_8dose[(1+max(age_max)):(2*max(age_max))]))
  annual_current_cases_averted_8dose_tot <- c(annual_current_cases_averted_8dose_tot, sum(current_cases_averted_8dose[(1+2*max(age_max)):(3*max(age_max))]))
  annual_current_cases_averted_8dose_asym <- c(annual_current_cases_averted_8dose_asym, sum(current_cases_averted_8dose[(1+3*max(age_max)):(4*max(age_max))]))
  
  annual_current_cases_averted_5dose_clin <- c(annual_current_cases_averted_5dose_clin, sum(current_cases_averted_5dose[1:max(age_max)]))
  annual_current_cases_averted_5dose_sev <- c(annual_current_cases_averted_5dose_sev, sum(current_cases_averted_5dose[(1+max(age_max)):(2*max(age_max))]))
  annual_current_cases_averted_5dose_tot <- c(annual_current_cases_averted_5dose_tot, sum(current_cases_averted_5dose[(1+2*max(age_max)):(3*max(age_max))]))
  annual_current_cases_averted_5dose_asym <- c(annual_current_cases_averted_5dose_asym, sum(current_cases_averted_5dose[(1+3*max(age_max)):(4*max(age_max))]))
  
  current_cases_8dose <- (PMC_impact_ppy_8dose %>% filter(area==area_names[i]))$value * age_distribution * area_population
  current_cases_5dose <- (PMC_impact_ppy_5dose %>% filter(area==area_names[i]))$value * age_distribution * area_population
  
  current_cases_8dose_clin <- c(current_cases_8dose_clin, current_cases_8dose[1:max(age_max)])
  current_cases_8dose_sev <- c(current_cases_8dose_sev, current_cases_8dose[(1+max(age_max)):(2*max(age_max))])
  current_cases_8dose_tot <- c(current_cases_8dose_tot, current_cases_8dose[(1+2*max(age_max)):(3*max(age_max))])
  current_cases_8dose_asym <- c(current_cases_8dose_asym, current_cases_8dose[(1+3*max(age_max)):(4*max(age_max))])
  
  current_cases_5dose_clin <- c(current_cases_5dose_clin, current_cases_5dose[1:max(age_max)])
  current_cases_5dose_sev <- c(current_cases_5dose_sev, current_cases_5dose[(1+max(age_max)):(2*max(age_max))])
  current_cases_5dose_tot <- c(current_cases_5dose_tot, current_cases_5dose[(1+2*max(age_max)):(3*max(age_max))])
  current_cases_5dose_asym <- c(current_cases_5dose_asym, current_cases_5dose[(1+3*max(age_max)):(4*max(age_max))])
  
  annual_current_cases_8dose_clin <- c(annual_current_cases_8dose_clin, sum(current_cases_8dose[1:max(age_max)]))
  annual_current_cases_8dose_sev <- c(annual_current_cases_8dose_sev, sum(current_cases_8dose[(1+max(age_max)):(2*max(age_max))]))
  annual_current_cases_8dose_tot <- c(annual_current_cases_8dose_tot, sum(current_cases_8dose[(1+2*max(age_max)):(3*max(age_max))]))
  annual_current_cases_8dose_asym <- c(annual_current_cases_8dose_asym, sum(current_cases_8dose[(1+3*max(age_max)):(4*max(age_max))]))

  annual_current_cases_5dose_clin <- c(annual_current_cases_5dose_clin, sum(current_cases_5dose[1:max(age_max)]))
  annual_current_cases_5dose_sev <- c(annual_current_cases_5dose_sev, sum(current_cases_5dose[(1+max(age_max)):(2*max(age_max))]))
  annual_current_cases_5dose_tot <- c(annual_current_cases_5dose_tot, sum(current_cases_5dose[(1+2*max(age_max)):(3*max(age_max))]))
  annual_current_cases_5dose_asym <- c(annual_current_cases_5dose_asym, sum(current_cases_5dose[(1+3*max(age_max)):(4*max(age_max))]))
  
  current_cases_reduction_8dose_clin <- c(current_cases_reduction_8dose_clin, (mean(current_cases_no_PMC[1:max(age_max)]) - mean(current_cases_8dose[1:max(age_max)])) / mean(current_cases_no_PMC[1:max(age_max)]) * 100)
  current_cases_reduction_8dose_sev <- c(current_cases_reduction_8dose_sev, (mean(current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))]) - mean(current_cases_8dose[(1+max(age_max)):(2*max(age_max))])) / mean(current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))]) * 100)
  current_cases_reduction_8dose_tot <- c(current_cases_reduction_8dose_tot, (mean(current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))]) - mean(current_cases_8dose[(1+2*max(age_max)):(3*max(age_max))])) / mean(current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))]) * 100)
  current_cases_reduction_8dose_asym <- c(current_cases_reduction_8dose_asym, (mean(current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))]) - mean(current_cases_8dose[(1+3*max(age_max)):(4*max(age_max))])) / mean(current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))]) * 100)
  
  current_cases_reduction_5dose_clin <- c(current_cases_reduction_5dose_clin, (mean(current_cases_no_PMC[1:max(age_max)]) - mean(current_cases_5dose[1:max(age_max)])) / mean(current_cases_no_PMC[1:max(age_max)]) * 100)
  current_cases_reduction_5dose_sev <- c(current_cases_reduction_5dose_sev, (mean(current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))]) - mean(current_cases_5dose[(1+max(age_max)):(2*max(age_max))])) / mean(current_cases_no_PMC[(1+max(age_max)):(2*max(age_max))]) * 100)
  current_cases_reduction_5dose_tot <- c(current_cases_reduction_5dose_tot, (mean(current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))]) - mean(current_cases_5dose[(1+2*max(age_max)):(3*max(age_max))])) / mean(current_cases_no_PMC[(1+2*max(age_max)):(3*max(age_max))]) * 100)
  current_cases_reduction_5dose_asym <- c(current_cases_reduction_5dose_asym, (mean(current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))]) - mean(current_cases_5dose[(1+3*max(age_max)):(4*max(age_max))])) / mean(current_cases_no_PMC[(1+3*max(age_max)):(4*max(age_max))]) * 100)
  
}


cases_no_PMC$clinical <- current_cases_no_PMC_clin
cases_no_PMC$severe <- current_cases_no_PMC_sev
cases_no_PMC$total <- current_cases_no_PMC_tot
cases_no_PMC$asymptomatic <- current_cases_no_PMC_asym

annual_cases_no_PMC$clinical <- annual_current_cases_no_PMC_clin
annual_cases_no_PMC$severe <- annual_current_cases_no_PMC_sev
annual_cases_no_PMC$total <- annual_current_cases_no_PMC_tot
annual_cases_no_PMC$asymptomatic <- annual_current_cases_no_PMC_asym

cases_averted_8dose$clinical <- current_cases_averted_8dose_clin
cases_averted_8dose$severe <- current_cases_averted_8dose_sev
cases_averted_8dose$total <- current_cases_averted_8dose_tot
cases_averted_8dose$asymptomatic <- current_cases_averted_8dose_asym

cases_averted_5dose$clinical <- current_cases_averted_5dose_clin
cases_averted_5dose$severe <- current_cases_averted_5dose_sev
cases_averted_5dose$total <- current_cases_averted_5dose_tot
cases_averted_5dose$asymptomatic <- current_cases_averted_5dose_asym

annual_cases_averted_8dose$clinical <- annual_current_cases_averted_8dose_clin
annual_cases_averted_8dose$severe <- annual_current_cases_averted_8dose_sev
annual_cases_averted_8dose$total <- annual_current_cases_averted_8dose_tot
annual_cases_averted_8dose$asymptomatic <- annual_current_cases_averted_8dose_asym

annual_cases_averted_5dose$clinical <- annual_current_cases_averted_5dose_clin
annual_cases_averted_5dose$severe <- annual_current_cases_averted_5dose_sev
annual_cases_averted_5dose$total <- annual_current_cases_averted_5dose_tot
annual_cases_averted_5dose$asymptomatic <- annual_current_cases_averted_5dose_asym

cases_8dose$clinical <- current_cases_8dose_clin
cases_8dose$severe <- current_cases_8dose_sev
cases_8dose$total <- current_cases_8dose_tot
cases_8dose$asymptomatic <- current_cases_8dose_asym

cases_5dose$clinical <- current_cases_5dose_clin
cases_5dose$severe <- current_cases_5dose_sev
cases_5dose$total <- current_cases_5dose_tot
cases_5dose$asymptomatic <- current_cases_5dose_asym

annual_cases_8dose$clinical <- annual_current_cases_8dose_clin
annual_cases_8dose$severe <- annual_current_cases_8dose_sev
annual_cases_8dose$total <- annual_current_cases_8dose_tot
annual_cases_8dose$asymptomatic <- annual_current_cases_8dose_asym

annual_cases_5dose$clinical <- annual_current_cases_5dose_clin
annual_cases_5dose$severe <- annual_current_cases_5dose_sev
annual_cases_5dose$total <- annual_current_cases_5dose_tot
annual_cases_5dose$asymptomatic <- annual_current_cases_5dose_asym

cases_reduction_8dose$clinical <- current_cases_reduction_8dose_clin
cases_reduction_8dose$severe <- current_cases_reduction_8dose_sev
cases_reduction_8dose$total <- current_cases_reduction_8dose_tot
cases_reduction_8dose$asymptomatic <- current_cases_reduction_8dose_asym

cases_reduction_5dose$clinical <- current_cases_reduction_5dose_clin
cases_reduction_5dose$severe <- current_cases_reduction_5dose_sev
cases_reduction_5dose$total <- current_cases_reduction_5dose_tot
cases_reduction_5dose$asymptomatic <- current_cases_reduction_5dose_asym

##################### REPEAT ALL OF THIS WITH 6 MONTH AGE GROUP RESULTS #########################


############################ 6 MONTH AGE GROUP RESULTS ###############

cases_no_PMC_sixmonth<- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                   age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                   NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                   units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_8dose_sixmonth<- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                  age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                  NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                  units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_5dose_sixmonth<- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                  age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                  NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                  units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_averted_8dose_sixmonth<- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                          age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                          NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                          units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_averted_5dose_sixmonth<- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                          age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                          NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                          units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_reduction_8dose_sixmonth <- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                             age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                             NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                             units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))

cases_reduction_5dose_sixmonth <- data.frame(age_group = rep(age_group_names_sixmonth[1:no_sixmonth_intervals], times=length(area_names)),
                                             age_in_days_midpoint = rep(sixmonth_intervals_midpoint[1:no_sixmonth_intervals],times=length(area_names)),
                                             NAME_2 = rep(area_names, each=length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])),
                                             units=rep("raw cases", length(area_names)*length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals])))



cases_no_PMC_clin <- cases_no_PMC_sev <- cases_no_PMC_tot <- cases_no_PMC_asym <- c()
cases_8dose_clin <- cases_8dose_sev <- cases_8dose_tot <- cases_8dose_asym <- c()
cases_5dose_clin <- cases_5dose_sev <- cases_5dose_tot <- cases_5dose_asym <- c()
cases_averted_8dose_clin <- cases_averted_8dose_sev <- cases_averted_8dose_tot <- cases_averted_8dose_asym <- c()
cases_averted_5dose_clin <- cases_averted_5dose_sev <- cases_averted_5dose_tot <- cases_averted_5dose_asym <- c()
cases_reduction_8dose_clin <- cases_reduction_8dose_sev <- cases_reduction_8dose_tot <- cases_reduction_8dose_asym <- c()
cases_reduction_5dose_clin <- cases_reduction_5dose_sev <- cases_reduction_5dose_tot <- cases_reduction_5dose_asym <- c()


# Cases, cases averted and % reduction in cases for 6 month age groups

for (i in 1:length(area_names)) {
  
  
  current_cases_no_PMC_clin <- (cases_no_PMC %>% filter(NAME_2 == area_names[i]))$clinical
  current_cases_no_PMC_sev <- (cases_no_PMC %>% filter(NAME_2 == area_names[i]))$severe
  current_cases_no_PMC_tot <- (cases_no_PMC %>% filter(NAME_2 == area_names[i]))$total
  current_cases_no_PMC_asym <- (cases_no_PMC %>% filter(NAME_2 == area_names[i]))$asymptomatic
  
  current_cases_8dose_clin <- (cases_8dose %>% filter(NAME_2 == area_names[i]))$clinical
  current_cases_8dose_sev <- (cases_8dose %>% filter(NAME_2 == area_names[i]))$severe
  current_cases_8dose_tot <- (cases_8dose %>% filter(NAME_2 == area_names[i]))$total
  current_cases_8dose_asym <- (cases_8dose %>% filter(NAME_2 == area_names[i]))$asymptomatic
  
  current_cases_5dose_clin <- (cases_5dose %>% filter(NAME_2 == area_names[i]))$clinical
  current_cases_5dose_sev <- (cases_5dose %>% filter(NAME_2 == area_names[i]))$severe
  current_cases_5dose_tot <- (cases_5dose %>% filter(NAME_2 == area_names[i]))$total
  current_cases_5dose_asym <- (cases_5dose %>% filter(NAME_2 == area_names[i]))$asymptomatic
  
  
  # split current area into each infection class
  split_current_cases_no_PMC_clin <- split(current_cases_no_PMC_clin, ceiling(seq_along(current_cases_no_PMC_clin)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_no_PMC_clin) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_no_PMC_sev <- split(current_cases_no_PMC_sev, ceiling(seq_along(current_cases_no_PMC_sev)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_no_PMC_sev) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_no_PMC_tot <- split(current_cases_no_PMC_tot, ceiling(seq_along(current_cases_no_PMC_tot)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_no_PMC_tot) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_no_PMC_asym <- split(current_cases_no_PMC_asym, ceiling(seq_along(current_cases_no_PMC_asym)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_no_PMC_asym) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  
  # Cases by 6 month age groups (8 dose PMC)
  split_current_cases_8dose_clin <- split(current_cases_8dose_clin, ceiling(seq_along(current_cases_8dose_clin)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_8dose_clin) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_8dose_sev <- split(current_cases_8dose_sev, ceiling(seq_along(current_cases_8dose_sev)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_8dose_sev) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_8dose_tot <- split(current_cases_8dose_tot, ceiling(seq_along(current_cases_8dose_tot)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_8dose_tot) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_8dose_asym <- split(current_cases_8dose_asym, ceiling(seq_along(current_cases_8dose_asym)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_8dose_asym) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  # Cases by 6 month age groups (5 dose PMC)
  split_current_cases_5dose_clin <- split(current_cases_5dose_clin, ceiling(seq_along(current_cases_5dose_clin)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_5dose_clin) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_5dose_sev <- split(current_cases_5dose_sev, ceiling(seq_along(current_cases_5dose_sev)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_5dose_sev) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_5dose_tot <- split(current_cases_5dose_tot, ceiling(seq_along(current_cases_5dose_tot)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_5dose_tot) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  split_current_cases_5dose_asym <- split(current_cases_5dose_asym, ceiling(seq_along(current_cases_5dose_asym)/(max(age_max)/length(sixmonth_intervals_midpoint[1:no_sixmonth_intervals]))))
  names(split_current_cases_5dose_asym) <- age_group_names_sixmonth[1:no_sixmonth_intervals]
  
  for (j in 1:length(sixmonth_intervals[1:no_sixmonth_intervals])) {
    cases_no_PMC_clin<- c(cases_no_PMC_clin, sum(split_current_cases_no_PMC_clin[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_no_PMC_sev<- c(cases_no_PMC_sev, sum(split_current_cases_no_PMC_sev[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_no_PMC_tot<- c(cases_no_PMC_tot, sum(split_current_cases_no_PMC_tot[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_no_PMC_asym<- c(cases_no_PMC_asym, sum(split_current_cases_no_PMC_asym[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    
    cases_8dose_clin<- c(cases_8dose_clin, sum(split_current_cases_8dose_clin[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_8dose_sev<- c(cases_8dose_sev, sum(split_current_cases_8dose_sev[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_8dose_tot<- c(cases_8dose_tot, sum(split_current_cases_8dose_tot[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_8dose_asym<- c(cases_8dose_asym, sum(split_current_cases_8dose_asym[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    
    cases_5dose_clin<- c(cases_5dose_clin, sum(split_current_cases_5dose_clin[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_5dose_sev<- c(cases_5dose_sev, sum(split_current_cases_5dose_sev[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_5dose_tot<- c(cases_5dose_tot, sum(split_current_cases_5dose_tot[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    cases_5dose_asym<- c(cases_5dose_asym, sum(split_current_cases_5dose_asym[[age_group_names_sixmonth[1:no_sixmonth_intervals][j]]]))
    
  }
  
}  

cases_averted_8dose_clin <- (cases_no_PMC_clin - cases_8dose_clin)
cases_averted_8dose_sev <- (cases_no_PMC_sev - cases_8dose_sev)
cases_averted_8dose_tot <- (cases_no_PMC_tot - cases_8dose_tot)
cases_averted_8dose_asym <- (cases_no_PMC_asym - cases_8dose_asym)

cases_averted_5dose_clin <- (cases_no_PMC_clin - cases_5dose_clin)
cases_averted_5dose_sev <- (cases_no_PMC_sev - cases_5dose_sev)
cases_averted_5dose_tot <- (cases_no_PMC_tot - cases_5dose_tot)
cases_averted_5dose_asym <- (cases_no_PMC_asym - cases_5dose_asym)

cases_reduction_8dose_clin <- (cases_no_PMC_clin - cases_8dose_clin)/cases_no_PMC_clin * 100
cases_reduction_8dose_sev <- (cases_no_PMC_sev - cases_8dose_sev)/cases_no_PMC_sev * 100
cases_reduction_8dose_tot <- (cases_no_PMC_tot - cases_8dose_tot)/cases_no_PMC_tot * 100
cases_reduction_8dose_asym <- (cases_no_PMC_asym - cases_8dose_asym)/cases_no_PMC_asym * 100

cases_reduction_5dose_clin <- (cases_no_PMC_clin - cases_5dose_clin)/cases_no_PMC_clin * 100
cases_reduction_5dose_sev <- (cases_no_PMC_sev - cases_5dose_sev)/cases_no_PMC_sev * 100
cases_reduction_5dose_tot <- (cases_no_PMC_tot - cases_5dose_tot)/cases_no_PMC_tot * 100
cases_reduction_5dose_asym <- (cases_no_PMC_asym - cases_5dose_asym)/cases_no_PMC_asym * 100

# append results to dataframe
cases_no_PMC_sixmonth$clinical <- cases_no_PMC_clin
cases_no_PMC_sixmonth$severe <- cases_no_PMC_sev
cases_no_PMC_sixmonth$total <- cases_no_PMC_tot
cases_no_PMC_sixmonth$asymptomatic <- cases_no_PMC_asym

cases_8dose_sixmonth$clinical<- cases_8dose_clin
cases_8dose_sixmonth$severe<- cases_8dose_sev
cases_8dose_sixmonth$total<- cases_8dose_tot
cases_8dose_sixmonth$asymptomatic<- cases_8dose_asym

cases_5dose_sixmonth$clinical<- cases_5dose_clin
cases_5dose_sixmonth$severe<- cases_5dose_sev
cases_5dose_sixmonth$total<- cases_5dose_tot
cases_5dose_sixmonth$asymptomatic<- cases_5dose_asym

cases_averted_8dose_sixmonth$clinical <- cases_averted_8dose_clin
cases_averted_8dose_sixmonth$severe <- cases_averted_8dose_sev
cases_averted_8dose_sixmonth$total <- cases_averted_8dose_tot
cases_averted_8dose_sixmonth$asymptomatic <- cases_averted_8dose_asym

cases_averted_5dose_sixmonth$clinical <- cases_averted_5dose_clin
cases_averted_5dose_sixmonth$severe <- cases_averted_5dose_sev
cases_averted_5dose_sixmonth$total <- cases_averted_5dose_tot
cases_averted_5dose_sixmonth$asymptomatic <- cases_averted_5dose_asym

cases_reduction_8dose_sixmonth$clinical <- cases_reduction_8dose_clin
cases_reduction_8dose_sixmonth$severe <- cases_reduction_8dose_sev
cases_reduction_8dose_sixmonth$total <- cases_reduction_8dose_tot
cases_reduction_8dose_sixmonth$asymptomatic <- cases_reduction_8dose_asym

cases_reduction_5dose_sixmonth$clinical <- cases_reduction_5dose_clin
cases_reduction_5dose_sixmonth$severe <- cases_reduction_5dose_sev
cases_reduction_5dose_sixmonth$total <- cases_reduction_5dose_tot
cases_reduction_5dose_sixmonth$asymptomatic <- cases_reduction_5dose_asym


############################################################################################################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


# TO PRODUCE GRAPHS FOR ANDRIA








##############################################################################
##############################################################################

# ggplot() + 
#   geom_line(aes(x=1:max(age_max), y=cases_no_PMC$Sud[1:max(age_max)], colour="no pmc"))+
#   geom_line(aes(x=1:max(age_max), y=cases_8dose$Sud[1:max(age_max)], colour="pmc"))
# 
# ggplot() + 
#   geom_line(aes(x=1:max(age_max), y=cases_no_PMC$Nord[1:max(age_max)], colour="no pmc"))+
#   geom_line(aes(x=1:max(age_max), y=cases_8dose$Nord[1:max(age_max)], colour="pmc"))
# 
# ggplot() + 
#   geom_line(aes(x=1:max(age_max), y=cases_no_PMC$Sud_Ouest[1:max(age_max)], colour="no pmc"))+
#   geom_line(aes(x=1:max(age_max), y=cases_8dose$Sud_Ouest[1:max(age_max)], colour="pmc"))


################################ GRAPHS #######################################

# CLINICAL INCIDENCE REDUCTION %

CMR_adm1<- sf::st_read("CMR_shapefiles/gadm41_CMR_1.shp")
CMR_adm1$NAME_2 <- stri_trans_general(str=gsub("-", "_", CMR_adm1$NAME_1), id = "Latin-ASCII")

# link data with shapefile
combined_df_8dose_reduction_CMR <- merge(cases_reduction_8dose,CMR_adm1, by="NAME_2") 
combined_df_5dose_reduction_CMR <- merge(cases_reduction_5dose,CMR_adm1, by="NAME_2") 

impact_graph_8dose <- ggplot(combined_df_8dose_reduction_CMR) + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle("8 DOSE") + labs(fill="% reduction") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

impact_graph_5dose <- ggplot(combined_df_5dose_reduction_CMR)  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle("5 DOSE") + labs(fill="% reduction") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# impact_graph_8dose <- ggplot(combined_df_8dose_reduction_CMR) + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#   theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#   geom_sf(data=CMR_adm1)+
#   geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#   geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#   ggtitle("8 DOSE") + labs(fill="% reduction") +
#   scale_fill_gradient(low="red", high="darkgreen", limits = c(0,15),breaks=seq(0, 50, by = 10) ) + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# impact_graph_5dose <- ggplot(combined_df_5dose_reduction_CMR)  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#   theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#   geom_sf(data=CMR_adm1)+
#   geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#   geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#   ggtitle("5 DOSE") + labs(fill="% reduction") +
#   scale_fill_gradient(low="red", high="darkgreen", limits = c(0,15),breaks=seq(0, 50, by = 10) ) + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())

plot_impact <- ggarrange(impact_graph_8dose, impact_graph_5dose, nrow=1, ncol=2, common.legend = TRUE, legend = "right")
annotate_figure(plot_impact, top = text_grob(paste0("Reduction in clinical cases in ", country_code, " for 8 and 5 dose PMC schedules")))

ggsave(paste0(getwd(), "/simulation_results/program_comparison_82cov_.jpeg"), units="in", width=8, height=5, dpi=300)


# repeat for 6 month intervals

# 8 dose

# link data with shapefile
combined_df_8dose_reduction_CMR_sixmonth <- merge(cases_reduction_8dose_sixmonth,CMR_adm1, by="NAME_2") 

impact_graph_8dose_0_182 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_0_182"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle(paste0("Age: 0-6 month")) + labs(fill="% reduction") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

if (no_sixmonth_intervals > 1) {
  impact_graph_8dose_183_364 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_183_364"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age: 6-12 months")) + labs(fill="% reduction") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 2) {
  impact_graph_8dose_365_547 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_365_547"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age: 12-30 months")) + labs(fill="% reduction") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 3) {
  impact_graph_8dose_548_729 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_548_729"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age: 12-30 months")) + labs(fill="% reduction") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 4) {
  impact_graph_8dose_730_912 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_730_912"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[5],"-", sixmonth_intervals[6], " days")) + labs(fill="% reduction") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

# impact_graph_8dose_0_182 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_0_182"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#   theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#   geom_sf(data=CMR_adm1)+
#   geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#   geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#   ggtitle(paste0("Age group: ", sixmonth_intervals[1],"-", sixmonth_intervals[2], " days")) + labs(fill="% reduction") +
#   scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# if (no_sixmonth_intervals > 1) {
#   impact_graph_8dose_183_364 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_183_364"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#     theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#     geom_sf(data=CMR_adm1)+
#     geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#     geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#     ggtitle(paste0("Age group: ", sixmonth_intervals[2],"-", sixmonth_intervals[3], " days")) + labs(fill="% reduction") +
#     scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
# }
# 
# if (no_sixmonth_intervals > 2) {
#   impact_graph_8dose_365_547 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_365_547"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#     theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#     geom_sf(data=CMR_adm1)+
#     geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#     geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#     ggtitle(paste0("Age group: ", sixmonth_intervals[3],"-", sixmonth_intervals[4], " days")) + labs(fill="% reduction") +
#     scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
# }
# 
# if (no_sixmonth_intervals > 3) {
#   impact_graph_8dose_548_729 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_548_729"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#     theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#     geom_sf(data=CMR_adm1)+
#     geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#     geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#     ggtitle(paste0("Age group: ", sixmonth_intervals[4],"-", sixmonth_intervals[5], " days")) + labs(fill="% reduction") +
#     scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
# }
# 
# if (no_sixmonth_intervals > 4) {
#   impact_graph_8dose_730_912 <- ggplot(combined_df_8dose_reduction_CMR_sixmonth %>% filter(age_group == "n_age_730_912"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
#     theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
#     geom_sf(data=CMR_adm1)+
#     geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
#     geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2), "%"), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
#     ggtitle(paste0("Age group: ", sixmonth_intervals[5],"-", sixmonth_intervals[6], " days")) + labs(fill="% reduction") +
#     scale_fill_gradient(low="red", high="darkgreen", limits = c(0,50),breaks=seq(0, 50, by = 10) ) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
# }

plot_impact <- ggarrange(impact_graph_8dose_0_182,impact_graph_8dose_183_364, impact_graph_8dose_548_729, nrow=1, ncol=3, common.legend = TRUE, legend = "right")
annotate_figure(plot_impact, top = text_grob(paste0("Reduction in clinical cases in ", country_code, " during 8 dose PMC schedule")))

ggsave(paste0(getwd(), "/simulation_results/program_comparison_by_age_groups_82cov.jpeg"), units="in", width=8, height=5, dpi=300)

# ANNUAL CLINICAL CASES AVERTED (RAW)

##link data with shapefile
combined_df_8dose_averted_CMR <- merge(annual_cases_averted_8dose,CMR_adm1, by="NAME_2") 
combined_df_5dose_averted_CMR <- merge(annual_cases_averted_5dose,CMR_adm1, by="NAME_2") 

impact_graph_8dose <- ggplot(combined_df_8dose_averted_CMR)  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle("8 DOSE") + labs(fill="Cases averted") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,100000),breaks=seq(0, 100000, by = 10000) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

impact_graph_5dose <- ggplot(combined_df_5dose_averted_CMR)  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle("5 DOSE") + labs(fill="Cases averted") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,100000),breaks=seq(0, 100000, by = 10000) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot_impact <- ggarrange(impact_graph_8dose, impact_graph_5dose, nrow=1, ncol=2, common.legend = TRUE, legend = "right")
annotate_figure(plot_impact, top = text_grob(paste0("Clinical cases averted in ", country_code, " during 8 dose PMC schedule")))


# repeat for 6 month intervals

# 8 dose

# link data with shapefile
combined_df_8dose_averted_CMR_sixmonth <- merge(cases_averted_8dose_sixmonth,CMR_adm1, by="NAME_2") 

impact_graph_8dose_0_182 <- ggplot(combined_df_8dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_0_182"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle(paste0("Age group: ", sixmonth_intervals[1],"-", sixmonth_intervals[2], " days")) + labs(fill="% averted") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_8dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_8dose_averted_CMR_sixmonth$clinical), by = 10000) ) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


if (no_sixmonth_intervals > 1) {
  impact_graph_8dose_183_364 <- ggplot(combined_df_8dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_183_364"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[2],"-", sixmonth_intervals[3], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_8dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_8dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 2) {
  impact_graph_8dose_365_547 <- ggplot(combined_df_8dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_365_547"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[3],"-", sixmonth_intervals[4], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_8dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_8dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 3) {
  impact_graph_8dose_548_729 <- ggplot(combined_df_8dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_548_729"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[4],"-", sixmonth_intervals[5], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_8dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_8dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 4) {
  impact_graph_8dose_730_912 <- ggplot(combined_df_8dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_730_912"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[5],"-", sixmonth_intervals[6], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_8dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_8dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

plot_impact <- ggarrange(impact_graph_8dose_0_182,impact_graph_8dose_183_364, impact_graph_8dose_365_547, impact_graph_8dose_548_729, impact_graph_8dose_730_912, nrow=2, ncol=3, common.legend = TRUE, legend = "right")
annotate_figure(plot_impact, top = text_grob(paste0("Clinical cases averted in ", country_code, " during 8 dose PMC schedule")))

##################################################################################################

# 5 dose

# link data with shapefile
combined_df_5dose_averted_CMR_sixmonth <- merge(cases_averted_5dose_sixmonth,CMR_adm1, by="NAME_2") 

impact_graph_5dose_0_182 <- ggplot(combined_df_5dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_0_182"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
  theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
  geom_sf(data=CMR_adm1)+
  geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
  geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
  ggtitle(paste0("Age group: ", sixmonth_intervals[1],"-", sixmonth_intervals[2], " days")) + labs(fill="% averted") +
  scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_5dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_5dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

if (no_sixmonth_intervals > 1) {
  impact_graph_5dose_183_364 <- ggplot(combined_df_5dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_183_364"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[2],"-", sixmonth_intervals[3], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_5dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_5dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 2) {
  impact_graph_5dose_365_547 <- ggplot(combined_df_5dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_365_547"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[3],"-", sixmonth_intervals[4], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_5dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_5dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 3) {
  impact_graph_5dose_548_729 <- ggplot(combined_df_5dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_548_729"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[4],"-", sixmonth_intervals[5], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_5dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_5dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

if (no_sixmonth_intervals > 4) {
  impact_graph_5dose_730_912 <- ggplot(combined_df_5dose_averted_CMR_sixmonth %>% filter(age_group == "n_age_730_912"))  + theme_bw()+  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), title = element_text(size=10))+
    theme(legend.key.height = unit(1.5, "cm"),legend.key.width = unit(0.5, "cm")) +
    geom_sf(data=CMR_adm1)+
    geom_sf(aes(fill=clinical, geometry=geometry), size=0.5)  +
    geom_sf_label(aes(label = paste0(NAME_1, ": ", round(clinical, 2)), geometry=geometry),fun.geometry = st_centroid, size=3, alpha = 0.2) +
    ggtitle(paste0("Age group: ", sixmonth_intervals[5],"-", sixmonth_intervals[6], " days")) + labs(fill="% averted") +
    scale_fill_gradient(low="red", high="darkgreen", limits = c(0,max(combined_df_5dose_averted_CMR_sixmonth$clinical)),breaks=seq(0, max(combined_df_5dose_averted_CMR_sixmonth$clinical), by = 10000) ) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}

plot_impact <- ggarrange(impact_graph_5dose_0_182,impact_graph_5dose_183_364, impact_graph_5dose_365_547, impact_graph_5dose_548_729, impact_graph_5dose_730_912, nrow=2, ncol=3, common.legend = TRUE, legend = "right")
annotate_figure(plot_impact, top = text_grob(paste0("Clinical cases averted in ", country_code, " during 5 dose PMC schedule")))
