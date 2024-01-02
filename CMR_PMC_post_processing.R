#source("CMR_model.R")
#source("running_model_CMR.R")
#source("CMR_merge_results.R")

###################### READ IN RESULTS FROM IMPERIAL MODEL ####################

#sim_output <- read.csv(paste0(getwd(), "/simulation_results/sim_output_merged.csv"))
incidence_ppy_df <- read.csv(paste0(getwd(), "/simulation_results/incidence_ppy_df_merged.csv"))
#incidence_df <- read.csv(paste0(getwd(), "/simulation_results/incidence_df_merged.csv"))
population_df <- read.csv(paste0(getwd(), "/simulation_results/population_df_merged.csv"))

########################### SET UP EFFICACY CURVES #############################

# Weibull scale parameters for each haplotype
lambda_trip<-59.57659
lambda_quadr<-33.05391
lambda_quint<-18.55328
lambda_sext<-12.81186
lambda_other<-23
lambda_VAGKAA<-18.55328 ### COMPLETE GUESS that it is the same as the QUINTUPLE
lambda_VAGKGS<-12.81186 ### COMPLETE GUESS that it is the same as the SEXTUPLE

# Weibull shape parameters for each haplotype
w_trip<- 8.435971
w_quadr<-4.862126
w_quint<-2.4840752
w_sext<-3.691953
w_other<-4.5
w_VAGKAA<-1.7
w_VAGKGS<- 4.1


# PMC schedules in CMR
# 8 doses for CMR: 10wk, 14wk, 6mo,  9mo, 12mo, 15mo, 18mo, 24mo
# 5 doses for CMR: 10wk, 14wk, 6mo,  9mo, 15mo
schedule_CMR_8doses<- c(10*7, 14*7, 6*30, 9*30, 12*30, 15*30, 18*30, 24*30)
schedule_CMR_5doses<- c(10*7, 14*7, 6*30, 9*30, 15*30)

# PMC coverages in CMR 
cov_CMR_8doses<- c(0.82, 0.73, 0.26, 0.52, 0.15 , 0.24 ,0.09, 0.03)
#cov_CMR_8doses<- c(0.82, 0.82, 0.82, 0.82, 0.82 , 0.82 ,0.82, 0.82)

cov_CMR_5doses<- c(0.82, 0.73,0.26, 0.52, 0.24)
#cov_CMR_5doses<- c(0.82, 0.82,0.82, 0.82, 0.82)


# set up df with the timesteps 
if (sim_length > length(age_in_days_midpoint)) {
  df<- data.frame(time=1:sim_length)
}

if (sim_length < length(age_in_days_midpoint)) {
  df<- data.frame(time=1:length(age_in_days_midpoint))
}

df$prot_trip<-NA
df$prot_quadr<-NA
df$prot_quint<-NA
df$prot_sext<-NA
df$prot_other<-NA
df$prot_VAGKAA<-NA
df$prot_VAGKGS<-NA

 
# Construct weibull curves to model the efficacy of SP through time for 
# each haplotype 

for (t in 1: schedule_CMR_8doses[1]) {  # day 0 to dose 1 on day 70 
  df$prot_trip[t]<-0
  df$prot_quadr[t]<-0
  df$prot_quint[t]<-0
  df$prot_sext[t]<- 0
  df$prot_other[t]<-0
  df$prot_VAGKAA[t]<- 0
  df$prot_VAGKGS[t]<- 0  
}


for (t in (schedule_CMR_8doses[1]+1) : schedule_CMR_8doses[2])  {   # day 71 to dose 2 on day 98
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_trip)^w_trip) * cov_CMR_8doses[1]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[1]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_quint)^w_quint)* cov_CMR_8doses[1]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_sext)^w_sext)* cov_CMR_8doses[1]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_other)^w_other)* cov_CMR_8doses[1]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[1]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[1]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[1]
}

for (t in (schedule_CMR_8doses[2]+1) : schedule_CMR_8doses[3])  {  #  day 99 to dose 3 on day 180
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_trip)^w_trip)* cov_CMR_8doses[2]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[2]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_quint)^w_quint)* cov_CMR_8doses[2]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_sext)^w_sext)* cov_CMR_8doses[2]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_other)^w_other)* cov_CMR_8doses[2]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[2]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[2]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[2]
}

for (t in (schedule_CMR_8doses[3]+1) : schedule_CMR_8doses[4])  {  # day 181 to dose 4 on day 270
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_trip)^w_trip)* cov_CMR_8doses[3]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[3]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_quint)^w_quint)* cov_CMR_8doses[3]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_sext)^w_sext)* cov_CMR_8doses[3]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_other)^w_other)* cov_CMR_8doses[3]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[3]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[3]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[3]
}

for (t in (schedule_CMR_8doses[4]+1) : schedule_CMR_8doses[5])  {  # day 271 to dose 5 on day 360 
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_trip)^w_trip)* cov_CMR_8doses[4]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[4]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_quint)^w_quint)* cov_CMR_8doses[4]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_sext)^w_sext)* cov_CMR_8doses[4]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_other)^w_other)* cov_CMR_8doses[4]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[4]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[4]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[4]
}

for (t in (schedule_CMR_8doses[5]+1) : schedule_CMR_8doses[6])  {  # day 361 to dose 6 on day 450 
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_trip)^w_trip)* cov_CMR_8doses[5]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[5]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_quint)^w_quint)* cov_CMR_8doses[5]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_sext)^w_sext)* cov_CMR_8doses[5]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_other)^w_other)* cov_CMR_8doses[5]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[5]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[5]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[5]
}

for (t in (schedule_CMR_8doses[6]+1) : schedule_CMR_8doses[7])  {  # day 451 to dose 7 on day 540 
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_trip)^w_trip)* cov_CMR_8doses[6]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[6]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_quint)^w_quint)* cov_CMR_8doses[6]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_sext)^w_sext)* cov_CMR_8doses[6]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_other)^w_other)* cov_CMR_8doses[6]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[6]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[6]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[6]
}

for (t in (schedule_CMR_8doses[7]+1) : schedule_CMR_8doses[8])  {  # day 541 to dose 8 on day 720 
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_trip)^w_trip)* cov_CMR_8doses[7]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[7]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_quint)^w_quint)* cov_CMR_8doses[7]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_sext)^w_sext)* cov_CMR_8doses[7]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_other)^w_other)* cov_CMR_8doses[7]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[7]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[7]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[7]
}

for (t in (schedule_CMR_8doses[8]+1) : nrow(df))  {  # day 721 to end of simulation   
  df$prot_trip[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_trip)^w_trip)* cov_CMR_8doses[8]
  df$prot_quadr[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_quadr)^w_quadr)* cov_CMR_8doses[8]
  df$prot_quint[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_quint)^w_quint)* cov_CMR_8doses[8]
  df$prot_sext[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_sext)^w_sext)* cov_CMR_8doses[8]
  df$prot_other[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_other)^w_other)* cov_CMR_8doses[8]
  df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_8doses[8]
  df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule_CMR_8doses[8]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_8doses[8]
}


##############################################

# Repeat for 5 dose schedule:

df$prot_trip_5doses<-NA
df$prot_quadr_5doses<-NA
df$prot_quint_5doses<-NA
df$prot_sext_5doses<-NA
df$prot_other_5doses<-NA
df$prot_VAGKAA_5doses<-NA
df$prot_VAGKGS_5doses<-NA

for (t in 1: schedule_CMR_5doses[1]) {  # day 0 to dose 1 on day 70 
  df$prot_trip_5doses[t]<-0
  df$prot_quadr_5doses[t]<-0
  df$prot_quint_5doses[t]<-0
  df$prot_sext_5doses[t]<- 0
  df$prot_other_5doses[t]<-0
  df$prot_VAGKAA_5doses[t]<- 0
  df$prot_VAGKGS_5doses[t]<- 0  
}



for (t in (schedule_CMR_5doses[1]+1) : schedule_CMR_5doses[2])  {   # day 71 to dose 2 on day 98
  df$prot_trip_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_trip)^w_trip) * cov_CMR_5doses[1]
  df$prot_quadr_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_quadr)^w_quadr)* cov_CMR_5doses[1]
  df$prot_quint_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_quint)^w_quint)* cov_CMR_5doses[1]
  df$prot_sext_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_sext)^w_sext)* cov_CMR_5doses[1]
  df$prot_other_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_other)^w_other)* cov_CMR_5doses[1]
  df$prot_VAGKAA_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_5doses[1]
  df$prot_VAGKGS_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[1]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_5doses[1]
}

for (t in (schedule_CMR_5doses[2]+1) : schedule_CMR_5doses[3])  {  # day 99 to dose 3 on day 180
  df$prot_trip_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_trip)^w_trip) * cov_CMR_5doses[2]
  df$prot_quadr_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_quadr)^w_quadr)* cov_CMR_5doses[2]
  df$prot_quint_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_quint)^w_quint)* cov_CMR_5doses[2]
  df$prot_sext_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_sext)^w_sext)* cov_CMR_5doses[2]
  df$prot_other_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_other)^w_other)* cov_CMR_5doses[2]
  df$prot_VAGKAA_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_5doses[2]
  df$prot_VAGKGS_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[2]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_5doses[2]
}

for (t in (schedule_CMR_5doses[3]+1) : schedule_CMR_5doses[4])  {  # day 181 to dose 270 on day 
  df$prot_trip_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_trip)^w_trip)* cov_CMR_5doses[3]
  df$prot_quadr_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_quadr)^w_quadr)* cov_CMR_5doses[3]
  df$prot_quint_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_quint)^w_quint)* cov_CMR_5doses[3]
  df$prot_sext_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_sext)^w_sext)* cov_CMR_5doses[3]
  df$prot_other_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_other)^w_other)* cov_CMR_5doses[3]
  df$prot_VAGKAA_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_5doses[3]
  df$prot_VAGKGS_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[3]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_5doses[3]
}

for (t in (schedule_CMR_5doses[4]+1) : schedule_CMR_5doses[5])  {  # day 271 to dose 5 on day 450
  df$prot_trip_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_trip)^w_trip)* cov_CMR_5doses[4]
  df$prot_quadr_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_quadr)^w_quadr)* cov_CMR_5doses[4]
  df$prot_quint_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_quint)^w_quint)* cov_CMR_5doses[4]
  df$prot_sext_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_sext)^w_sext)* cov_CMR_5doses[4]
  df$prot_other_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_other)^w_other)* cov_CMR_5doses[4]
  df$prot_VAGKAA_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_5doses[4]
  df$prot_VAGKGS_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[4]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_5doses[4]
}

for (t in (schedule_CMR_5doses[5]+1) : nrow(df))  {  # day 451 to end of simulation 
  df$prot_trip_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_trip)^w_trip)* cov_CMR_5doses[5]
  df$prot_quadr_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_quadr)^w_quadr)* cov_CMR_5doses[5]
  df$prot_quint_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_quint)^w_quint)* cov_CMR_5doses[5]
  df$prot_sext_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_sext)^w_sext)* cov_CMR_5doses[5]
  df$prot_other_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_other)^w_other)* cov_CMR_5doses[5]
  df$prot_VAGKAA_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_VAGKAA)^w_VAGKAA)* cov_CMR_5doses[5]
  df$prot_VAGKGS_5doses[t]<- exp(-(df$time[t-schedule_CMR_5doses[5]]/lambda_VAGKGS)^w_VAGKGS)* cov_CMR_5doses[5]
}


colors <- c("triple" = "#00C094", "quadruple" = "#00B6EB", "quintuple" = "#FFA500", 
            "VAGKAA"= "gold1", "sextuple" = "#F8766D" , "VAGKGS"= "#FB61D7", "other"= "#A58AFF")

# plot graphs for the probability of drug protection for each genotype through
# time during the SP dosing schedule given varying levels of coverage. The
# graph shapes will be the same but since coverage varies, the height differs

ggplot(data=df)+ theme_bw() +
  geom_line(aes(x=time, y=prot_trip,color="I_AKA_")) + 
  geom_line(aes(x=time, y=prot_quadr,color="I_GKA_")) + 
  geom_line(aes(x=time, y=prot_quint,color="I_GEA_")) + 
  geom_line(aes(x=time, y=prot_VAGKAA,color= "V_GKA_")) + 
  geom_line(aes(x=time, y=prot_sext, color="I_GEG_")) + 
  geom_line(aes(x=time, y=prot_VAGKGS, color="V_GKG_")) + 
  geom_line(aes(x=time, y=prot_other, color="other")) + 
  geom_vline(xintercept = schedule_CMR_8doses, color="blue", linetype="dashed") +
  geom_vline(xintercept = schedule_CMR_5doses, color="aquamarine4", linetype="dashed") +
  labs(color="Resistance genotype") +
  ylab("Probability of drug protection")+xlab("Age(days)")

ggplot(data=df)+ theme_bw() +
  geom_line(aes(x=time, y=prot_trip_5doses,color="I_AKA_")) + 
  geom_line(aes(x=time, y=prot_quadr_5doses,color="I_GKA_")) + 
  geom_line(aes(x=time, y=prot_quint_5doses,color="I_GEA_")) + 
  geom_line(aes(x=time, y=prot_VAGKAA_5doses,color= "V_GKA_")) + 
  geom_line(aes(x=time, y=prot_sext_5doses, color="I_GEG_")) + 
  geom_line(aes(x=time, y=prot_VAGKGS_5doses, color="V_GKG_")) + 
  geom_line(aes(x=time, y=prot_other_5doses, color="other")) + 
  geom_vline(xintercept = schedule_CMR_5doses, color="aquamarine4", linetype="dashed") +
  labs(color="Resistance genotype") +
  ylab("Probability of drug protection")+xlab("Age(days)")


# Read in proportions of each haplotype

full_prop_data<-read_excel("C:/Users/IshanaAtkinson/OneDrive - London School of Hygiene and Tropical Medicine/Documents/malariasimulation work/genetic_data.xlsx", sheet = "proportions")

# weighted means to calculate overall protection by area

for (i in 1:length(area_names)){
  proportions <- full_prop_data %>% filter(name_2 == area_names[i])
  df[paste0("prot_overall_", area_names[i])] <- proportions$prop_trip*df$prot_trip +
    proportions$prop_quadr*df$prot_quadr +
    proportions$prop_quint*df$prot_quint +
    proportions$prop_sext*df$prot_sext +
    proportions$prop_VAGKAA*df$prot_VAGKAA +
    proportions$prop_VAGKGS*df$prot_VAGKGS +
    proportions$prop_other*df$prot_other
  
  df[paste0("prot_overall_5doses_", area_names[i])]<- proportions$prop_trip*df$prot_trip_5doses +
    proportions$prop_quadr*df$prot_quadr_5doses +
    proportions$prop_quint*df$prot_quint_5doses +
    proportions$prop_sext*df$prot_sext_5doses +
    proportions$prop_VAGKAA*df$prot_VAGKAA_5doses +
    proportions$prop_VAGKGS*df$prot_VAGKGS_5doses +
    proportions$prop_other*df$prot_other_5doses

}


colours_sites<-c("South"= "blue", "SouthWest"="darkgreen", "North"="maroon") 

graph_8dose <- ggplot()+ theme_bw() +
  geom_line(data=df, aes(x=time, y=prot_overall_Sud,color="South")) +
  #geom_line(data=df, aes(x=time, y=prot_overall_Sud_Ouest,color="SouthWest")) +
  geom_line(data=df, aes(x=time, y=prot_overall_Nord,color="North")) +
  scale_color_manual("Site", values = colours_sites) +
  ylab("Prob. of drug protection")+xlab("Age (days)")

graph_5dose <- ggplot()+ theme_bw() +
  geom_line(data=df, aes(x=time, y=prot_overall_5doses_Sud,color="South")) +
  #geom_line(data=df, aes(x=time, y=prot_overall_5doses_Sud_Ouest,color="SouthWest")) +
  geom_line(data=df, aes(x=time, y=prot_overall_5doses_Nord,color="North")) +
  scale_color_manual("Site", values = colours_sites) +
  ylab("Prob. of drug protection")+xlab("Age (days)")


ggarrange(graph_8dose, graph_5dose, ncol=1, nrow=2)

##############################################################################
##################################### CLINICAL #################################
##############################################################################

PMC_impact_ppy_8dose <- data.frame()
PMC_impact_ppy_5dose <- data.frame()

# convert this to match new format of dataframes 

for (i in 1:length(area_names)){
  
  # PMC impact on incidence (ppy)
  new_PMC_impact_ppy_8dose_df <- (incidence_ppy_df %>% filter(area == area_names[i]))   
  new_PMC_impact_ppy_5dose_df <- (incidence_ppy_df %>% filter(area == area_names[i]))   
  
  PMC_impact_ppy_8dose_cases <- new_PMC_impact_ppy_8dose_df$value * rep((1-df[paste0("prot_overall_", area_names[i])][min(age_min):max(age_max),]),4) 
  PMC_impact_ppy_5dose_cases <- new_PMC_impact_ppy_5dose_df$value * rep((1-df[paste0("prot_overall_5doses_", area_names[i])][min(age_min):max(age_max),]),4) 
  
  new_PMC_impact_ppy_8dose_df$value <- PMC_impact_ppy_8dose_cases
  new_PMC_impact_ppy_5dose_df$value <- PMC_impact_ppy_5dose_cases
  
  PMC_impact_ppy_8dose <- rbind(PMC_impact_ppy_8dose, new_PMC_impact_ppy_8dose_df)
  PMC_impact_ppy_5dose <- rbind(PMC_impact_ppy_5dose, new_PMC_impact_ppy_5dose_df)
}

# Average clinical incidence per person per year for each age group
# Graph shows only the ages we are looking at and the PMC doses that impact 
# that age group

ggplot() +
  geom_line(aes(x=unique(PMC_impact_ppy_8dose$age_in_days_midpoint), y=as.double((PMC_impact_ppy_8dose %>% filter(area=="Sud" & infection_class == "clinical"))$value), colour = "Sud"))+
  geom_line(aes(x=unique(PMC_impact_ppy_8dose$age_in_days_midpoint), y=as.double((PMC_impact_ppy_8dose %>% filter(area=="Nord" & infection_class == "clinical"))$value), colour = "Nord"))+
  geom_vline(xintercept = schedule_CMR_8doses[schedule_CMR_8doses < max(age_max)], color="aquamarine4", linetype="dashed") +
  geom_text(aes(x=schedule_CMR_8doses[schedule_CMR_8doses < max(age_max)], label=cov_CMR_8doses[1:length(schedule_CMR_8doses[schedule_CMR_8doses < max(age_max)])], y=2),
            colour="black", angle=90, size = 2.25, nudge_x = -20) +
  labs(title="Average clinical incidence each year per age group with PMC (8 doses)", x="Age (days)", y="New clinical infections (ppy)",
       colour = "Level-1 admin")

# Average clinical incidence per person per year for each age group
# Graph shows only the ages we are looking at and the PMC doses that impact 
# that age group

ggplot() +
  geom_line(aes(x=unique(PMC_impact_ppy_5dose$age_in_days_midpoint), y=as.double((PMC_impact_ppy_5dose %>% filter(area=="Sud" & infection_class == "clinical"))$value), colour = "Sud"))+
  geom_line(aes(x=unique(PMC_impact_ppy_5dose$age_in_days_midpoint), y=as.double((PMC_impact_ppy_5dose %>% filter(area=="Nord" & infection_class == "clinical"))$value), colour = "Nord"))+
  geom_vline(xintercept = schedule_CMR_5doses[schedule_CMR_5doses < max(age_max)], color="aquamarine4", linetype="dashed") +
  geom_text(aes(x=schedule_CMR_5doses[schedule_CMR_5doses < max(age_max)], label=cov_CMR_5doses[1:length(schedule_CMR_5doses[schedule_CMR_5doses < max(age_max)])], y=2),
            colour="black", angle=90, size = 2.25, nudge_x = -20) +
  labs(title="Average clinical incidence each year per age group with PMC (5 doses)", x="Age (days)", y="New clinical infections (ppy)",
       colour = "Level-1 admin")



######################### SAVE PMC IMPACT DATAFRAMES #################

# 8 dose PMC impact clinical cases 
write.csv(PMC_impact_ppy_8dose, paste0(getwd(), "/simulation_results/PMC_impact_ppy_8dose.csv"), row.names=FALSE)

# 8 dose PMC impact clinical cases 
write.csv(PMC_impact_ppy_5dose, paste0(getwd(), "/simulation_results/PMC_impact_ppy_5dose.csv"), row.names=FALSE)

