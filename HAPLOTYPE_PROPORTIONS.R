library(pacman)
p_load(malariasimulation, foreSIGHT, malariaEquilibrium, tidyr, dplyr, ggplot2,
       reshape2, ggpubr, gridExtra, readxl, stringi, scene, XML, maps, readr,
       here, sf)

# set your working directory to where the haplotype data is found
setwd(paste0("C:/Users/IshanaAtkinson/OneDrive - London School of Hygiene and Tropical Medicine/Documents/malariasimulation work"))

# can add more countries when I clean excel data further
country_names <- c("Benin", "Mozambique", "Cameroon", "Cote d'Ivoire", "Zambia", "Ghana")

# set this to the name of your haplotype data file
haplotype_data <- read_xlsx("haplotypes_combined_forIshana2.xlsx")

# converts admin1 names into clean format that R can understand 
haplotype_data$NAME_2 <- stri_trans_general(str=gsub("-", "_", haplotype_data$NAME_1), id = "Latin-ASCII")
haplotype_data$NAME_2 <- stri_trans_general(str=gsub(" ", "_", haplotype_data$NAME_2), id = "Latin-ASCII")

haplotype_proportions <- data.frame()

I_AKA_ <- I_GKA_ <- I_GEA_ <- I_GEG_ <- V_GKA_ <- V_GKG_ <- c()

# loops over each country and admin1 area and calculates proportions for areas 
# which have 6 codon data and appends to data frame

for (i in 1:length(country_names)) {
  
  country_haplotype_data <- haplotype_data %>% filter(country == country_names[i])
  admin1 <- unique(country_haplotype_data$NAME_2)

  for (j in 1:length(admin1)) {
    
    # area haplotype data
    hap <- country_haplotype_data %>% filter(NAME_2 == admin1[j])
    
    denom <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS,
                              hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS,
                              hap$ISGEAA, hap$IAGEAA, hap$ISGEAS,
                              hap$ISGEGA, hap$IAGEGA, 
                              hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS,
                              hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA,
                              hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA
                              )), na.rm=TRUE)
    
    I_AKA_ <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    I_GKA_ <- sum(as.integer(c(hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    I_GEA_ <- sum(as.integer(c(hap$ISGEAA, hap$IAGEAA, hap$ISGEAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    I_GEG_ <- sum(as.integer(c(hap$ISGEGA, hap$IAGEGA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    V_GKA_ <- sum(as.integer(c(hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    V_GKG_ <- sum(as.integer(c(hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    other <- sum(as.integer(c(hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
    
    props_5sf <- round(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other), digits=5)
    sum_of_prop <- sum(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other))
    
    # check proportions add to 1 
    if (is.nan(sum_of_prop) == TRUE) {
    print(paste0("Summation of proportions for ", admin1[j], " in ", country_names[i], " cannot be calculated due to lack of 6-codon data"))
    } else {
      print(paste0("Summation of proportions for ", admin1[j], " in ", country_names[i], " is: ", sum_of_prop))
      
      haplotype_proportions <- rbind(haplotype_proportions, c(country_names[i], admin1[j], props_5sf, sum_of_prop))
    }
    
    
    
  }
}

# change column names 
colnames(haplotype_proportions) <- c("Country", "NAME_2", "I_AKA_", "I_GKA_", "I_GEA_", "I_GEG_", "V_GKA_", "V_GKG_", "other", "Sum of proportions" )


haplotype_proportions

