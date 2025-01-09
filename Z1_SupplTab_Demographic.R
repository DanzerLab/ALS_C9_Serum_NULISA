##### HEAD #####################################################################
###                                                                          ###
### *Z1_SupplTab_Demographic.R                                               ###
### *v1.0                                                                    ###
### *VG20.12.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Boxplots of Target NPQs with LOD                           ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

SampleInfo <- qread(
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)




  ### 2.0 AUX function definition ----------------------------------------------

only_value <- function(x){
  if (any(is.null(x)) | length(x)==0) stop("Input not valid") 
  
  if(length(unique(x)==1)){
    return(x[1])
  } else {
    if(length(unique(x)[!is.na(unique(x))])==1){
      warning("The input contains NAs! ")
      return(unique(x)[!is.na(unique(x))][1])
    }
    stop("The input contains more than one unique value!")
  }
}




  ### 3.0 Calculate NPQ - LOD matrix -------------------------------------------

SampleInfo %>% 
  group_by(Case) %>% 
  mutate(Disease_progression = Disease_progression/12) %>% 
  summarize(
    Age_mean = Age_at_collection %>% mean %>% round(1) %>% format(nsmall=1), 
    Age_sd = Age_at_collection %>% sd %>% round(1) %>% format(nsmall=1), 
    Sex_m = (Sex == "m") %>% sum,
    Sex_f = (Sex == "f") %>% sum, 
    Age_at_onset_mean = Age_at_onset %>% mean %>% round(1) %>% format(nsmall=1), 
    Age_at_onset_sd = Age_at_onset %>% sd %>% round(1) %>% format(nsmall=1), 
    Age_at_onset_min = Age_at_onset %>% min %>% round(1) %>% format(nsmall=1), 
    Age_at_onset_max = Age_at_onset %>% max %>% round(1) %>% format(nsmall=1), 
    ALSFRS_mean = ALSFRS %>% mean %>% round(1) %>% format(nsmall=1), 
    ALSFRS_sd = ALSFRS %>% sd %>% round(1) %>% format(nsmall=1), 
    ALSFRS_min = ALSFRS %>% min %>% round(1) %>% format(nsmall=1), 
    ALSFRS_max = ALSFRS %>% max %>% round(1) %>% format(nsmall=1), 
    Disease_duration_mean = Disease_duration %>% mean %>% round(1) %>% format(nsmall=1), 
    Disease_duration_sd = Disease_duration %>% sd %>% round(1) %>% format(nsmall=1), 
    Disease_duration_min = Disease_duration %>% min %>% round(1) %>% format(nsmall=1), 
    Disease_duration_max = Disease_duration %>% max %>% round(1) %>% format(nsmall=1), 
    Disease_progression_mean = Disease_progression %>% mean %>% round(1) %>% format(nsmall=1), 
    Disease_progression_sd = Disease_progression %>% sd %>% round(1) %>% format(nsmall=1), 
    Disease_progression_min = Disease_progression %>% min %>% round(1) %>% format(nsmall=1), 
    Disease_progression_max = Disease_progression %>% max %>% round(1) %>% format(nsmall=1)
  ) %>% 
  mutate(
    Age = paste0(Age_mean, "±", Age_sd), 
    Sex_f_m = paste0(" ", Sex_f, "/", Sex_m), 
    Age_at_onset = paste0(
      Age_at_onset_mean, "±", Age_at_onset_sd, 
      "\n", 
      "(", Age_at_onset_min, " - ", Age_at_onset_max, ")"
    ),
    Disease_duration = paste0(
      Disease_duration_mean, "±", Disease_duration_sd, 
      "\n", 
      "(", Disease_duration_min, " - ", Disease_duration_max, ")"
    ),
    Disease_Progression_Rate = paste0(
      Disease_progression_mean, "±", Disease_progression_sd, 
      "\n", 
      "(", Disease_progression_min, " - ", Disease_progression_max, ")"
    ), 
    Disease_severity = paste0(
      ALSFRS_mean, "±", ALSFRS_sd, 
      "\n", 
      "(", ALSFRS_min, " - ", ALSFRS_max, ")"
    )
   
  ) %>% 
  select(Case, Age, Sex_f_m, Age_at_onset, Disease_duration, Disease_Progression_Rate, Disease_severity) %>% 
  apply(FUN=str_replace_all, MARGIN=1, pattern="NA±NA", repl="") %>% 
  apply(FUN=str_replace_all, MARGIN=2, pattern="\\(NA - NA\\)", repl="") %>% 
  as.data.frame() %>% 
  mutate(Val=c(
    "", 
    "Age [y] \nmean±sd \n(range)", 
    "Sex [n] \nf/m", 
    "Age at onset [y] \nmean±sd \n(range)", 
    "Disease duration [y] \nmean±sd \n(range)", 
    "Disease progression rate [ΔALS-FRS/m] \nmean±sd \n(range)", 
    "Disase severity [ALS-FRS] \nmean±sd \n(range)"
  )) %>% 
  select(c(4,2,1,3)) %>%
  write_excel_csv("../Data/Figures/SupplTab_1.csv", na=" - ", col_names = FALSE)





