##### HEAD #####################################################################
###                                                                          ###
### *Z2_SupplTab_Targets.R                                                   ###
### *v1.0                                                                    ###
### *VG09.01.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Suppl. Table Targets                                       ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

LOD <- qread(
  "../Data/NPQ/LOD.qrds", 
  nthr = nthr
)



  ### 3.0 Calculate NPQ - LOD matrix -------------------------------------------

LOD %>%
  write_excel_csv("../Data/Figures/SupplTab_2.csv", na=" - ", col_names = TRUE)


