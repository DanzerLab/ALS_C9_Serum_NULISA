##### HEAD #####################################################################
###                                                                          ###
### *A4_Linear_Models.R                                                      ###
### *v1.0                                                                    ###
### *VG09.12.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Linear models for targeted differential proteomics         ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)
library(data.table)
library(NULISAseqR)
library(patchwork)

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

NPQ <- qread(
  "../Data/NPQ/NPQ_Wide_Filtered.qrds"
)

SampleInfo <- qread(
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)


SampleInfo$Case <- factor(SampleInfo$Case, levels = c("control", "Carrier", "fALS_C9")) 




  ### 2.0 Comparison Carrier - HC ----------------------------------------------


res_Carrier_HC_Age_Sex <- lmerNULISAseq(
  data = as.matrix(NPQ), 
  sampleInfo = SampleInfo, 
  sampleName_var = "SampleName", 
  modelFormula_fixed = "Case + Sex + Age_at_collection", 
  modelFormula_random = "(1|Collection_date)", 
  exclude_samples = c(
    SampleInfo$SampleName[SampleInfo$Disease_state=="symptomatic"]
  )
)

qsave(
  res_Carrier_HC_Age_Sex, 
  "../Data/Models/res_Carrier_HC_Age_Sex.qrds", 
  nthr=nthr
)




  ### 4.0 Comparison ALS - HC --------------------------------------------------

res_ALS_HC_Age_Sex <- lmerNULISAseq(
  data = as.matrix(NPQ), 
  sampleInfo = SampleInfo, 
  sampleName_var = "SampleName", 
  modelFormula_fixed = "Case + Sex + Age_at_collection", 
  modelFormula_random = "(1|Collection_date)", 
  exclude_samples = c(
    SampleInfo$SampleName[SampleInfo$Case=="Carrier"]
  )
)

qsave(
  res_ALS_HC_Age_Sex, 
  "../Data/Models/res_ALS_HC_Age_Sex.qrds", 
  nthr=nthr
)




  ### 5.0 Comparison ALS - Carrier ---------------------------------------------

res_ALS_Carrier_Age_Sex <- lmerNULISAseq(
  data = as.matrix(NPQ), 
  sampleInfo = SampleInfo, 
  sampleName_var = "SampleName", 
  modelFormula_fixed = "Case + Sex + Age_at_collection", 
  modelFormula_random = "(1|Collection_date)", 
  exclude_samples = c(
    SampleInfo$SampleName[SampleInfo$Case=="control"]
  )
) 

qsave(
  res_ALS_Carrier_Age_Sex, 
  "../Data/Models/res_ALS_Carrier_Age_Sex.qrds", 
  nthr = nthr
)




  ### 6.0 Comparison Symptomatic - NonSymptomatic (All) ------------------------

SampleInfo$Symptomatic <- factor(SampleInfo$Symptomatic, levels=c("NonSymptomatic", "Symptomatic"))

res_Symptomatic_Nonsymptomatic_Age_Sex <- lmerNULISAseq(
  data = as.matrix(NPQ), 
  sampleInfo = SampleInfo, 
  sampleName_var = "SampleName", 
  modelFormula_fixed = "Symptomatic + Sex + Age_at_collection", 
  modelFormula_random = "(1|Collection_date)", 
  exclude_samples = c()
)

qsave(
  res_Symptomatic_Nonsymptomatic_Age_Sex, 
  "../Data/Models/res_Symptomatic_Nonsymptomatic_Age_Sex.qrds", 
  nthr = nthr
)




  ### 7.0 Comparison C9orf72-HRE Positive vs. Negative -------------------------

SampleInfo$C9ORF72 <- factor(SampleInfo$Mutation, levels = c("no", "C9orf72")) 

res_C9_NonC9_Age_Sex <- lmerNULISAseq(
  data = as.matrix(NPQ), 
  sampleInfo = SampleInfo, 
  sampleName_var = "SampleName", 
  modelFormula_fixed = "C9ORF72 + Sex + Age_at_collection", 
  modelFormula_random = "(1|Collection_date)", 
  exclude_samples = c()
)

qsave(
  res_C9_NonC9_Age_Sex, 
  "../Data/Models/res_C9_NonC9_Age_Sex.qrds", 
  nthr = nthr
)












