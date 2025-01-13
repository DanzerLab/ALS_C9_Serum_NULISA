##### HEAD #####################################################################
###                                                                          ###
### *A3_Age_Models.R                                                         ###
### *v1.0                                                                    ###
### *VG09.01.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Generate age and sex linear models and adjusted data       ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)


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




  ### 2.0 Subset Controls ------------------------------------------------------

Control_Info <- SampleInfo %>% 
  filter(
    Case == "control"
  ) 

Control_NPQ <- NPQ %>% 
  select(all_of(Control_Info$SampleName)) %>% 
  t() %>% 
  as.data.frame()

setequal(Control_Info$SampleName, rownames(Control_NPQ))
all(
  rownames(Control_Info) == 
    rownames(Control_NPQ)
)




  ### 3.0 Subset Unsymptomatic cases -------------------------------------------

Unsymptomatic_Info <- SampleInfo %>% 
  filter(
    Symptomatic == "NonSymptomatic"
  ) 

Unsymptomatic_NPQ <- NPQ %>% 
  select(all_of(Unsymptomatic_Info$SampleName)) %>% 
  t() %>% 
  as.data.frame()

setequal(Unsymptomatic_Info$SampleName, rownames(Unsymptomatic_NPQ))
all(
  rownames(Unsymptomatic_Info) == 
    rownames(Unsymptomatic_NPQ)
)




  ### 4.0 Check Correlation with Age in Controls -------------------------------

Correlations_Targets_Age_Controls <- lapply(
  Control_NPQ, 
  function(x){
    cor.test(
      Control_Info$Age_at_collection, 
      x, 
      method = "pearson"
    )
  }
)

qsave(
  Correlations_Targets_Age_Controls, 
  "../Data/Models/Correlations_Targets_Age_Controls.qrds", 
  nthr=nthr
)




  ### 5.0 Check Correlation with Age in Unsymptomatic --------------------------

Correlations_Targets_Age_Unsymptomatic <- lapply(
  Unsymptomatic_NPQ, 
  function(x){
    cor.test(
      Unsymptomatic_Info$Age_at_collection, 
      x, 
      method = "pearson"
    )
  }
)

qsave(
  Correlations_Targets_Age_Unsymptomatic, 
  "../Data/Models/Correlations_Targets_Age_Unsymptomatic.qrds", 
  nthr=nthr
)




  ### 6.0 Check Association with Sex in Controls -------------------------------

Associations_Targets_Sex_Controls <- lapply(
  Control_NPQ, 
  function(x){
    wilcox.test(
      x[Control_Info$Sex=="f"], 
      x[Control_Info$Sex=="m"]
    )
  }
)

Associations_Targets_Sex_Controls_FC <- sapply(
  Control_NPQ, 
  function(x){
    mean(x[Control_Info$Sex=="f"])/mean(x[Control_Info$Sex=="m"])
  }
)

qsave(
  Associations_Targets_Sex_Controls, 
  "../Data/Models/Associations_Targets_Sex_Controls.qrds", 
  nthr=nthr
) 

qsave(
  Associations_Targets_Sex_Controls_FC, 
  "../Data/Models/Associations_Targets_Sex_Controls_FC.qrds", 
  nthr=nthr
) 




  ### 7.0 Check Association with Sex in Unsymptomatic --------------------------

Associations_Targets_Sex_Unsymptomatic <- lapply(
  Unsymptomatic_NPQ, 
  function(x){
    wilcox.test(
      x[Unsymptomatic_Info$Sex=="f"], 
      x[Unsymptomatic_Info$Sex=="m"]
    )
  }
)

Associations_Targets_Sex_Unsymptomatic_FC <- sapply(
  Unsymptomatic_NPQ, 
  function(x){
    mean(x[Unsymptomatic_Info$Sex=="f"])/mean(x[Unsymptomatic_Info$Sex=="m"])
  }
)


qsave(
  Associations_Targets_Sex_Unsymptomatic, 
  "../Data/Models/Associations_Targets_Sex_Unsymptomatic.qrds", 
  nthr=nthr
) 

qsave(
  Associations_Targets_Sex_Unsymptomatic_FC, 
  "../Data/Models/Associations_Targets_Sex_Unsymptomatic_FC.qrds", 
  nthr=nthr
) 




  ### 8.0 Generate models with corrected values --------------------------------



    ## 8.1 Build Matrix adjusted for Age and Sex from NonSymptomatic -----------

data <- NPQ %>% 
  t() %>% 
  as.data.frame() 

all(
  rownames(SampleInfo)==
    rownames(data)
)

data$Symptomatic <- SampleInfo$Symptomatic
data_NonSymptomatic <- subset(data, Symptomatic == "NonSymptomatic")
data_Symptomatic <- subset(data, Symptomatic == "Symptomatic")
data_NonSymptomatic <- data_NonSymptomatic %>% 
  select(-Symptomatic)
data_Symptomatic <- data_Symptomatic %>% 
  select(-Symptomatic)

SampleInfo_NonSymptomatic <- subset(SampleInfo, Symptomatic == "NonSymptomatic")
SampleInfo_Symptomatic <- subset(SampleInfo, Symptomatic == "Symptomatic")

all(
  rownames(SampleInfo_Symptomatic) == 
    rownames(data_Symptomatic)
)

all(
  rownames(SampleInfo_NonSymptomatic) == 
    rownames(data_NonSymptomatic)
)



models_Age_Sex_NonSymptomatic <- lapply(
  data_NonSymptomatic, 
  function(x){
    lm(x ~ SampleInfo_NonSymptomatic$Age_at_collection + SampleInfo_NonSymptomatic$Sex)
  }
)

Age_Sex_Corrected_NonSymptomatic <- data.frame(
  SampleName = SampleInfo_NonSymptomatic$SampleName
)
Age_Sex_Corrected_Symptomatic <- data.frame(
  SampleName = SampleInfo_Symptomatic$SampleName
)

for (i in 1:length(models_Age_Sex_NonSymptomatic)){
  df = data.frame(
    Age = SampleInfo_NonSymptomatic$Age_at_collection, 
    Sex = SampleInfo_NonSymptomatic$Sex, 
    Target = data_NonSymptomatic[,i]
  )
  mod.tmp <- lm(Target ~ Age + Sex, data = df)
  
  df = df %>% 
    select(Age, Sex)
  
  Age_Sex_Corrected_NonSymptomatic[[names(models_Age_Sex_NonSymptomatic)[i]]] <- predict(mod.tmp, newdata = df)
  
  df2 = data.frame(
    Age  = SampleInfo_Symptomatic$Age_at_collection, 
    Sex = SampleInfo_Symptomatic$Sex
  )
  
  Age_Sex_Corrected_Symptomatic[[names(models_Age_Sex_NonSymptomatic)[i]]] <- predict(mod.tmp, newdata = df2)
  rm(df, df2)
  
}

all(colnames(Age_Sex_Corrected_NonSymptomatic)==colnames(Age_Sex_Corrected_Symptomatic))
Age_Sex_Corrected <- rbind(
  Age_Sex_Corrected_NonSymptomatic, 
  Age_Sex_Corrected_Symptomatic
)

Age_Sex_Corrected <- Age_Sex_Corrected %>% 
  column_to_rownames("SampleName")

Age_Sex_Corrected <- Age_Sex_Corrected %>% 
  t() %>% 
  as.data.frame()

setequal(
  colnames(Age_Sex_Corrected), 
    SampleInfo$SampleName
)

Age_Sex_Corrected <- Age_Sex_Corrected[,match(SampleInfo$SampleName, colnames(Age_Sex_Corrected))]

all(
  colnames(Age_Sex_Corrected) == 
    SampleInfo$SampleName
)  

rm(data, data_NonSymptomatic, data_Symptomatic, SampleInfo_NonSymptomatic, SampleInfo_Symptomatic, mod.tmp, Age_Sex_Corrected_NonSymptomatic, Age_Sex_Corrected_Symptomatic, i)


qsave(
  models_Age_Sex_NonSymptomatic, 
  "../Data/Models/models_Age_Sex_NonSymptomatic.qrds", 
  nthr = nthr
)

qsave(
  Age_Sex_Corrected, 
  "../Data/Models/Age_Sex_Yhat.qrds", 
  nthr = nthr
)

cor.test(
  as.matrix(Age_Sex_Corrected), 
  as.matrix(NPQ)  
)



    ## 8.2 Build residuals matrix ----------------------------------------------

all(colnames(Age_Sex_Corrected)==colnames(NPQ))
all(rownames(Age_Sex_Corrected)==rownames(NPQ))
Residuals <- NPQ - Age_Sex_Corrected

qsave(
  Residuals, 
  "../Data/NPQ/Age_Sex_Corrected.qrds", 
  nthr = nthr
)
