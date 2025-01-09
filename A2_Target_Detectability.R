##### HEAD #####################################################################
###                                                                          ###
### *A2_Target_Detectability.R                                               ###
### *v1.0                                                                    ###
### *VG12.12.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Calculate target detectability and filter targets          ###
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
  "../Data/NPQ/NPQ_Long.qrds"
)

NPQ$Target <- str_replace_all(
  NPQ$Target, 
  pattern="β", 
  replacement="Beta"
)

NPQ_Wide <- qread(
  "../Data/NPQ/NPQ_Wide.qrds", 
  nthr = nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

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

  ### 3.0 Calculate target detectability ---------------------------------------

LOD <- NPQ %>% 
  group_by(Target) %>% 
  summarize(LoD = only_value(LOD)) %>% 
  mutate(LoD = as.numeric(LoD))

summary(as.numeric(LOD$LoD))
LOD$LoD[is.na(LOD$LoD)] <- 0

setequal(
  rownames(NPQ_Wide), 
  LOD$Target
)
LOD <- LOD[match(rownames(NPQ_Wide), LOD$Target),]
all(
  rownames(NPQ_Wide) == 
    LOD$Target
)


LODs.tmp <- c()
for (i in 1:nrow(NPQ_Wide)){
  if(any(is.na(NPQ_Wide[i,]))) warning(paste0("Warning! Target ", rownames(NPQ_Wide)[i], " contains NA values"))
  LODs.tmp[i] <- mean(NPQ_Wide[i,] > LOD$LoD[i])*100 
  if(rownames(NPQ_Wide)[i] == LOD$Target[i]){
    names(LODs.tmp)[i] <- LOD$Target[i]
  } else{
    warning(paste0("Warning! NPQ_Wide and LOD Target names mismatch! "))
  }
}
all(names(LODs.tmp)==LOD$Target)
LOD$Detectability <- LODs.tmp
rm(LODs.tmp, i)

summary(LOD$Detectability < 50)
LOD$Target[LOD$Detectability < 50]
LOD$Target[LOD$Detectability < 10]

NPQ_Wide <- NPQ_Wide %>% 
  rownames_to_column("Target") %>% 
  filter(! Target %in% LOD$Target[LOD$Detectability < 10]) %>% 
  column_to_rownames("Target")


  ### 4.0 Save Data ------------------------------------------------------------

qsave(
  NPQ_Wide, 
  "../Data/NPQ/NPQ_Wide_Filtered.qrds", 
  nthr=nthr
)

qsave(
  LOD, 
  "../Data/NPQ/LOD.qrds", 
  nthr = nthr
)


