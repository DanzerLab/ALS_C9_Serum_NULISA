##### HEAD #####################################################################
###                                                                          ###
### *A1_Read_Data.R                                                          ###
### *v1.0                                                                    ###
### *VG12.12.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Read-in and formatting of the NULISA data                  ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse) 
library(lubridate)
library(data.table)
library(NULISAseqR)

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 44




  ### 1.0 Load data ------------------------------------------------------------



    ## 1.1 Intra-plate normalized, log2+1-transformed counts -------------------

NPQ <- readxl::read_xlsx(
  path = paste0(
    "../Data/Input/", 
    "Alamar_Input_data.xlsx"
  ), 
  sheet = readxl::excel_sheets(
    paste0(
      "../Data/Input/", 
      "Alamar_Input_data.xlsx"
    )
  )[1]
)

NPQ$Target <- str_replace_all(
  NPQ$Target, 
  pattern="Î²", 
  replacement="β"
)

P000572_Target_List <- fread(
  paste0(
    "../Data/Input/", 
    "P_000572_Target_List",
    ".txt"
  ), 
  data.table = FALSE
)

colnames(P000572_Target_List) <- str_replace_all(
  colnames(P000572_Target_List), 
  pattern = " ", 
  replacement = "_"
)

colnames(P000572_Target_List) <- str_replace_all(
  colnames(P000572_Target_List), 
  pattern = "/", 
  replacement = "_"
)




    ## 1.2 Sample Info ---------------------------------------------------------

SampleInfo <- fread(
  paste0(
    "../Data/Input/", 
    "Metafile_NULISA.csv"
  ), 
  data.table=FALSE
)

SampleInfo$Sample_ID <- str_replace_all(
  SampleInfo$Sample_ID, 
  "#", 
  ""
)
SampleInfo$SampleName <- paste0(
  "Sample", 
  SampleInfo$Sample_ID
)

SampleInfo$Symptomatic <- "Symptomatic"
SampleInfo$Symptomatic[SampleInfo$Case %in% c("Carrier", "control")] <- "NonSymptomatic"




  ### 2.0 Sanity checks --------------------------------------------------------

      

    ## 2.1 Targets and UniProt IDs ---------------------------------------------

unique(NPQ$UniProtID)[
  ! unique(NPQ$UniProtID) %in% unique(P000572_Target_List$Uniprot_ID)
] 

unique(P000572_Target_List$Uniprot_ID)[
  ! unique(P000572_Target_List$Uniprot_ID) %in% unique(NPQ$UniProtID) 
] 

NPQ$Target[
  NPQ$UniProtID == unique(NPQ$UniProtID)[
    ! unique(NPQ$UniProtID) %in% unique(P000572_Target_List$Uniprot_ID)
  ] 
] %>% unique() 



P000572_Target_List$Gene[
  P000572_Target_List$Uniprot_ID== unique(P000572_Target_List$Uniprot_ID)[
    ! unique(P000572_Target_List$Uniprot_ID) %in% unique(NPQ$UniProtID) 
  ] 
]

      # The only UniProt mismatch is a change in the IGF1R UniProt ID -> OK 

table(
  unique(NPQ$Target) %in% 
    unique(P000572_Target_List$Gene)
)

table(
  unique(P000572_Target_List$Gene) %in% 
    unique(NPQ$Target)
)

unique(NPQ$Target)[
  ! unique(NPQ$Target) %in% unique(P000572_Target_List$Gene)
]

unique(P000572_Target_List$Gene)[
  ! unique(P000572_Target_List$Gene) %in% unique(NPQ$Target)
]

     # The mismatching Targets are synonymous -> OK 




  ### 3.0 Format data ----------------------------------------------------------

NPQ_Wide <- reshape(
  as.data.frame(NPQ[,c(4,6,12)]), 
  direction = "wide", 
  idvar = "Target", 
  timevar = "SampleName"
)
rownames(NPQ_Wide) <- str_replace_all(NPQ_Wide$Target, "β", "Beta")

NPQ_Wide <- NPQ_Wide %>% 
  select(-Target)

colnames(NPQ_Wide) <- apply(
  str_split(
    colnames(NPQ_Wide), 
    "_", 
    simplify = TRUE
  )[,c(3,4)], 
  FUN=function(x){paste0(x, collapse="")}, 
  MARGIN = 1
)


sum(is.na(SampleInfo$Collection_date) | SampleInfo$Collection_date=="") == 
  sum(is.na(dmy(SampleInfo$Collection_date))) 

SampleInfo <- SampleInfo %>% 
  group_by(Case_ID) %>% 
  mutate(LatestTimepoint = dmy(Collection_date) == max(dmy(Collection_date), na.rm = TRUE)) %>% 
  arrange(as.numeric(Sample_ID)) %>% 
  ungroup() %>% 
  mutate(LatestTimepoint = replace(LatestTimepoint, is.na(LatestTimepoint), TRUE))


SampleInfo_AllSamples <- SampleInfo 
NPQ_Wide_AllSamples <- NPQ_Wide

SampleInfo <- SampleInfo %>% 
  filter(Include & LatestTimepoint) %>% 
  select(- c(Include, LatestTimepoint)) %>% 
  column_to_rownames("SampleName") %>% 
  mutate(SampleName = rownames(.))

NPQ_Wide <- NPQ_Wide %>% 
  select(all_of(rownames(SampleInfo)))




  ### 4.0 Export data ----------------------------------------------------------

qsave(
  NPQ, 
  "../Data/NPQ/NPQ_Long.qrds", 
  nthr = nthr
)

qsave(
  NPQ_Wide_AllSamples, 
  "../Data/NPQ/NPQ_Wide_AllSamples.qrds", 
  nthr = nthr
)

qsave(
  NPQ_Wide, 
  "../Data/NPQ/NPQ_Wide.qrds", 
  nthr = nthr
)


qsave(
  SampleInfo_AllSamples, 
  "../Data/Samples/SampleInfo_AllSamples.qrds", 
  nthr = nthr
)

qsave(
  SampleInfo, 
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)
