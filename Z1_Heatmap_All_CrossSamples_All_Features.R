##### HEAD #####################################################################
###                                                                          ###
### *Z1_Heatmap_All_CrossSamples_All_Features                                ###
### *v1.0                                                                    ###
### *VG05.12.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Heatmap of all cross-sectional samples & all features      ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)
library(pheatmap)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

NPQ_Wide <- qread(
  "../Data/NPQ/NPQ_Wide.qrds"
)

SampleInfo <- qread(
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)

Control_Samples <- c(
  "SCRep01", 
  "SCRep02", 
  "SCRep03"
)




  ### 2.0 Define AUX functions -------------------------------------------------

print_dims <- function(x){
  tryCatch({print(message(paste0("Dimensions: ", paste0(dim(x), collapse = ", "))))}, error = function(e){message("Dimentions could not be determined...")}) 
  return(x)
}




  ### 3.0 Generate heatmap -----------------------------------------------------

mat <- NPQ_Wide %>% 
  select(
    all_of(
      SampleInfo %>% 
        filter(LatestTimepoint) %>% 
        filter(Include) %>% 
        pull(SampleName)
    )
  ) %>% 
  print_dims %>% 
  t() %>% scale(center=TRUE, scale=TRUE) %>% t() 


col <- SampleInfo %>% 
  filter(LatestTimepoint) %>% 
  select(SampleName, Case) %>% 
  column_to_rownames("SampleName") 


p1 <- pheatmap(
  mat = mat, 
  annotation_col = col, 
  annotation_colors = list(Case=c(Carrier="#5a92f2", control="purple", fALS_C9= "darkolivegreen2")),
  scale = "none", 
  cluster_cols = hclust(
    dist(t(mat), method="euclidean", p=2), 
    method = "ward.D2"
  ),
  
  color = colorRampPalette(c( "blue", "blue", "white", "red", "red"))(100),
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  treeheight_row = 0
) 



col_dend <- p1$tree_col
col_dend <- dendextend::rotate(col_dend, order=c(1:37, 61:49,38:42, 44:48,43,63,77:76, 74:75, 73:64,62))

pheatmap::pheatmap(
  mat = mat, 
  annotation_col = col, 
  annotation_colors = list(Case=c(Carrier="#5a92f2", control="darkolivegreen2", fALS_C9= "#b31405")),
  scale = "none", , 
  cluster_cols = col_dend,
  color = colorRampPalette(c( "blue", "blue", "white", "red", "red"))(100),
  show_rownames = FALSE, 
  show_colnames = FALSE, 
  treeheight_row = 0
  
)




