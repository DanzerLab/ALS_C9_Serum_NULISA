##### HEAD #####################################################################
###                                                                          ###
### *Z4_Heatmap_AllSamples_AllFeatures_Unadjusted.R                          ###
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

NPQ <- qread(
  "../Data/NPQ/NPQ_Wide_Filtered.qrds", 
  nthr = nthr
)

SampleInfo <- qread(
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

print_dims <- function(x){
  tryCatch({print(message(paste0("Dimensions: ", paste0(dim(x), collapse = ", "))))}, error = function(e){message("Dimentions could not be determined...")}) 
  return(x)
}




  ### 3.0 Generate heatmap -----------------------------------------------------

mat <- NPQ %>% 
  print_dims %>% 
  t() %>% scale(center=TRUE, scale=TRUE) %>% t() 


col <- SampleInfo %>% 
  select(Case, Age_at_collection, Sex) 

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
col_dend <- dendextend::rotate(col_dend, order=c(1:56, 58:67, 74:73, 71:72, 70, 75:77, 69:68, 57))

pheatmap::pheatmap(
  mat = mat, 
  annotation_col = col, 
  annotation_colors = list(
    Case=c(Carrier="#5a92f2", control="darkolivegreen2", fALS_C9= "#b31405"), 
    Sex = c(f="#f7d9ba", m="#7d6cb8")
  ),
  scale = "none", , 
  cluster_cols = col_dend,
  color = colorRampPalette(c( "blue", "blue", "white", "red", "red"))(100),
  show_rownames = FALSE, 
  show_colnames = FALSE, 
  treeheight_row = 0
  
)

dev.copy2pdf(
  file = "../Data/Figures/SupplFig_2/Heatmap_AllCases_AllTargets_Unadjusted.pdf", 
  width = 12.17, 
  height = 6.23
); dev.off()



