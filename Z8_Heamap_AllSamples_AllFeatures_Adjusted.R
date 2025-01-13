##### HEAD #####################################################################
###                                                                          ###
### *Z8_Heamap_AllSamples_AllFeatures_Adjusted.R                             ###
### *v1.0                                                                    ###
### *VG10.01.24                                                              ###
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

Residuals <- qread(
  "../Data/NPQ/Age_Sex_Corrected.qrds", 
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

mat <- Residuals %>% 
  print_dims %>% 
  t() %>% scale(center=TRUE, scale=TRUE) %>% t() 


col <- SampleInfo %>% 
  select(Case) 

p1 <- pheatmap(
  mat = mat, 
  annotation_col = col, 
  annotation_colors = list(Case=c(Carrier="#5a92f2", control="purple", fALS_C9= "darkolivegreen2")),
  scale = "none", 
  cluster_cols = hclust(
    dist(t(mat), method="euclidean"), 
    method = "ward.D2"
  ),
  
  color = colorRampPalette(c( "blue", "blue", "white", "red", "red"))(100),
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  treeheight_row = 0
) 



col_dend <- p1$tree_col
col_dend <- dendextend::rotate(col_dend, order=c(46:26, 1:2, 25:13, 3:7, 10:12, 9:8, 64:68, 70:77, 74, 69, 63:47))

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
  file = "../Data/Figures/Fig_1/Heatmap_AllCases_AllTargets_Adjusted_For_Age_Sex.pdf", 
  width = 12.17, 
  height = 6.23
); dev.off()


