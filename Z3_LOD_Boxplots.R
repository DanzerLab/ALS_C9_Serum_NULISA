##### HEAD #####################################################################
###                                                                          ###
### *Z3_LOD_Boxplots.R                                                       ###
### *v1.0                                                                    ###
### *VG12.12.24                                                              ###
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

NPQ_Wide <- qread(
  "../Data/NPQ/NPQ_Wide.qrds", 
  nthr=nthr
)

SampleInfo <- qread(
  "../Data/Samples/SampleInfo.qrds", 
  nthr = nthr
)

LOD <- qread(
  "../Data/NPQ/LOD.qrds", 
  nthr=nthr
)




  ### 2.0 Calculate NPQ - LOD matrix -------------------------------------------

all(rownames(NPQ_Wide)==LOD$Target)
NPQ_tall <- as.data.frame(t(NPQ_Wide))
all(colnames(NPQ_tall) == LOD$Target)

NPQ_delLOD <- sweep(NPQ_tall, MARGIN = 2, STATS = LOD$LoD, FUN="-")
NPQ_delLOD <- NPQ_delLOD %>% select(
  apply(NPQ_delLOD, MARGIN = 2, FUN=median) %>% 
    as.numeric() %>% order(decreasing = TRUE)
)
rm(NPQ_tall)




  ### 3.0 Generate Boxplots ----------------------------------------------------

NPQ_delLOD$Case <- SampleInfo$Case[match(rownames(NPQ_delLOD), SampleInfo$SampleName)]
df <- NPQ_delLOD %>% 
  pivot_longer(
    cols = - Case, 
    names_to = "Target", 
    values_to = "NPQ"
  ) %>% 
  mutate(Target = factor(Target, levels = unique(Target))) %>% 
  mutate(Detectability = LOD$Detectability[match(Target, LOD$Target)]) 
  
boxplot_stats <- df %>%
  group_by(Target) %>%
  summarise(
    lower = quantile(NPQ, 0.25) - 1.5 * IQR(NPQ),
    upper = quantile(NPQ, 0.75) + 1.5 * IQR(NPQ),
    lower_whisker = quantile(NPQ, 0.25) - 3 * IQR(NPQ),
    upper_whisker = quantile(NPQ, 0.75) + 3 * IQR(NPQ)
  )

outliers <- df %>%
  left_join(boxplot_stats, by = "Target") %>%
  filter(NPQ < lower | NPQ > upper)

ggplot(df) + 
  aes(Target, NPQ) + 
  geom_hline(yintercept = 0, col="orangered3")  +
  geom_boxplot(aes(fill=Detectability), col="#000000", fatten=2, outlier.shape = NA) + 
  scale_fill_viridis_c(name = "Detectability", option="mako") +  
  ggnewscale::new_scale("fill") + 
  geom_point(data = outliers, aes(x = Target, y = NPQ, fill = Case), 
             pch = 21, color = "black") + 
  scale_fill_manual(values=c(Carrier="#5a92f2", control="darkolivegreen2", fALS_C9= "#b31405")) + 
  ylab("NPQ - LoD\n") + 
  theme_classic() + 
  theme(
    panel.grid.major.y = element_line(color="#CCCCCC", linetype=2), 
    axis.line = element_line("#000000"),  
    axis.title.x = element_blank(), 
    axis.text.x = element_text(size=8, face = "bold", color="#000000", angle = 90, hjust = 1.0), 
    axis.title.y = element_text(size=10, color="#000000", face = "bold.italic"), 
    axis.text.y = element_text(size=9, color="#000000", face= "bold.italic")
  )

ggsave(
  "../Data/Figures/SupplFig_LoD/NPQ_LoD_Boxplots.pdf", 
  width = 18.20, 
  height = 10.30, 
  units = "in"
)


