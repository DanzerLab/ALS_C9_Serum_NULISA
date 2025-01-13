##### HEAD #####################################################################
###                                                                          ###
### *Z5_Age_Correlations_DotPlots.R                                          ###
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

Correlations_Targets_Age_Controls <- qread(
  "../Data/Models/Correlations_Targets_Age_Controls.qrds", 
  nthr=nthr
)

Correlations_Targets_Age_Unsymptomatic <- qread(
  "../Data/Models/Correlations_Targets_Age_Unsymptomatic.qrds", 
  nthr=nthr
)




  ### 2.0 Prepare data frames --------------------------------------------------

Cor_Controls <- data.frame(
  Target = names(Correlations_Targets_Age_Controls),
  R=sapply(Correlations_Targets_Age_Controls, function(x){return(x$estimate)}), 
  PVAL = sapply(Correlations_Targets_Age_Controls, function(x){return(x$p.value)})
) %>% 
  mutate(QVAL=p.adjust(PVAL, method="fdr")) %>% 
  mutate(SIGNIF=QVAL<0.05) %>% 
  mutate(DIRECTION=ifelse(SIGNIF, ifelse(R<0, "blue", "red"), "grey")) %>% 
  arrange(R) 
Cor_Controls$Target=factor(Cor_Controls$Target, levels=Cor_Controls$Target)
rownames(Cor_Controls) <- Cor_Controls$Target 


Cor_Unsymptomatic <- data.frame(
  Target = names(Correlations_Targets_Age_Unsymptomatic),
  R=sapply(Correlations_Targets_Age_Unsymptomatic, function(x){return(x$estimate)}), 
  PVAL = sapply(Correlations_Targets_Age_Unsymptomatic, function(x){return(x$p.value)})
) %>% 
  mutate(QVAL=p.adjust(PVAL, method="fdr")) %>% 
  mutate(SIGNIF=QVAL<0.05) %>% 
  mutate(DIRECTION=ifelse(SIGNIF, ifelse(R<0, "blue", "red"), "grey")) %>% 
  arrange(R) 
Cor_Unsymptomatic$Target=factor(Cor_Unsymptomatic$Target, levels=Cor_Unsymptomatic$Target)
rownames(Cor_Unsymptomatic) <- Cor_Unsymptomatic$Target  




  ### 3.0 Plot dotplots --------------------------------------------------------

Cor_Controls %>% 
  ggplot() + 
  aes(x=R, y=Target, size=-log10(QVAL), fill=DIRECTION, alpha=-log10(QVAL)) + 
  geom_vline(xintercept = 0, col="#CCCCCC") + 
  geom_point(pch=21) + 
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_discrete(expand=c(0.02,0.02)) + 
  scale_fill_manual(values=c("grey"="grey", "blue"="blue", "red"="red")) + 
  scale_size_continuous(limits=c(0,5)) + 
  scale_alpha_continuous(limits=c(0,5)) + 
  xlab("Pearson correlation coefficient (r)") + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(face="bold.italic", colour="#000000"), 
    axis.title.y = element_blank(), 
    axis.text = element_text(face="bold.italic", colour="#000000")
  )

ggsave(
  "../Data/Figures/SupplFig_3/DotPlot_Age_Correlations_Controls.pdf", 
  width = 6.20, 
  height = 15.30, 
  units = "in"
)


Cor_Unsymptomatic %>% 
  ggplot() + 
  aes(x=R, y=Target, size=-log10(QVAL), fill=DIRECTION, alpha=-log10(QVAL)) + 
  geom_vline(xintercept = 0, col="#CCCCCC") + 
  geom_point(pch=21) + 
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_discrete(expand=c(0.02,0.02)) + 
  scale_fill_manual(values=c("grey"="grey", "blue"="blue", "red"="red")) + 
  scale_size_continuous(limits=c(0,5)) + 
  scale_alpha_continuous(limits=c(0,5)) + 
  xlab("Pearson correlation coefficient (r)") + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(face="bold.italic", colour="#000000"), 
    axis.title.y = element_blank(), 
    axis.text = element_text(face="bold.italic", colour="#000000")
  ) 

ggsave(
  "../Data/Figures/SupplFig_3/DotPlot_Age_Correlations_Unsymptomatic.pdf", 
  width = 6.20, 
  height = 15.30, 
  units = "in"
)

Cor_Controls %>% 
  filter(QVAL<0.05) %>% 
  select(Target)

Cor_Unsymptomatic %>% 
  filter(QVAL<0.05) %>% 
  select(Target)

table(
  {Cor_Controls %>% 
    filter(QVAL<0.05) %>% 
    pull(Target) %>% 
    as.character()} %in% 
  {Cor_Unsymptomatic%>% 
    filter(QVAL<0.05) %>% 
    pull(Target) %>% 
    as.character()}
)

setdiff(
  Cor_Controls %>% 
    filter(QVAL<0.05) %>% 
    pull(Target) %>% 
    as.character(), 
  Cor_Unsymptomatic%>% 
      filter(QVAL<0.05) %>% 
      pull(Target) %>% 
      as.character()
)

table(
  {Cor_Unsymptomatic%>% 
      filter(QVAL<0.05) %>% 
      pull(Target) %>% 
      as.character()} %in% 
  {Cor_Controls %>% 
      filter(QVAL<0.05) %>% 
      pull(Target) %>% 
      as.character()}
)

setdiff(
  Cor_Unsymptomatic%>% 
    filter(QVAL<0.05) %>% 
    pull(Target) %>% 
    as.character(), 
  Cor_Controls %>% 
    filter(QVAL<0.05) %>% 
    pull(Target) %>% 
    as.character()
)
