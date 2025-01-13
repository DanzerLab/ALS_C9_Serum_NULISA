##### HEAD #####################################################################
###                                                                          ###
### *Z6_Sex_Associations_DotPlots.R                                          ###
### *v1.0                                                                    ###
### *VG10.01.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: DotPlots of targets' associations with sex                 ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

Associations_Targets_Sex_Controls <- qread(
  "../Data/Models/Associations_Targets_Sex_Controls.qrds", 
  nthr=nthr
) 

Associations_Targets_Sex_Unsymptomatic <- qread(
  "../Data/Models/Associations_Targets_Sex_Unsymptomatic.qrds", 
  nthr=nthr
)

Associations_Targets_Sex_Controls_FC <- qread(
  "../Data/Models/Associations_Targets_Sex_Controls_FC.qrds", 
  nthr=nthr
) 

Associations_Targets_Sex_Unsymptomatic_FC <- qread(
  "../Data/Models/Associations_Targets_Sex_Unsymptomatic_FC.qrds", 
  nthr=nthr
)




  ### 2.0 Prepare data frames --------------------------------------------------

Cor_Controls <- data.frame(
  Target = names(Associations_Targets_Sex_Controls),
  FC=Associations_Targets_Sex_Controls_FC, 
  PVAL = sapply(Associations_Targets_Sex_Controls, function(x){return(x$p.value)})
) %>% 
  mutate(QVAL=p.adjust(PVAL, method="fdr")) %>% 
  mutate(SIGNIF=QVAL<0.05) %>% 
  mutate(DIRECTION=ifelse(SIGNIF, ifelse(FC<1, "blue", "red"), "grey")) %>% 
  arrange(FC) 
Cor_Controls$Target=factor(Cor_Controls$Target, levels=Cor_Controls$Target)
rownames(Cor_Controls) <- Cor_Controls$Target 


Cor_Unsymptomatic <- data.frame(
  Target = names(Associations_Targets_Sex_Unsymptomatic),
  FC=Associations_Targets_Sex_Unsymptomatic_FC, 
  PVAL = sapply(Associations_Targets_Sex_Unsymptomatic, function(x){return(x$p.value)})
) %>% 
  mutate(QVAL=p.adjust(PVAL, method="fdr")) %>% 
  mutate(SIGNIF=QVAL<0.05) %>% 
  mutate(DIRECTION=ifelse(SIGNIF, ifelse(FC<1, "blue", "red"), "grey")) %>% 
  arrange(FC) 
Cor_Unsymptomatic$Target=factor(Cor_Unsymptomatic$Target, levels=Cor_Unsymptomatic$Target)
rownames(Cor_Unsymptomatic) <- Cor_Unsymptomatic$Target  




  ### 3.0 Plot dotplots --------------------------------------------------------

Cor_Controls %>% 
  ggplot() + 
  aes(x=FC, y=Target, size=-log10(QVAL), fill=DIRECTION, alpha=-log10(QVAL)) + 
  geom_vline(xintercept = 1, col="#CCCCCC") + 
  geom_point(pch=21) + 
  scale_x_continuous(limits = c(-0.5, 2.5)) + 
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
  "../Data/Figures/SupplFig_3/DotPlot_Sex_Associations_Controls.pdf", 
  width = 6.20, 
  height = 15.30, 
  units = "in"
)


Cor_Unsymptomatic %>% 
  ggplot() + 
  aes(x=FC, y=Target, size=-log10(QVAL), fill=DIRECTION, alpha=-log10(QVAL)) + 
  geom_vline(xintercept = 1, col="#CCCCCC") + 
  geom_point(pch=21) + 
  scale_x_continuous(limits = c(-0.5, 2.5)) + 
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
  "../Data/Figures/SupplFig_3/DotPlot_Sex_Associations_Unsymptomatic.pdf", 
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


