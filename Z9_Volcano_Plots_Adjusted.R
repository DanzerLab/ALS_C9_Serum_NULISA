##### HEAD #####################################################################
###                                                                          ###
### *Z9_Volcano_Plots_Adjusted.R                                             ###
### *v1.0                                                                    ###
### *VG10.01.24                                                              ###
### *Project: NULISA_ALS_Serum                                               ###
###                                                                          ###
### *Description: Volcano plots of differential analysis with Age, Sex       ###
###                                                                          ###  
##### ENDHEAD ##################################################################




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)  
library(tidyverse)
library(ggrepel)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

nthr = 46




  ### 1.0 Load data ------------------------------------------------------------

res_Carrier_HC_Age_Sex <- qread(
  "../Data/Models/res_Carrier_HC_Age_Sex.qrds", 
  nthr = nthr
)


res_ALS_HC_Age_Sex <- qread(
  "../Data/Models/res_ALS_HC_Age_Sex.qrds", 
  nthr = nthr
)


res_ALS_Carrier_Age_Sex <- qread(
  "../Data/Models/res_ALS_Carrier_Age_Sex.qrds", 
  nthr = nthr
) 


res_Symptomatic_Nonsymptomatic_Age_Sex <- qread(
  "../Data/Models/res_Symptomatic_Nonsymptomatic_Age_Sex.qrds", 
  nthr = nthr
) 


res_C9_NonC9_Age_Sex <- qread(
  "../Data/Models/res_C9_NonC9_Age_Sex.qrds", 
  nthr = nthr
)





  ### 2.0 Defind AUX functions -------------------------------------------------

VolcanoPlot <- function(data){
  
  ggplot(df) + 
    aes(FC, LOG10QVAL, fill=DIRECTION, label = Target) + 
    geom_vline(xintercept = 0, col="#CCCCCC") + 
    geom_vline(xintercept = -1, col="#CCCCCC", linetype = "dashed") + 
    geom_vline(xintercept = 1, col="#CCCCCC", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "#CCCCCC", linetype = "dashed") + 
    geom_point(pch=21) + 
    geom_text_repel(data = df[df$SIGN,], size = 3, fontface = "bold.italic") + 
    scale_fill_manual(values=c("Grey"="#CCCCCC", "Blue"="#1133AA", "Red"="#DD1122")) + 
    xlab("\nFold change [ log2 ]") + 
    ylab("q-value [ -log10 ]\n") + 
    theme_classic() + 
    theme(
      axis.title = element_text(face = "bold.italic", colour = "#000000"), 
      axis.text = element_text(face = "bold.italic", colour = "#000000"), 
      legend.position = "Null"
    ) %>% 
    return()
  
}




  ### 3.0 Generate Volcano Plots -----------------------------------------------



    ## 3.1 Carrier vs. HC ------------------------------------------------------

res = res_Carrier_HC_Age_Sex

df <- data.frame(
  Target = res$modelStats$target, 
  FC = res$modelStats$CaseCarrier_coef,
  QVAL = res$modelStats$CaseCarrier_pval_FDR
) %>% 
  mutate(SIGN = QVAL < 0.05) %>% 
  mutate(DIRECTION = ifelse(SIGN, ifelse(FC < 0, "Blue", "Red"), "Grey")) %>% 
  mutate(LOG10QVAL = -log10(QVAL))


VolcanoPlot(df) + 
  scale_x_continuous(limits = c(-3,3)) + 
  scale_y_continuous(limits = c(0, 15))


ggsave(
  "../Data/Figures/Fig_1/VolcanoPlot_Car_HC_AgeSex.pdf", 
  width = 3.20, 
  height = 3.20, 
  units = "in"
) 


rm(res, df) 



    ## 3.2 ALS vs. HC ----------------------------------------------------------

res = res_ALS_HC_Age_Sex

df <- data.frame(
  Target = res$modelStats$target, 
  FC = res$modelStats$CasefALS_C9_coef,
  QVAL = res$modelStats$CasefALS_C9_pval_FDR
) %>% 
  mutate(SIGN = QVAL < 0.05) %>% 
  mutate(DIRECTION = ifelse(SIGN, ifelse(FC < 0, "Blue", "Red"), "Grey")) %>% 
  mutate(LOG10QVAL = -log10(QVAL)) 

VolcanoPlot(df) + 
  scale_x_continuous(limits = c(-3,3)) + 
  scale_y_continuous(limits = c(0, 15))


ggsave(
  "../Data/Figures/Fig_1/VolcanoPlot_ALS_HC_AgeSex.pdf", 
  width = 3.20, 
  height = 3.20, 
  units = "in"
) 


rm(res, df)



    ## 3.3 Symtomatic vs. NonSymptomatic ---------------------------------------

res = res_Symptomatic_Nonsymptomatic_Age_Sex

df <- data.frame(
  Target = res$modelStats$target, 
  FC = res$modelStats$SymptomaticSymptomatic_coef, 
  QVAL = res$modelStats$SymptomaticSymptomatic_pval_FDR
) %>% 
  mutate(SIGN = QVAL < 0.05) %>% 
  mutate(DIRECTION = ifelse(SIGN, ifelse(FC < 0, "Blue", "Red"), "Grey")) %>% 
  mutate(LOG10QVAL = -log10(QVAL)) 

VolcanoPlot(df) + 
  scale_x_continuous(limits = c(-3,3)) + 
  scale_y_continuous(limits = c(0, 15))


ggsave(
  "../Data/Figures/Fig_1/VolcanoPlot_Symptomatic_NonSymptomatic_AgeSex.pdf", 
  width = 3.20, 
  height = 3.20, 
  units = "in"
) 

rm(res, df)



    ## 3.4 C9ORF72 Positive vs. Negative ---------------------------------------

res = res_C9_NonC9_Age_Sex

df <- data.frame(
  Target = res$modelStats$target, 
  FC = res$modelStats$C9ORF72C9orf72_coef, 
  QVAL = res$modelStats$C9ORF72C9orf72_pval_FDR
) %>% 
  mutate(SIGN = QVAL < 0.05) %>% 
  mutate(DIRECTION = ifelse(SIGN, ifelse(FC < 0, "Blue", "Red"), "Grey")) %>% 
  mutate(LOG10QVAL = -log10(QVAL)) 

VolcanoPlot(df) + 
  scale_x_continuous(limits = c(-3,3)) + 
  scale_y_continuous(limits = c(0, 15)) 

ggsave(
  "../Data/Figures/Fig_1/VolcanoPlot_C9Positive_Negative_AgeSex.pdf", 
  width = 3.20, 
  height = 3.20, 
  units = "in"
) 

rm(res, df)

