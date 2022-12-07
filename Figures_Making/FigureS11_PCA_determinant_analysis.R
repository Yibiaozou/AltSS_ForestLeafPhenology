rm(list=ls())
library(tidyverse)

load("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/AltSS_LP_RF_Caret_GMP_3GroupVar_3PC.RData")
load("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/AltSS_LP_RF_Caret_BiGrid_3GroupVar_3PC.RData")
load("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/AltSS_LP_RF_Caret_EvGrid_3GroupVar_3PC.RData")
load("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/AltSS_LP_RF_Caret_DeGrid_3GroupVar_3PC.RData")

# extract relative permutation importance
VarTmpDf_GMP <- as.data.frame(varImp(AltSS_LP_RF_Caret_GMP)$importance)
rownames(VarTmpDf_GMP)
VarTmpDf_GMP$VarNames <- rownames(VarTmpDf_GMP)
VarTmpDf_GMP$VarNames <- factor(VarTmpDf_GMP$VarNames, levels = rev(rownames(VarTmpDf_GMP))) # factorize the varaible names and fix the order

# extract relative permutation importance
VarTmpDf_BiGrid <- as.data.frame(varImp(AltSS_LP_RF_Caret_BiGrid)$importance)
rownames(VarTmpDf_BiGrid)
VarTmpDf_BiGrid$VarNames <- rownames(VarTmpDf_BiGrid)
VarTmpDf_BiGrid$VarNames <- factor(VarTmpDf_BiGrid$VarNames, levels = rev(rownames(VarTmpDf_BiGrid)))

# extract relative permutation importance
VarTmpDf_EvGrid <- as.data.frame(varImp(AltSS_LP_RF_Caret_EvGrid)$importance)
rownames(VarTmpDf_EvGrid)
VarTmpDf_EvGrid$VarNames <- rownames(VarTmpDf_EvGrid)
VarTmpDf_EvGrid$VarNames <- factor(VarTmpDf_EvGrid$VarNames, levels = rev(rownames(VarTmpDf_EvGrid)))

# extract relative permutation importance
VarTmpDf_DeGrid <- as.data.frame(varImp(AltSS_LP_RF_Caret_DeGrid)$importance)
rownames(VarTmpDf_DeGrid)
VarTmpDf_DeGrid$VarNames <- rownames(VarTmpDf_DeGrid)
VarTmpDf_DeGrid$VarNames <- factor(VarTmpDf_DeGrid$VarNames, levels = rev(rownames(VarTmpDf_DeGrid)))

## Figure S11A: global analysis
VarTmpDf_GMP$MID <- as.factor(1*VarTmpDf_GMP$Overall==max(VarTmpDf_GMP$Overall)) # highlight the most important variable with red color
ggplot(VarTmpDf_GMP, aes(x=VarNames, y=Overall, fill=MID))+ 
  geom_bar(stat="identity", position="dodge")+coord_flip()+
  scale_fill_manual( values = c( "grey30", "red"), guide = "none" )+
  ylab("Variable Importance")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_text(size = 12, face ="bold"))

## Figure S11B left: bimodal clusters
VarTmpDf_BiGrid$MID <- as.factor(1*VarTmpDf_BiGrid$Overall==max(VarTmpDf_BiGrid$Overall))
ggplot(VarTmpDf_BiGrid, aes(x=VarNames, y=Overall, fill=MID))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  scale_fill_manual( values = c( "grey30", "red"), guide = "none" )+
  ylab("Variable Importance")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_text(size = 12, face ="bold"))

## Figure S11B middle: evergreen-dominated clusters
VarTmpDf_EvGrid$MID <- as.factor(1*VarTmpDf_EvGrid$Overall==max(VarTmpDf_EvGrid$Overall))
ggplot(VarTmpDf_EvGrid, aes(x=VarNames, y=Overall, fill=MID))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  scale_fill_manual( values = c( "grey30", "red"), guide = "none" )+
  ylab("Variable Importance")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_text(size = 12, face ="bold"))

## Figure S11B right: deciduous-dominated clusters
VarTmpDf_DeGrid$MID <- as.factor(1*VarTmpDf_DeGrid$Overall==max(VarTmpDf_DeGrid$Overall))
ggplot(VarTmpDf_DeGrid, aes(x=VarNames, y=Overall, fill=MID))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  scale_fill_manual( values = c( "grey30", "red"), guide = "none" )+
  ylab("Variable Importance")+
  xlab("")+
  theme_classic()+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_text(size = 12, face ="bold"))