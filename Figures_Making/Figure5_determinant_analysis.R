rm(list=ls())
library(tidyverse)
library(tictoc)
library(data.table)
library(plotrix)

## load dataset which contain information of variable importance from fitted random forest models
varImp_GMP_Df <- readRDS("**/Global_Analysis/Model_Output/varImp_GMP_GFBI_Df_FSD.rds")
varImp_BiGrid_Df <- readRDS("**/Global_Analysis/Model_Output/varImp_BiGrid_GFBI_Df_FSD.rds")
varImp_EvGrid_Df <- readRDS("**/Global_Analysis/Model_Output/varImp_EvGrid_GFBI_Df_FSD.rds")
varImp_DeGrid_Df <- readRDS("**/Global_Analysis/Model_Output/varImp_DeGrid_GFBI_Df_FSD.rds")

## expand the short names to the full variable names for all four dataset
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="MAT"] <- "Mean annual temperature"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="MAP"] <- "Annual precipitation"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="Tco"] <- "Temperature coldest quarter"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="Pdr"] <- "Precipitation driest quarter"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="Soil_pH"] <- "Soil pH 0 to 100cm"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="Soil_N"] <- "Soil Nitrogen"
varImp_GMP_Df$VarName[varImp_GMP_Df$VarName=="Soil_CNratio"] <- "Soil C:N ratio"

varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="MAT"] <- "Mean annual temperature"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="MAP"] <- "Annual precipitation"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="Tco"] <- "Temperature coldest quarter"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="Pdr"] <- "Precipitation driest quarter"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="Soil_pH"] <- "Soil pH 0 to 100cm"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="Soil_N"] <- "Soil Nitrogen"
varImp_BiGrid_Df$VarName[varImp_BiGrid_Df$VarName=="Soil_CNratio"] <- "Soil C:N ratio"

varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="MAT"] <- "Mean annual temperature"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="MAP"] <- "Annual precipitation"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="Tco"] <- "Temperature coldest quarter"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="Pdr"] <- "Precipitation driest quarter"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="Soil_pH"] <- "Soil pH 0 to 100cm"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="Soil_N"] <- "Soil Nitrogen"
varImp_EvGrid_Df$VarName[varImp_EvGrid_Df$VarName=="Soil_CNratio"] <- "Soil C:N ratio"

varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="MAT"] <- "Mean annual temperature"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="MAP"] <- "Annual precipitation"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="Tco"] <- "Temperature coldest quarter"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="Pdr"] <- "Precipitation driest quarter"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="Soil_pH"] <- "Soil pH 0 to 100cm"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="Soil_N"] <- "Soil Nitrogen"
varImp_DeGrid_Df$VarName[varImp_DeGrid_Df$VarName=="Soil_CNratio"] <- "Soil C:N ratio"

## make sure the values of variable importance are numeric
varImp_GMP_Df$VarImp <- as.numeric(varImp_GMP_Df$VarImp)
varImp_BiGrid_Df$VarImp <- as.numeric(varImp_BiGrid_Df$VarImp)
varImp_EvGrid_Df$VarImp <- as.numeric(varImp_EvGrid_Df$VarImp)
varImp_DeGrid_Df$VarImp <- as.numeric(varImp_DeGrid_Df$VarImp)

# aggregate the dataset to calculate mean, standard deviation and standard error of variable importance across bootstrapping repetition
varImp_GMP_Df_scaled <- varImp_GMP_Df%>%group_by(VarName)%>%summarize(VarImp_Mean=mean(VarImp), # mean
                                                                      VarImp_SD=sd(VarImp), # standard deviation
                                                                      VarImp_SE=std.error(VarImp)) # standard error
varImp_GMP_Df_scaled$order <- c(4, 7, 5, 3, 1, 2, 6) # to make sure when plotting, climatic variables locate together, whereas soil variables locate together
varImp_GMP_Df_scaled$VarName <- fct_reorder(varImp_GMP_Df_scaled$VarName, varImp_GMP_Df_scaled$order, .desc = F) # factorize variable order

# do the same for the rest three dataset
varImp_BiGrid_Df_scaled <- varImp_BiGrid_Df%>%group_by(VarName)%>%summarize(VarImp_Mean=mean(VarImp),
                                                                            VarImp_SD=sd(VarImp), 
                                                                            VarImp_SE=std.error(VarImp))
varImp_BiGrid_Df_scaled$VarName <- factor(varImp_BiGrid_Df_scaled$VarName, levels = levels(varImp_GMP_Df_scaled$VarName)) # use the same order of variable as 'varImp_GMP_Df_scaled'


varImp_EvGrid_Df_scaled <- varImp_EvGrid_Df%>%group_by(VarName)%>%summarize(VarImp_Mean=mean(VarImp),
                                                                            VarImp_SD=sd(VarImp), 
                                                                            VarImp_SE=std.error(VarImp))
varImp_EvGrid_Df_scaled$VarName <- factor(varImp_EvGrid_Df_scaled$VarName, levels = levels(varImp_GMP_Df_scaled$VarName))


varImp_DeGrid_Df_scaled <- varImp_DeGrid_Df%>%group_by(VarName)%>%summarize(VarImp_Mean=mean(VarImp),
                                                                            VarImp_SD=sd(VarImp), 
                                                                            VarImp_SE=std.error(VarImp))
varImp_DeGrid_Df_scaled$VarName <- factor(varImp_DeGrid_Df_scaled$VarName, levels = levels(varImp_GMP_Df_scaled$VarName))

## Figure 5A: global analysis
ggplot(varImp_GMP_Df_scaled, aes(x=VarName, y=VarImp_Mean, fill=VarImp_Mean))+
  geom_bar(stat="identity", position="dodge", width=0.7)+coord_flip()+
  scale_fill_gradient(low="blue", high="red")+
  geom_errorbar(aes(ymin=VarImp_Mean-VarImp_SD, ymax=VarImp_Mean+VarImp_SD), width=.12,
                position=position_dodge(0.05))+
  xlab("Variables")+
  ylab("Importance")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_blank(),
        legend.position="none")

## Figure 5B left: deciduous-dominated clusters
ggplot(varImp_DeGrid_Df_scaled, aes(x=VarName, y=VarImp_Mean, fill=VarImp_Mean))+
  geom_bar(stat="identity", position="dodge", width=0.7)+coord_flip()+
  scale_fill_gradient(low="blue", high="red")+
  geom_errorbar(aes(ymin=VarImp_Mean-VarImp_SD, ymax=VarImp_Mean+VarImp_SD), width=.12,
                position=position_dodge(0.05))+
  ylim(c(0,0.06))+
  xlab("Variables")+
  ylab("Importance")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=16, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_blank(),
        legend.position="none")

## Figure 5B middle: bimodal clusters
ggplot(varImp_BiGrid_Df_scaled, aes(x=VarName, y=VarImp_Mean, fill=VarImp_Mean))+
  geom_bar(stat="identity", position="dodge", width=0.7)+coord_flip()+
  scale_fill_gradient(low="blue", high="red")+
  geom_errorbar(aes(ymin=VarImp_Mean-VarImp_SD, ymax=VarImp_Mean+VarImp_SD), width=.12,
                position=position_dodge(0.05))+
  ylim(c(0,0.08))+
  xlab("Variables")+
  ylab("Importance")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_blank(),
        legend.position="none")

## Figure 5B right: evergreen-dominated clusters
ggplot(varImp_EvGrid_Df_scaled, aes(x=VarName, y=VarImp_Mean, fill=VarImp_Mean))+
  geom_bar(stat="identity", position="dodge", width=0.7)+coord_flip()+
  scale_fill_gradient(low="blue", high="red")+
  geom_errorbar(aes(ymin=VarImp_Mean-VarImp_SD, ymax=VarImp_Mean+VarImp_SD), width=.12,
                position=position_dodge(0.05))+
  ylim(c(0,0.06))+
  xlab("Variables")+
  ylab("Importance")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_blank(),
        legend.position="none")
