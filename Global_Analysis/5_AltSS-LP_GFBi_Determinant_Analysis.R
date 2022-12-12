rm(list=ls())
library(tidyverse)
library(caret)
library(tictoc)
library(data.table)
library(plotrix)
library(diptest)


####-----Boostrapping determinant analysis on the seven selected variables---####
# load the plot-level GFBi dataset
GFBI_Df2_Aggregated_Full_2 <- readRDS("**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")
# load the cluster-level GFBi dataset
GFBI_Df3_Scaled_10min_Full <- readRDS("**/Global_Analysis/Data/GFBI_Df3_Scaled_10min_FullStandardized.rds")

GMP_Df_Full =read.csv("**/Global_Analysis/Data/GEE_RF/Input_FSD/GMP_Df_FSD.csv")

# retain cluster with sample size >= 50, to make sure BI can accurately depict the distribution
GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min_Full[GFBI_Df3_Scaled_10min_Full$SampleSize>=50, ]%>%drop_na()

# index of bimodal clusters
Bimodal_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI<0.22 & GFBI_Df3_Scaled_10min$BI> -0.22]
# index of evergreen-dominated clusters
EV_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI>0.22]
# index of deciduous-dominated clusters
DE_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI< -0.22]

#  only keep plots for which information on leaf phenology strategies was available for > 90% of trees (weighted by basal area) 
GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$LP_BA/GFBI_Df2_Aggregated_Full_2$BA>=0.9, ]

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2$count>=10 & GFBI_Df2_Aggregated_Full_2$spp.count>=2,]

# forest plot data within bimodal clusters
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%Bimodal_Cluster_index,]%>%drop_na()
# forest plot data within evergreen-dominated clusters
GFBI_Df2_EV_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%EV_Cluster_index,]%>%drop_na()
# forest plot data within deciduous-dominated clusters
GFBI_Df2_DE_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%DE_Cluster_index,]%>%drop_na()

# target variable for global anlaysis
target_globe <- "BI"
# target variable for plot-level analysis in three types of clusters
target_clusters <- "relEV"

## the original variable names in the dataset
covariables_RfLst_uniVari <- c("CHELSA_Annual_Mean_Temperature", 
                               "CHELSA_Annual_Precipitation", 
                               "CHELSA_Mean_Temperature_of_Coldest_Quarter",
                               "CHELSA_Precipitation_of_Driest_Quarter",
                               "SG_Soil_pH_H2O_0_100cm",
                               "Nitrogen",
                               "cnRatio"
)

## short names for random forest training
Predictors_Lst <- c("MAT", 
                    "MAP", 
                    "Tco",
                    "Pdr",
                    "Soil_pH",
                    "Soil_N",
                    "Soil_CNratio"
)

# change original variable names to corresponding short names
setnames(GMP_Df_Full, covariables_RfLst_uniVari, Predictors_Lst)
setnames(GFBI_Df2_Bimodal_0, covariables_RfLst_uniVari, Predictors_Lst)
setnames(GFBI_Df2_EV_0, covariables_RfLst_uniVari, Predictors_Lst)
setnames(GFBI_Df2_DE_0, covariables_RfLst_uniVari, Predictors_Lst)

# only select predictors, target variable and spatial clusters
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Bimodal_0 %>%dplyr::select(Predictors_Lst, target_grids, Cluster_10min)
GFBI_Df2_EV_0 <- GFBI_Df2_EV_0 %>%dplyr::select(Predictors_Lst, target_grids, Cluster_10min)
GFBI_Df2_DE_0 <- GFBI_Df2_DE_0 %>%dplyr::select(Predictors_Lst, target_grids, Cluster_10min)
GMP_Df_Full <- GMP_Df_Full%>%dplyr::select(Predictors_Lst, target_globe)

### formula for random forest
# model BI as a function of the 7 selected covariates
formula_Rf_GMP_UniVar <- as.formula(paste0("BI~",paste(Predictors_Lst, collapse="+")))
# model plot-level relative evergreen abundance as a function of the 7 selected covariates
formula_Rf_Grid_UniVar <- as.formula(paste0("relEV~",paste(Predictors_Lst, collapse="+")))

# number of bootstrapping repetition
Nrep <- 100

varImp_GMP_Df <- NULL
varImp_BiGrid_Df <- NULL
varImp_EvGrid_Df <- NULL
varImp_DeGrid_Df <- NULL

for (i in 1:Nrep){
  set.seed(i)
  print(i)
  tic()
  
  # get bootstrapping cluster index
  Bimodal_Cluster_index_sub <- Bimodal_Cluster_index[sample(1:length(Bimodal_Cluster_index),size=round(length(Bimodal_Cluster_index)), replace=T)]
  EV_Cluster_index_sub <- EV_Cluster_index[sample(1:length(EV_Cluster_index),size=round(length(EV_Cluster_index)), replace=T)]
  DE_Cluster_index_sub <- DE_Cluster_index[sample(1:length(DE_Cluster_index),size=round(length(DE_Cluster_index)), replace=T)]
  GMP_index_sub <- sample(1:nrow(GMP_Df_Full), size=round(nrow(GMP_Df_Full)), replace=T)
  
  # get the bootstrapping dataset based on the bootstrapping index
  GFBI_Df2_Bimodal_sub <- GFBI_Df2_Bimodal_0[GFBI_Df2_Bimodal_0$Cluster_10min%in%Bimodal_Cluster_index_sub,]
  GFBI_Df2_EV_sub <- GFBI_Df2_EV_0[GFBI_Df2_EV_0$Cluster_10min%in%EV_Cluster_index_sub,]
  GFBI_Df2_DE_sub <- GFBI_Df2_DE_0[GFBI_Df2_DE_0$Cluster_10min%in%DE_Cluster_index_sub,]
  GMP_Df_sub <- GMP_Df_Full[GMP_index_sub,]
  
  ## Global BI analysis
  # train random forest to model BI as a function of the 7 selected covariates
  AltSS_LP_RF_Caret_GMP_UniVar <- train(
    form = formula_Rf_GMP_UniVar,
    data = GMP_Df_sub,
    method = "ranger",
    tuneGrid = expand.grid( .mtry = c(4),
                            .min.node.size = c(2),
                            .splitrule = c("maxstat")),
    trControl = trainControl(method="none"),
    importance = "permutation" # use permutation importance as metric of variable importance
  )
  toc()
  
  # extract information of variable importance from the fitted random forest model
  varImp_GMP_Df <- rbind(varImp_GMP_Df,
                         as.data.frame(cbind(Rep=rep(i, length(Predictors_Lst)), # id of the ith repetition
                                             VarName=rownames(varImp(AltSS_LP_RF_Caret_GMP_UniVar, scale = FALSE)$importance), # variable names
                                             VarImp=varImp(AltSS_LP_RF_Caret_GMP_UniVar, scale = FALSE)$importance$Overall))) # variable permutation importance
  
  ## Bimodal cluster analysis
  tic()
  
  # train random forest to model plot-level relative evergreen abundance as a function of the 7 selected covariates
  AltSS_LP_RF_Caret_BiGrid_UniVar <- train(
    form = formula_Rf_Grid_UniVar, 
    data = GFBI_Df2_Bimodal_sub, 
    method = "ranger",
    tuneGrid = expand.grid( .mtry = c(4),
                            .min.node.size = c(2),
                            .splitrule = c("maxstat")),
    #c("variance", "maxstat")),
    trControl = trainControl(method="none"),
    importance = "permutation"
  )
  toc()
  
  varImp_BiGrid_Df <- rbind(varImp_BiGrid_Df, 
                            as.data.frame(cbind(Rep=rep(i, length(Predictors_Lst)),
                                                VarName=rownames(varImp(AltSS_LP_RF_Caret_BiGrid_UniVar, scale = FALSE)$importance),
                                                VarImp=varImp(AltSS_LP_RF_Caret_BiGrid_UniVar, scale = FALSE)$importance$Overall))) 
  
  # ## Evergren-dominated cluster analysis
  tic()
  # train random forest to model plot-level relative evergreen abundance as a function of the 7 selected covariates
  AltSS_LP_RF_Caret_EvGrid_UniVar <- train(
    form = formula_Rf_Grid_UniVar,
    data = GFBI_Df2_EV_sub,
    method = "ranger",
    tuneGrid = expand.grid( .mtry = c(4),
                            .min.node.size = c(2),
                            .splitrule = c("maxstat")),
    #c("variance", "maxstat")),
    trControl = trainControl(method="none"),
    importance = "permutation"
  )
  toc()
  #
  varImp_EvGrid_Df <- rbind(varImp_EvGrid_Df,
                            as.data.frame(cbind(Rep=rep(i, length(Predictors_Lst)),
                                                VarName=rownames(varImp(AltSS_LP_RF_Caret_EvGrid_UniVar, scale = FALSE)$importance),
                                                VarImp=varImp(AltSS_LP_RF_Caret_EvGrid_UniVar, scale = FALSE)$importance$Overall)))
  # # 
  # # ## Deciduous-dominated cluster analysis
  tic()
  # train random forest to model plot-level relative evergreen abundance as a function of the 7 selected covariates
  AltSS_LP_RF_Caret_DeGrid_UniVar <- train(
    form = formula_Rf_Grid_UniVar,
    data = GFBI_Df2_DE_sub,
    method = "ranger",
    tuneGrid = expand.grid( .mtry = c(4),
                            .min.node.size = c(2),
                            .splitrule = c("maxstat")),
    #c("variance", "maxstat")),
    trControl = trainControl(method="none"),
    importance = "permutation"
  )


  varImp_DeGrid_Df <- rbind(varImp_DeGrid_Df,
                            as.data.frame(cbind(Rep=rep(i, length(Predictors_Lst)),
                                                VarName=rownames(varImp(AltSS_LP_RF_Caret_DeGrid_UniVar, scale = FALSE)$importance),
                                                VarImp=varImp(AltSS_LP_RF_Caret_DeGrid_UniVar, scale = FALSE)$importance$Overall)))

  saveRDS(varImp_GMP_Df, file="**/Global_Analysis/Model_Output/varImp_GMP_GFBI_Df_FSD.rds")
  saveRDS(varImp_BiGrid_Df, file="**/Global_Analysis/Model_Output/varImp_BiGrid_GFBI_Df_FSD.rds")
  saveRDS(varImp_EvGrid_Df, file="**/Global_Analysis/Model_Output/varImp_EvGrid_GFBI_Df_FSD.rds")
  saveRDS(varImp_DeGrid_Df, file="**/Global_Analysis/Model_Output/varImp_DeGrid_GFBI_Df_FSD.rds")

  toc()
}

####-----Determinant analysis on the three leading principle components of climate, soil and topography---####
# load the plot-level GFBi dataset
GFBI_Df2_Aggregated_Full_2 <- readRDS("**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")
# load the cluster-level GFBi dataset
GFBI_Df3_Scaled_10min_Full <- readRDS("**/Global_Analysis/Data/GFBI_Df3_Scaled_10min_FullStandardized.rds")

GMP_Df_Full =read.csv("**/Global_Analysis/Data/GEE_RF/Input_Full/GMP_Df_Full.csv")

GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min_Full[GFBI_Df3_Scaled_10min_Full$SampleSize>=10, ]%>%drop_na()

Bimodal_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI<0.22 & GFBI_Df3_Scaled_10min$BI> -0.22]
# index of evergreen-dominated clusters
EV_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI>0.22]
# index of deciduous-dominated clusters
DE_Cluster_index <- GFBI_Df3_Scaled_10min$Cluster_10min[GFBI_Df3_Scaled_10min$BI< -0.22]

# forest plot data within bimodal clusters
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%Bimodal_Cluster_index,]%>%drop_na()
# forest plot data within evergreen-dominated clusters
GFBI_Df2_EV_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%EV_Cluster_index,]%>%drop_na()
# forest plot data within deciduous-dominated clusters
GFBI_Df2_DE_0 <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$Cluster_10min%in%DE_Cluster_index,]%>%drop_na()

# get names of climatic variables
Covariates_Climate <- colnames(GMP_Df_Full)[c(4:27,41:42, 50, 52, 78, 81)]
# get names of soil variables
Covariates_Soil <- colnames(GMP_Df_Full)[c(40, 53:57, 64:65, 79:80)]
# get names of topographic variables
Covariates_Topography <- colnames(GMP_Df_Full)[28:39]

###---run PCA to get the leading 3 PCs for global determinant analysis
## climate
tic()
env.pca <- prcomp(GMP_Df_Full[,Covariates_Climate], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three climatic PCs 
GMP_Df_Full <- GMP_Df_Full %>% mutate(PC1_Clim = env.pca$x[,1], PC2_Clim = env.pca$x[,2],
                                      PC3_Clim = env.pca$x[,3])

## Soil 
tic()
env.pca <- prcomp(GMP_Df_Full[,Covariates_Soil], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three soil PCs 
GMP_Df_Full <- GMP_Df_Full %>% mutate(PC1_Soil = env.pca$x[,1], PC2_Soil = env.pca$x[,2],
                                      PC3_Soil = env.pca$x[,3])

## Topography 
tic()
env.pca <- prcomp(GMP_Df_Full[,Covariates_Topography], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three topographic PCs 
GMP_Df_Full <- GMP_Df_Full %>% mutate(PC1_Topo = env.pca$x[,1], PC2_Topo = env.pca$x[,2],
                                      PC3_Topo = env.pca$x[,3])

###---run PCA to get the leading 3 PCs for determinant analysis in bimodal clusters
## climate
tic()
env.pca <- prcomp(GFBI_Df2_Bimodal_0[,Covariates_Climate], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three climatic PCs 
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Bimodal_0 %>% mutate(PC1_Clim = env.pca$x[,1], PC2_Clim = env.pca$x[,2],
                                      PC3_Clim = env.pca$x[,3])

## Soil 
tic()
env.pca <- prcomp(GFBI_Df2_Bimodal_0[,Covariates_Soil], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three soil PCs 
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Bimodal_0 %>% mutate(PC1_Soil = env.pca$x[,1], PC2_Soil = env.pca$x[,2],
                                      PC3_Soil = env.pca$x[,3])

## Topography 
tic()
env.pca <- prcomp(GFBI_Df2_Bimodal_0[,Covariates_Topography], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three topographic PCs 
GFBI_Df2_Bimodal_0 <- GFBI_Df2_Bimodal_0 %>% mutate(PC1_Topo = env.pca$x[,1], PC2_Topo = env.pca$x[,2],
                                      PC3_Topo = env.pca$x[,3])

###---run PCA to get the leading 3 PCs for determinant analysis in evergreen-dominated clusters
## climate
tic()
env.pca <- prcomp(GFBI_Df2_EV_0[,Covariates_Climate], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three climatic PCs 
GFBI_Df2_EV_0 <- GFBI_Df2_EV_0 %>% mutate(PC1_Clim = env.pca$x[,1], PC2_Clim = env.pca$x[,2],
                                                    PC3_Clim = env.pca$x[,3])

## Soil 
tic()
env.pca <- prcomp(GFBI_Df2_EV_0[,Covariates_Soil], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three soil PCs 
GFBI_Df2_EV_0 <- GFBI_Df2_EV_0 %>% mutate(PC1_Soil = env.pca$x[,1], PC2_Soil = env.pca$x[,2],
                                                    PC3_Soil = env.pca$x[,3])

## Topography 
tic()
env.pca <- prcomp(GFBI_Df2_EV_0[,Covariates_Topography], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three topographic PCs 
GFBI_Df2_EV_0 <- GFBI_Df2_EV_0 %>% mutate(PC1_Topo = env.pca$x[,1], PC2_Topo = env.pca$x[,2],
                                                    PC3_Topo = env.pca$x[,3])

###---run PCA to get the leading 3 PCs for determinant analysis in deciduous-dominated clusters
## climate
tic()
env.pca <- prcomp(GFBI_Df2_DE_0[,Covariates_Climate], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three climatic PCs 
GFBI_Df2_DE_0 <- GFBI_Df2_DE_0 %>% mutate(PC1_Clim = env.pca$x[,1], PC2_Clim = env.pca$x[,2],
                                          PC3_Clim = env.pca$x[,3])

## Soil 
tic()
env.pca <- prcomp(GFBI_Df2_DE_0[,Covariates_Soil], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three soil PCs 
GFBI_Df2_DE_0 <- GFBI_Df2_DE_0 %>% mutate(PC1_Soil = env.pca$x[,1], PC2_Soil = env.pca$x[,2],
                                          PC3_Soil = env.pca$x[,3])

## Topography 
tic()
env.pca <- prcomp(GFBI_Df2_DE_0[,Covariates_Topography], scale = T)
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,2),'% of variation capture in by the first 3 principle components./n'))
toc()
# add the leading three topographic PCs 
GFBI_Df2_DE_0 <- GFBI_Df2_DE_0 %>% mutate(PC1_Topo = env.pca$x[,1], PC2_Topo = env.pca$x[,2],
                                          PC3_Topo = env.pca$x[,3])

#### Set up random forest models for variable importance analysis
covariables_RfLst_3PC <- c("PC1_Clim","PC2_Clim","PC3_Clim",
                           "PC1_Soil","PC2_Soil","PC3_Soil",
                           "PC1_Topo","PC2_Topo","PC3_Topo") 

## model BI as a function of the 9 PCs
formula_Rf_GMP <- as.formula(paste0("BI~",paste(covariables_RfLst_3PC, collapse="+")))

## model plot-level relative evergreen abundance as a function of the 9 PCs
formula_Rf_Grid <- as.formula(paste0("relEV~",paste(covariables_RfLst_3PC, collapse="+")))

# train random forest to model BI as a function of the 9PCs
tic()
AltSS_LP_RF_Caret_GMP <- train(
  form = formula_Rf_GMP, 
  data = GMP_Df_Full, 
  method = "ranger",
  tuneGrid = expand.grid( .mtry = c(4),
                          .min.node.size = c(2),
                          .splitrule = c("maxstat")),
  importance = "permutation"
)
toc()

# bimodal clusters                                        
# train random forest to model plot-level relative evergreen abundance as a function of the 9 PCs
tic()
AltSS_LP_RF_Caret_BiGrid <- train(
  form = formula_Rf_BiGrid, 
  data = GFBI_Df2_Bimodal_0, 
  method = "ranger",
  tuneGrid = expand.grid( .mtry = c(4),
                          .min.node.size = c(2),
                          .splitrule = c("maxstat")),
  trControl = trainControl(method = "cv", number = 5),
  importance = "permutation"
)
toc()

# evergreen-dominated clusters 
# train random forest to model plot-level relative evergreen abundance as a function of the 9 PCs
formula_Rf_EvGrid <- as.formula(paste0("relEV~",paste(covariables_RfLst_3PC, collapse="+")))
tic()
AltSS_LP_RF_Caret_EvGrid <- train(
  form = formula_Rf_EvGrid, 
  data = GFBI_Df2_EV_0, 
  method = "ranger",
  tuneGrid = expand.grid( .mtry = c(4),
                          .min.node.size = c(2),
                          .splitrule = c("maxstat")),
  trControl = trainControl(method = "cv", number = 5),
  importance = "permutation"
)
toc()

# deciduous-dominated clusters 
# train random forest to model plot-level relative evergreen abundance as a function of the 9 PCs
formula_Rf_DeGrid <- as.formula(paste0("relEV~",paste(covariables_RfLst_3PC, collapse="+")))
tic()
AltSS_LP_RF_Caret_DeGrid <- train(
  form = formula_Rf_DeGrid, 
  data = GFBI_Df2_DE_0, 
  method = "ranger",
  tuneGrid = expand.grid( .mtry = c(4),
                          .min.node.size = c(2),
                          .splitrule = c("maxstat")),
  trControl = trainControl(method = "cv", number = 5),
  importance = "permutation"
)
toc()

save(AltSS_LP_RF_Caret_GMP, file="**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_GMP_3GroupVar_3PC.RData")
save(AltSS_LP_RF_Caret_BiGrid, file="**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_BiGrid_3GroupVar_3PC.RData")
save(AltSS_LP_RF_Caret_EvGrid, file="**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_EvGrid_3GroupVar_3PC.RData")
save(AltSS_LP_RF_Caret_DeGrid, file="**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_DeGrid_3GroupVar_3PC.RData")


load("**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_GMP_3GroupVar_3PC.RData")
load("**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_BiGrid_3GroupVar_3PC.RData")
load("**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_EvGrid_3GroupVar_3PC.RData")
load("**/Global_Analysis/Model_Output/AltSS_LP_RF_Caret_DeGrid_3GroupVar_3PC.RData")

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