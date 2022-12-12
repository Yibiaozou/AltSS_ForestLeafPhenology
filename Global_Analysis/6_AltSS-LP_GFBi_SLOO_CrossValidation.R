rm(list=ls())
library(tidyverse)
library(geosphere)
library(ranger)
library(tictoc)
library(caret)

# load the cluster-level dataset which was used to train the spatial random forest model
GMP_Df_Full =read.csv("**/Global_Analysis/Data/GEE_RF/Input_FSD/GMP_Df_FSD.csv")
AltSS_LP_RF_Ranger_FSD <- readRDS("**/Global_Analysis/Data/AltSS_LP_RF_Caret_FSD_AllCovariates.rds")

covariateList <- colnames(GMP_Df_Full)
covariateList <- covariateList[-which(covariateList%in%c("BI", "latitude", "longitude", "WWF_Biome"))]
formula_Rf <- as.formula(paste0("BI~",paste(covariateList, collapse="+")))

# range of buffer size in meter
Buffer_Size_List <- c(0, 10000, 100000, 1000000, 250000,  500000, 750000)

SLOOCV_Df <- NULL

# analysis across the given range of buffer size
for(i in 1:length(Buffer_Size_List)){
  predict_SLOOCV <- as.data.frame(cbind(pred=rep(NA, nrow(GMP_Df_Full)), obs=rep(NA, length(valSampleID))))
  predict_SLOOCV$latitude <- GMP_Df_Full$latitude
  predict_SLOOCV$longitude <- GMP_Df_Full$longitude
  predict_SLOOCV$Buffer_Size <- as.numeric(Buffer_Size_List[i])
  ## implement spatial leave-one-out cross validation for each row of data
  for(j in 1:nrow(GMP_Df_Full)){
    GMP_Df_Test <-  GMP_Df_Full[j,]
    DisAll <- distm(GMP_Df_Test[,c('longitude','latitude')], GMP_Df_Full[,c('longitude','latitude')], fun = distHaversine)
    GMP_Df_Train <-  GMP_Df_Full[DisAll>buffer_size,]
    
    AltSS_LP_RF_Ranger <- ranger(
      formula=formula_Rf,
      data=GMP_Df_Train,
      num.trees=250,
      mtry=20, # best tuned
      min.node.size = 2, # best tuned
      splitrule = "maxstat" # best tuned
    )
    
    predict_SLOOCV$pred[j] <- predict(AltSS_LP_RF_Ranger, data=GMP_Df_Test)$predictions
    predict_SLOOCV$obs[j] <- GMP_Df_Test$BI
  }
  
  SLOOCV_Df <- rbind(SLOOCV_Df, predict_SLOOCV)
}

saveRDS(SLOOCV_Df, file="**/Global_Analysis/Model_Output/SLOOCV_Df_FSD.rds")
