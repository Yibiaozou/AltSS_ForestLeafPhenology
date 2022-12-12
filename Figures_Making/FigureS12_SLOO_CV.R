rm(list=ls())
library(tidyverse)
library(Metrics)

#### function to calculate R2
rsq <- function(x, y) summary(lm(y~x))$r.squared

####----Visualize results of SLOOCV----####
# load the results of SLOO-CV
SLOOCV_Df_All <- readRDS("**/Global_Analysis/Model_Output/SLOOCV_Df_FSD.rds")

Buffer_Size_Name_List <- c(0, 10000, 100000, 1000000, 250000,  500000, 750000)
SLOOCV_Df_Results <- as.data.frame(cbind(buffer_Size=Buffer_Size_Name_List, R_2=NA, RMSE=NA))
for (i in 1:length(Buffer_Size_Name_List)){
  SLOOCV_Df_Sub <- SLOOCV_Df_All[SLOOCV_Df_All$Buffer_Size==Buffer_Size_Name_List[i],]
  SLOOCV_Df_Results$R_2[i] <- rsq(SLOOCV_Df_Sub$obs, SLOOCV_Df_Sub$pred) # compute R^2 between observation and model prediction
  SLOOCV_Df_Results$RMSE[i] <- rmse(SLOOCV_Df_Sub$obs, SLOOCV_Df_Sub$pred) # compute RMSE of model prediction
} 

SLOOCV_Df_Results$buffer_Size_factor <- as.factor(SLOOCV_Df_Results$buffer_Size/1000)
# fix the order of buffer size as decreasing order
SLOOCV_Df_Results$buffer_Size_factor <- fct_reorder(SLOOCV_Df_Results$buffer_Size_factor , SLOOCV_Df_Results$buffer_Size , .desc = F) 

## Figure S12A
ggplot(SLOOCV_Df_Results, aes(x=buffer_Size/1000, y=R_2))+
  geom_point(size=4, alpha=0.7)+
  geom_line(size=2, alpha=0.5)+
  ylim(c(0,0.75))+
  xlab("Buffer radius (km)")+
  ylab(expression(paste(italic(R)^2)))+
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

## Figure S12B
ggplot(SLOOCV_Df_Results, aes(x=buffer_Size/1000, y=RMSE))+
  geom_point(size=4, alpha=0.7)+
  geom_line(size=2, alpha=0.5)+
  ylim(c(0,0.47))+
  xlab("Buffer radius (km)")+
  ylab(expression(paste(italic(RMSE))))+
  theme_classic()+
  # scale_x_log10()+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))

####----Random 10-fold CV----####
GMP_Df_Full =read.csv("**/Global_Analysis/Data/GEE_RF/Input_FSD/GMP_Df_FSD.csv")


covariateList <- colnames(GMP_Df_Full)
covariateList <- covariateList[-which(covariateList%in%c("BI", "latitude", "longitude", "WWF_Biome"))]
formula_Rf <- as.formula(paste0("BI~",paste(covariateList, collapse="+")))

## randomly allocate cross-validation folds
flds_index <- sample(1:10, size=nrow(GMP_Df_Full), prob = rep(1/10, 10), replace = T)
predict_RandCV <- as.data.frame(cbind(ID=flds_index, 
                                      pred=rep(NA, length(flds_index)), 
                                      obs=rep(NA, length(flds_index)), 
                                      latitude=rep(NA, length(flds_index)), 
                                      longitude=rep(NA, length(flds_index))
))

# random cross-validation
for (i in 1:10){
  testing_set <- GMP_Df_Full[flds_index==i,]
  training_set <- GMP_Df_Full[flds_index!=i,]
  
  tic()
  
  AltSS_LP_RF_Ranger <- ranger(
    formula=formula_Rf,
    data=training_set,
    num.trees=250,
    mtry=20, # best tuned
    min.node.size = 2, # best tuned
    splitrule = "maxstat" # best tuned
  )
  
  predict_RandCV$pred[flds_index==i] <- predict(AltSS_LP_RF_Ranger, data=testing_set)$predictions
  predict_RandCV$obs[flds_index==i] <- testing_set$BI
  predict_RandCV$latitude[flds_index==i] <- testing_set$latitude
  predict_RandCV$longitude[flds_index==i] <- testing_set$longitude
  
  toc()
}

####---Residual semi-variogram---####
library(gstat)
library(sp)
library(mgcv)

# compute residuals
SLOOCV_Df_All$residuals <- SLOOCV_Df_All$obs - SLOOCV_Df_All$pred
predict_RandCV$residuals <- predict_RandCV$obs - predict_RandCV$pred

# compute semi-variance for residuals of all models
Variogram_Df <- NULL
coordinates(predict_RandCV) <- ~ longitude + latitude
proj4string(predict_RandCV)  <- CRS("+init=epsg:4326")
vario       <-     variogram(residuals~ 1, data = predict_RandCV, cutoff=1000, width=10) 
Variogram_Df <- rbind(Variogram_Df, as.data.frame(cbind(group="Random CV", value=vario$gamma, dist=vario$dist)))

for (i in 1:length(Buffer_Size_Name_List)){
  SLOOCV_Vario_Df_Sub <- SLOOCV_Df_All[SLOOCV_Df_All$Buffer_Size==Buffer_Size_Name_List[i],] %>% dplyr::select(residuals, latitude, longitude)
  SLOOCV_Vario_Df_Sub <- SLOOCV_Vario_Df_Sub%>%drop_na()
  coordinates(SLOOCV_Vario_Df_Sub) <- ~ longitude + latitude
  proj4string(SLOOCV_Vario_Df_Sub)  <- CRS("+init=epsg:4326")
  vario       <-     variogram(residuals~ 1, data = SLOOCV_Vario_Df_Sub, cutoff=1000, width=10)           #get semivariance values # cutoff=1000, width=5
  Variogram_Df <- rbind(Variogram_Df, as.data.frame(cbind(group=paste0("SLOO-CV (",Buffer_Size_Name_List[i]/1000,"km)"), value=vario$gamma, dist=vario$dist)))
}

Variogram_Df$radius <- c(rep(Buffer_Size_Name_List, each=100), rep(0, 100))
Variogram_Df$group <- factor(Variogram_Df$group, levels=c("Random CV", "SLOO-CV (0km)", "SLOO-CV (10km)",
                                                          "SLOO-CV (100km)", "SLOO-CV (250km)", "SLOO-CV (500km)",
                                                          "SLOO-CV (750km)", "SLOO-CV (1000km)"))

## Figure S12C
library("viridis")  
ggplot(Variogram_Df, aes(x=dist, y=value, group=group, color=group))+
  geom_line(size=2, alpha=0.8,)+
  scale_color_manual(values=viridis(255)[seq(1,250, length=8)])+
  ylim(c(0,0.18))+
  xlab("Distance (km)")+
  ylab("Semivariance")+
  theme_classic()+
  geom_vline(xintercept=500, color="red", size=1.5, linetype = "dashed")+
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"))


