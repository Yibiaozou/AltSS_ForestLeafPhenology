rm(list=ls())

## load packages
library(tidyverse)
library(tictoc)
library(geosphere)
library(mousetrap)
library(diptest)
library(moments)
library(raster)
library(BimodalIndex)
library(LaplacesDemon)
library(MASS)
library(ggpubr)
library(rgdal)
library(data.table)
# library("rnaturalearth")
# library("rnaturalearthdata")
# library(rgeos)
library(matrixStats)
library(factoextra)
library(ggfortify)
library(caret)

# function to get density of one variable along another variable, for visualization
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## load tree-level GFBI dataset
GFBI_Df1 <- read.csv("**/Global_Analysis/Data/GFBI_Tree_Level.csv") 

#################################################################################################
####---------Extract leaf phenology and morphological type to each species in GFBI dataset---####
## load species and genus-level information of leaf phenology and leaf form  based on TRY dataset
load("**/Global_Analysis/Data/Leaf_Phenology_Df.RData")
load("**/Global_Analysis/Data/Leaf_Type_Df.RData")

GFBI_Df1 <- GFBI_Df1%>%mutate(Leaf_Phenology=NA, Leaf_Type=NA)

# get all unique tree species in GFBi data
SppLst1 <- unique(GFBI_Df1$SPCD)
Len1 <- length(SppLst1)

# data frame to store leaf phenology information (evergreen vs. deciduous) of each species
LP_explore_Df <- as.data.frame(cbind(SppName = SppLst1, GenName = rep(NA, Len1), 
                                     Leaf_Phenology = rep(NA, Len1)))
# data frame to store leaf form information (needleleaf vs. broadleaf) of each species
LT_explore_Df <- as.data.frame(cbind(SppName = SppLst1, GenName = rep(NA, Len1), 
                                     Leaf_Type = rep(NA, Len1)))

# get genus name for each species
LP_explore_Df$GenName <- sapply(1:Len1, function(x) strsplit(LP_explore_Df$SppName[x], split=" ")[[1]][1])
LT_explore_Df$GenName <- sapply(1:Len1, function(x) strsplit(LT_explore_Df$SppName[x], split=" ")[[1]][1])

LP_Name <- c("Deciduous", "Evergreen")
LT_Name <- c("Broadleaf", "Needleleaf")

tic()

# Species level search
for(i in 1:Len1){
  LP_ts_Vector <- Leaf_Phenology_Df$Leaf_Phenology[Leaf_Phenology_Df$AccSpeciesName==LP_explore_Df$SppName[i]]
  LP_ts_Table <- table(LP_ts_Vector) 
  if(length(LP_ts_Table)==2){  
    LP_explore_Df$Leaf_Phenology[i] <- LP_Name[which(LP_ts_Table==max(LP_ts_Table))[1]]
  }else if(length(LP_ts_Table)==1){
    LP_explore_Df$Leaf_Phenology[i] <- LP_ts_Vector[1]
  }
  
  LT_ts_Vector <- Leaf_Type_Df$Leaf_Type[Leaf_Type_Df$AccSpeciesName==LT_explore_Df$SppName[i]]
  LT_ts_Table <- table(LT_ts_Vector) 
  if(length(LT_ts_Table)==2){  
    LT_explore_Df$Leaf_Type[i] <- LT_Name[which(LT_ts_Table==max(LT_ts_Table))[1]]
  }else if(length(LT_ts_Table)==1){
    LT_explore_Df$Leaf_Type[i] <- LT_ts_Vector[1]
  }
}

toc()

tic()
# Genus level search
for(i in 1:Len1){
  if(is.na(LP_explore_Df$Leaf_Phenology[i])){
    LP_ts_Vector <- Leaf_Phenology_Df$Leaf_Phenology[Leaf_Phenology_Df$GenName==LP_explore_Df$GenName[i]]
    LP_ts_Table <- table(LP_ts_Vector) 
    if(length(LP_ts_Table)==2){  
      LP_explore_Df$Leaf_Phenology[i] <- LP_Name[which(LP_ts_Table==max(LP_ts_Table))[1]]
    }else if(length(LP_ts_Table)==1){
      LP_explore_Df$Leaf_Phenology[i] <- LP_ts_Vector[1]
    }
  }
  
  if(is.na(LT_explore_Df$Leaf_Type[i])){
    LT_ts_Vector <- Leaf_Type_Df$Leaf_Type[Leaf_Type_Df$GenName==LT_explore_Df$GenName[i]]
    LT_ts_Table <- table(LT_ts_Vector)
    if(length(LT_ts_Table)==2){  
      LT_explore_Df$Leaf_Type[i] <- LT_Name[which(LT_ts_Table==max(LT_ts_Table))[1]]
    }else if(length(LT_ts_Table)==1){
      LT_explore_Df$Leaf_Type[i] <- LT_ts_Vector[1]
    }
  }
}
toc()

## add all species/genus level information back to GFBi dataset
tic()
for(i in 1:Len1){
  if(!is.na(LP_explore_Df$Leaf_Phenology[i])){
    GFBI_Df1$Leaf_Phenology1[GFBI_Df1$SPCD==LP_explore_Df$SppName[i]] <- LP_explore_Df$Leaf_Phenology[i]
  }
  if(!is.na(LT_explore_Df$Leaf_Type[i])){
    GFBI_Df1$Leaf_Type1[GFBI_Df1$SPCD==LT_explore_Df$SppName[i]] <- LT_explore_Df$Leaf_Type[i]
  }
}
toc()


GFBI_Df1 <- GFBI_Df1%>%mutate(BA=pi*(DBH/2)^2, # basal area of each individual tree
                              LP_known=ifelse(is.na(Leaf_Phenology1),0,1), # binary value to show whether leaf phenology information is available for the tree
                              EV_Status=ifelse(Leaf_Phenology1=="Evergreen", 1, 0), # binary value to show whether the tree is evergreen or not
                              DE_Status=ifelse(Leaf_Phenology1=="Deciduous", 1, 0), # binary value to show whether the tree is deciduous or not
                              LT_known=ifelse(is.na(Leaf_Type1),0,1), # binary value to show whether leaf form information is available for the tree
                              BL_Status=ifelse(Leaf_Type1=="Broadleaf", 1, 0), # binary value to show whether the tree is broadleaf or not
                              NL_Status=ifelse(Leaf_Type1=="Needleleaf", 1, 0)) # binary value to show whether the tree is needleleaf or not

GFBI_Df1 <- GFBI_Df1%>%mutate(LP_BA=BA*LP_known, # basal area of the tree with known leaf phenology information
                              EV_BA=BA*EV_Status, # basal area of evergreen tree
                              DE_BA=BA*DE_Status, # basal area of deciduous tree
                              LT_BA=BA*LT_known, # basal area of the tree with known leaf form information
                              BL_BA=BA*BL_Status, # basal area of broadleaf tree
                              NL_BA=BA*NL_Status) # basal area of needleleaf tree
save(GFBI_Df1, file="**/Global_Analysis/Data/GFBI_Lab_Full.R")


##################################################################
####---------Aggregating and calculating plot-level metrics---####
tic()
## define unique plot id for each plot in GFBi dataset
GFBI_Df1$PLTN <- sapply(1:nrow(GFBI_Df1), function(x)paste0(GFBI_Df1$PLT[x],"_",GFBI_Df1$LON[x],"_",GFBI_Df1$LAT[x]))
toc()



## function to get relative basal area of the tree species with the largest cumulative basal area in a plot
spp_dominance <- function(x){
  TT <- table(na.omit(x))
  maxProp <- max(TT)/sum(TT)
  return(maxProp)
}

## aggregate tree-level GFBi data to plot-level metrics, such as total basal area, basal area of evergreen trees. 
tic()
GFBI_Df2_aggregated <- GFBI_Df1%>%group_by(PLTN)%>%summarise(longitude=mean(LON,na.rm=T), # plot longitude
                                                             latitude=mean(LAT,na.rm=T), # plot latitude
                                                             BA=sum(BA, na.rm=T), # plot total basal area
                                                             LP_BA=sum(LP_BA, na.rm=T), # plot basal area of all trees with known leaf phenology information
                                                             LT_BA=sum(LP_BA, na.rm = T), # plot basal area of all trees with known leaf form information
                                                             # Myco_BA=sum(Myco_BA, na.rm = T),
                                                             EV_BA=sum(EV_BA, na.rm = T), # plot basal area of evergreen trees
                                                             BL_BA=sum(BL_BA, na.rm = T), # plot basal area of broadleaf trees
                                                             # EM_BA=sum(EM_BA, na.rm = T),
                                                             spp.count=length(unique(SPCD)), # number of species in the plot
                                                             spp_maxprop=spp_dominance(SPCD), # the relative basal area of the tree species with the largest cumulative basal area in the plot
                                                             count=sum(count, na.rm=T), # number of trees in the plot
                                                             EV_Density=sum(EV_Density,na.rm=T), # number of evergreen trees in the plot
                                                             DE_Density=sum(DE_Density,na.rm=T)) # number of deciduous trees in the plot
toc()

# set a simpler numeric plot id
GFBI_Df2_aggregated$Plot_id <- 1:nrow(GFBI_Df2_aggregated)

# extract coordinates of the plot-level GFBi data
GFBI_Df2_coordinates <- GFBI_Df2_aggregated[,c("Plot_id","longitude", "latitude")]%>%drop_na()
write.csv(GFBI_Df2_coordinates, file="**/Global_Analysis/Data/GFBI_Df2_coordinates_new.csv", row.names = F)

### spatial covariates extracted from Crowther Lab Composite based on the provided coordinated of the dataset. The extracting code is provided in '.py'.
GFBI_Df2_covariates <- read.csv("**/Global_Analysis/Data/20221013_GFBI_Yibiao.csv")
## add spatial covariates to the plot-level GFBi data
GFBI_Df2_Aggregated_Full_2 <- left_join(GFBI_Df2_covariates, GFBI_Df2_covariates, by="Plot_id")
GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2 %>% mutate(relEV=EV_BA/LP_BA, # calculate plot-level relative evergreen abundance weighted by basal area
                                                                    SampleSize = 1) # set sample size of each plot as 1


#######################################################################################################################
####---------Spatial Clustering: Aggregating and calculating cluster-level metrics based on spatial partitioning---####
### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree, i.e., 10 arc-min. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot, which later used for 
GFBI_Df2_Aggregated_Full_2$Cluster_10min <- raster::extract(fishNet, GFBI_Df2_Aggregated_Full_2[,c("longitude", "latitude")])
### save the plot-level GFBi data
saveRDS(GFBI_Df2_Aggregated_Full_2, file="**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")

## get the names of spatial covariates
GMP_Df_Full =read.csv("**/Global_Analysis/Data/GEE_RF/Input_Full/GMP_Df_Full.csv")
CoVaList <- c(colnames(GMP_Df_Full)[4:81])

#  only keep plots for which information on leaf phenology strategies was available for > 90% of trees (weighted by basal area) 
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$LP_BA/GFBI_Df2_Aggregated_Full_2$BA>=0.9, ]

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2_filtered[GFBI_Df2_Aggregated_Full_2_filtered$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2_filtered$count>=10 & GFBI_Df2_Aggregated_Full_2_filtered$spp.count>=2,]



## function to compute the mode of a vector
mode <- function(x,na.rm=T){
  if(na.rm==T){
    tb=table(na.omit(x))
  }else{tb=table(x)}
  
  out=as.numeric(names(tb)[which.max(tb)])
  out=out[1]
  if(length(out)==0){
    out=NA
  }
  return(out)
}

## aggregate the plot-level dataset to cluster-level dataset
tic()
GFBI_Df3_Scaled_10min <- GFBI_Df2_Aggregated_Full_2_filtered%>%
  group_by(Cluster_10min)%>%summarise(latitude=mean(latitude), # cluster-mean latitude
                                      longitude=mean(longitude), # cluster-mean longitude
                                      dip_Stat=dip.test(relEV)$statistic, # dip-test statistics
                                      dip_P=dip.test(relEV)$p.value, # dip-test p-value
                                      skewness=skewness(relEV), # skewness of distribution of plot-level relEV
                                      kurtosis=kurtosis(relEV), # kurtosis of distribution of plot-level relEV
                                      SampleSize=sum(SampleSize), # total sample size in the cluster
                                      relEV_mu=mean(relEV, na.rm=T), # cluster-mean relEV
                                      WWF_Biome = mode(WWF_Biome), # the majority biome group of the plots within the cluster belong to
                                      across(CoVaList, mean, na.rm=T) # cluster-mean values of spatial environmental covariates
  )
toc()

GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min %>% mutate(dip_Stat_AD=sqrt(SampleSize)*dip_Stat) # adapted dip-test statistics, calculated based on Hartigan & Hatigan (1985)

# compute the bimodality index (BI). The BI ranges from -1 to 1, 
# and we empirically derived bimodality cutoffs (see methods in the paper), 
# with BIs < -0.22 representing deciduous-dominated clusters, 
# BIs > 0.22 representing evergreen-dominated clusters, 
# and BIs of -0.22 – 0.22 representing bimodal clusters.
GFBI_Df3_Scaled_10min <- GFBI_Df3_Scaled_10min %>% mutate(BI=-exp(-(dip_Stat_AD)^2*6)*sign(skewness)) %>% drop_na() 
saveRDS(GFBI_Df3_Scaled_10min, file="**/Global_Analysis/Data/GFBI_Df3_Scaled_10min_FullStandardized.rds")

GFBI_Df3_scaled_10min_RF <- GFBI_Df3_Scaled_10min%>%dplyr::select(-Cluster_10min, -dip_Stat, 
                                                                  -dip_P, -skewness, -SampleSize, -dip_Stat_AD, -density, -skewness, -kurtosis, -relEV_mu)
GFBI_Df3_scaled_10min_RF <- GFBI_Df3_scaled_10min_RF %>%drop_na()
write.csv(GFBI_Df3_scaled_10min_RF, file="**/Global_Analysis/Data/GEE_RF/Input_FSD/GMP_Df_FSD.csv", row.names = F)


##### Hyperparamter tuning
### training and hyperparameter tuning
library(caret)
library(ranger)

covariateList <- colnames(GFBI_Df3_scaled_10min_RF)
covariateList <- covariateList[-which(covariateList%in%c("BI", "latitude", "longitude", "WWF_Biome"))]

## model BI as a function of all selected environmental covariates
formula_Rf <- as.formula(paste0("BI~",paste(covariateList, collapse="+")))

# set seed for reproducibility 
set.seed(2022)

tic()
AltSS_LP_RF_Caret_FSD_AllCovariates <- train(
  form = formula_Rf, 
  data = GFBI_Df3_scaled_10min_RF, 
  method = "ranger",
  # ranges for hyperparameter tuning
  tuneGrid = expand.grid( .mtry = c(1, 2, 4, 5, 8, 10, 15, 20), # number of variables sampled at each split
                          .min.node.size = c(1, 2, 5, 10, 20), # the minimum sample size at the end of the nodes
                          .splitrule = c("maxstat")),
  #c("variance", "maxstat")),
  trControl = trainControl(method = "cv", number = 10), # 10-fold cross-validation
  importance = 'impurity'
  
)
toc()
saveRDS(AltSS_LP_RF_Caret_FSD_AllCovariates, file="**/Global_Analysis/Data/AltSS_LP_RF_Caret_FSD_AllCovariates.rds")


###################################################################################################################################
####---------Environmental Clustering: Aggregating and calculating cluster-level metrics based on environmental partitioning---####
GFBI_Df2_Aggregated_Full_2 <- readRDS("**/Global_Analysis/Data/GFBI_Df2_aggregated_New.rds")

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2_Env_Kmeans <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2$count>=10 & GFBI_Df2_Aggregated_Full_2$spp.count>=2,]

GFBI_Df2_Aggregated_Full_2_Env_Kmeans <- GFBI_Df2_Aggregated_Full_2_Env_Kmeans[,c(covariateList, 
                                                                                    "relEV", "SampleSize",
                                                                                    "latitude", "longitude", "WWF_Biome")]%>%drop_na()

tic()
## PCA across all spatial environmental covariates
env.pca <- prcomp(GFBI_Df2_Aggregated_Full_2_Env_Kmeans[,covariateList], scale = T)
## get the variation captured by the leading 3 PCs
var.3 <- sum(summary(env.pca)$importance[2,1:3]) * 100
cat(paste0(round(var.3,3),'% of variation capture in by the first 3 principle components./n'))
toc()

## add the leading 3 PCs to the dataset
GFBI_Df2_Aggregated_Full_2_Env_Kmeans <- GFBI_Df2_Aggregated_Full_2_Env_Kmeans %>% mutate(PC1 = env.pca$x[,1], 
                                                                                          PC2 = env.pca$x[,2], 
                                                                                          PC3 = env.pca$x[,3] 
                                                                                            # PC4 = env.pca$x[,4], PC5 = env.pca$x[,5],
                                                                                            # PC6 = env.pca$x[,6], PC7 = env.pca$x[,7],
                                                                                            # PC8 = env.pca$x[,8], PC9 = env.pca$x[,9],
                                                                                            # PC10 = env.pca$x[,10]
                                                                                            )
# extract the columns of the leading 3 PCs
xy_GFBI <- GFBI_Df2_Aggregated_Full_2_Env_Kmeans%>%dplyr::select(PC1, PC2, PC3)

# implement K means clustering on the leading 3 PCs
tic()
k1 <- kmeans(xy_GFBI, centers = 15000, nstart = 5, iter.max=30)
toc()

## add the environmental clusters to the dataset
GFBI_Df2_Aggregated_Full_2_Env_Kmeans$Cluster_PCA <- k1$cluster

## save the dataset
saveRDS(GFBI_Df2_Aggregated_Full_2_Clim_Kmeans, file="**/Global_Analysis/Data/GFBI_Df2_Aggregated_Full_2_PCA_Kmeans_FSD.rds")

tic()
## aggregate the plot-level dataset to cluster-level dataset
GFBI_Df3_Scaled_PCA_Full_Kmeans <- GFBI_Df2_Aggregated_Full_2_Clim_Kmeans%>%
  group_by(Cluster_PCA)%>%summarise(
    latitude=mean(latitude), # cluster-mean latitude
    longitude=mean(longitude), # cluster-mean longitude
    dip_Stat=dip.test(relEV)$statistic, # dip-test statistics
    dip_P=dip.test(relEV)$p.value, # dip-test p-value
    skewness=skewness(relEV), # skewness of distribution of plot-level relEV
    WWF_Biome = mode(WWF_Biome), # the majority biome group of the plots within the cluster belong to
    SampleSize=sum(SampleSize), # total sample size in the cluster
    PC1 = mean(PC1), # cluster-mean PC1
    PC2 = mean(PC2), # cluster-mean PC2
    PC3 = mean(PC3), # cluster-mean PC3
    across(covariateList, mean, na.rm=T) # cluster-mean values of environmental covariates
  )
toc()


write.csv(GFBI_Df3_Scaled_PCA_Full_Kmeans, 
          file="**/Global_Analysis/Data/GEE_RF/Input_PCA/GMP_Df_PCA_Kmeans_FSD.csv",
          row.names=F)

##### Hyperparamter tuning
### training and hyperparameter tuning
library(caret)
library(ranger)

GFBI_Df3_Scaled_PCA_Full_Kmeans <- readRDS("**/Global_Analysis/Data/GFBI_Df3_Scaled_PCA_Full_Kmeans_FSD.rds")

## model BI as a function of all selected environmental covariates
formula_Rf <- as.formula(paste0("BI~",paste(covariateList, collapse="+")))

GFBI_Df3_Scaled_PCA_Full_Kmeans <- GFBI_Df3_Scaled_PCA_Full_Kmeans%>%drop_na()

# set seed for reproducibility 
set.seed(2022)
tic()
AltSS_LP_RF_PCA_Full_Kmeans <- train(
  form = formula_Rf, 
  data = GFBI_Df3_Scaled_PCA_Full_Kmeans,
  method = "ranger",
  # ranges for hyperparameter tuning
  tuneGrid = expand.grid( .mtry = c(1, 2, 4, 5, 8, 10), # number of variables sampled at each split
                          .min.node.size = c(1, 2, 5, 10, 20), # the minimum sample size at the end of the nodes
                          .splitrule = c("maxstat")),
  trControl = trainControl(method = "cv", number = 10), # 10-fold cross validation
  importance = 'impurity'
  
)

toc()
saveRDS(AltSS_LP_RF_PCA_Full_Kmeans, file="**/Global_Analysis/Data/AltSS_LP_RF_PCA_Full_Kmeans_FSD.rds")

