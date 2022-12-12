rm(list=ls())
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
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library(matrixStats)
library(factoextra)
library(ggfortify)
library(caret)
library(rowr)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


####-----------Preprocess FunDivEUROPE data------####

#### dataset for Germany, Switzerland, Spain and Sweden
FunDivEUROPE_plots_DE_SW_ES <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_DE_SW_ES.csv")
FunDivEUROPE_species <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_species.csv")
FunDivEUROPE_trees_DE_SW_ES <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_trees_DE_SW_ES.csv")

# rename the column of species id
colnames(FunDivEUROPE_species)[1] <- "speciesid"
# remove space in the string of plot id
FunDivEUROPE_trees_DE_SW_ES$plotcode <- gsub("[[:space:]]", "", FunDivEUROPE_trees_DE_SW_ES$plotcode)

# subset the plot-level data to the same plots of tree-level data
FunDivEUROPE_plots_DE_SW_ES <- FunDivEUROPE_plots_DE_SW_ES[FunDivEUROPE_plots_DE_SW_ES$plotcode%in%unique(FunDivEUROPE_trees_DE_SW_ES$plotcode),]

# add plot information to the tree-level dataset
tic()
FunDivEUROPE_trees_DE_SW_ES_Full <- left_join(FunDivEUROPE_trees_DE_SW_ES, FunDivEUROPE_plots_DE_SW_ES,
                                              by="plotcode")

toc()

# add species information to the tree-level dataset
tic()
FunDivEUROPE_trees_DE_SW_ES_Full <- left_join(FunDivEUROPE_trees_DE_SW_ES_Full, FunDivEUROPE_species,
                                              by="speciesid")

toc()


#### dataset for Finland and Belgium
FunDivEUROPE_plots_FI_WA <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_FI_WA.csv")
FunDivEUROPE_trees_FI_WA <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_trees_FI_WA.csv")

# add plot information to the tree-level dataset
FunDivEUROPE_trees_FI_WA_Full <- left_join(FunDivEUROPE_trees_FI_WA, FunDivEUROPE_plots_FI_WA,
                                           by="plotcode")
# add species information to the tree-level dataset
FunDivEUROPE_trees_FI_WA_Full <- left_join(FunDivEUROPE_trees_FI_WA_Full, FunDivEUROPE_species,
                                           by="speciesid")

# merge both dataset
FunDivEUROPE_trees_EU_Full <- rbind(FunDivEUROPE_trees_DE_SW_ES_Full, FunDivEUROPE_trees_FI_WA_Full)
saveRDS(FunDivEUROPE_trees_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_trees_EU_Full.rds")

# create a column of species full name
FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full %>% mutate(species_name=paste(genus, species))


setwd("D:/Zeus/ETH_zurich_MSc/ETHz_S2/Eco_Evo_ETH-D_assistantship/Leaf_Phenology_Project/Data") 
## load species and genus-level information of leaf phenology and leaf form  based on TRY dataset
load("Leaf_Phenology_Df.RData")
load("Leaf_Type_Df.RData")

FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full %>% mutate(Leaf_Phenology=NA, Leaf_Type=NA)

# get all unique tree species in FunDivEUROPE data
SppLst1 <- unique(FunDivEUROPE_trees_EU_Full$species_name)
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

## add all species/genus level information back to FunDivEUROPE dataset
for(i in 1:Len1){
  if(!is.na(LP_explore_Df$Leaf_Phenology[i])){
    FunDivEUROPE_trees_EU_Full$Leaf_Phenology[FunDivEUROPE_trees_EU_Full$species_name==LP_explore_Df$SppName[i]] <- LP_explore_Df$Leaf_Phenology[i]
  }
  if(!is.na(LT_explore_Df$Leaf_Type[i])){
    FunDivEUROPE_trees_EU_Full$Leaf_Type[FunDivEUROPE_trees_EU_Full$species_name==LT_explore_Df$SppName[i]] <- LT_explore_Df$Leaf_Type[i]
  }
}


FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full%>%mutate(BA=pi*(dbh2/2)^2, # basal area of each individual tree
                                                                  LP_known=ifelse(is.na(Leaf_Phenology),0,1), # binary value to show whether leaf phenology information is available for the tree
                                                                  EV_Status=ifelse(Leaf_Phenology=="Evergreen", 1, 0), # binary value to show whether the tree is evergreen or not
                                                                  DE_Status=ifelse(Leaf_Phenology=="Deciduous", 1, 0), # binary value to show whether the tree is deciduous or not
                                                                  BL_Status=ifelse(Leaf_Type=="Broadleaf", 1, 0), # binary value to show whether the tree is broadleaf or not
                                                                  NL_Status=ifelse(Leaf_Type=="Needleleaf", 1, 0)) # binary value to show whether the tree is needleleaf or not
saveRDS(FunDivEUROPE_trees_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_trees_EU_Full.rds")

# we define recruits of each leaf phenology type as all corresponding individuals in the current census, which were not present in the previous census.
FunDivEUROPE_trees_EU_Full$recruit[FunDivEUROPE_trees_EU_Full$dbh1==0 & FunDivEUROPE_trees_EU_Full$dbh2>0] <- 1
# compute recruits for both evergreen and deciduous trees
FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full%>% mutate(recruit.ev=recruit*EV_Status, recruit.de=recruit*DE_Status)
# we define recruits of each leaf phenology type as all corresponding individuals present in the previous census, which were not present in the current census.
FunDivEUROPE_trees_EU_Full$mortality[FunDivEUROPE_trees_EU_Full$dbh1>0 & FunDivEUROPE_trees_EU_Full$dbh2==0] <- 1
# compute the basal area of trees with known information of leaf phenology
FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full%>% mutate(BA_LP=BA*LP_known)

## function to get relative basal area of the tree species with the largest cumulative basal area in a plot
spp_dominance <- function(x){
  TT <- table(x)
  maxProp <- max(TT)/sum(TT)
  return(maxProp)
}


tic()
FunDivEUROPE_plots_EU_Full <- FunDivEUROPE_trees_EU_Full%>%group_by(plotcode)%>%summarise(
  latitude=mean(latitude, na.rm=T), # plot longitude
  longitude=mean(longitude,na.rm=T), # plot latitude
  spp_maxprop=spp_dominance(speciesid), # the relative basal area of the tree species with the largest cumulative basal area in the plot
  spp.count=length(unique(speciesid)), # number of species in the plot
  stem.density=length(speciesid),  # number of trees in the plot
  BA=sum(BA, na.rm=T), # plot total basal area
  BA_LP=sum(BA_LP, na.rm=T),  # plot basal area of all trees with known leaf phenology information
  BA_EV=sum(BA_EV, na.rm=T), # plot basal area of evergreen trees
  recruit.ev=sum(recruit.ev,na.rm=T), # total recruits of evergreen trees in the plot
  recruit.de=sum(recruit.de,na.rm=T) # total recruits of deciduous trees in the plot
)
toc()

saveRDS(FunDivEUROPE_plots_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_EU_Full.rds")

#  only keep plots for which information on leaf phenology strategies was available for > 90% of trees (weighted by basal area) 
FunDivEUROPE_plots_EU_Full <- FunDivEUROPE_plots_EU_Full[FunDivEUROPE_plots_EU_Full$BA_LP/FunDivEUROPE_plots_EU_Full$BA>0.9,]
# compute plot-level relative evergreen abundance
FunDivEUROPE_plots_EU_Full <- FunDivEUROPE_plots_EU_Full%>%mutate(relEV=BA_EV/BA_LP)%>%drop_na()

## get coordinates of FunDivEUROPE, which is later used to extract spatial environmental covariates from the lab composite
FunDivEUROPE_Coordinates <- FunDivEUROPE_plots_EU_Full%>% dplyr::select(plotcode, latitude, longitude) 
colnames(FunDivEUROPE_Coordinates) <- c("plotcode", "latitude", "longitude")
write.csv(FunDivEUROPE_Coordinates,row.names = F, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_Coordinates.csv")

# load spatial environmental covariates from the lab composite
FunDivEUROPE_covariates <- read.csv("**Continental_Analysis/FunDivEurope_EU/Data/20220925_FunDivEUROPE.csv")

tic()
# add spatial environmental covariates to plot-level FunDivEUROPE
FunDivEUROPE_plots_EU_Full <- left_join(FunDivEUROPE_plots_EU_Full,FunDivEUROPE_covariates, by="plotcode")
toc()
saveRDS(FunDivEUROPE_plots_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_EU_Full.rds")

tic()
# add spatial environmental covariates to tree-level FunDivEUROPE
FunDivEUROPE_trees_EU_Full <- left_join(FunDivEUROPE_trees_EU_Full,FunDivEUROPE_covariates, by="plotcode")
toc()

## get names of spatial environmental covariates
covariateList <- colnames(FunDivEUROPE_covariates)

# extract columns of spatial environmental covariates
composite <- FunDivEUROPE_plots_EU_Full[,covariateList]

tic()
## PCA across all spatial environmental covariates
env.pca <- prcomp(composite, scale = T)
## get the variation captured by the leading 10 PCs
var.10 <- sum(summary(env.pca)$importance[2,1:10]) * 100
cat(paste0(round(var.10,2),'% of variation capture in by the first 10 principle components./n'))
toc()

# add the leading 10 PCs to the original plot-level dataset
FunDivEUROPE_plots_EU_Full <- cbind(FunDivEUROPE_plots_EU_Full, as.data.frame(env.pca$x[,1:10]))

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
FunDivEUROPE_plots_EU_Full$Cluster <- raster::extract(fishNet, FunDivEUROPE_plots_EU_Full[,c("longitude", "latitude")])



FunDivEUROPE_trees_EU_Full <- FunDivEUROPE_trees_EU_Full[FunDivEUROPE_trees_EU_Full$plotcode%in%FunDivEUROPE_plots_EU_Full$plotcode,]
# add the 10 PCs to the tree-level dataset
FunDivEUROPE_trees_EU_Full <- left_join(FunDivEUROPE_trees_EU_Full, 
                                        FunDivEUROPE_plots_EU_Full[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")],
                                        by="plotcode")

saveRDS(FunDivEUROPE_trees_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_trees_EU_Full.rds")
saveRDS(FunDivEUROPE_plots_EU_Full, file="**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_EU_Full.rds")

