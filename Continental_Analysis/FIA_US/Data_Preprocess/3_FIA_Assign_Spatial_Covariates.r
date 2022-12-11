#PCA of environmental factors.

rm(list = ls())

storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')

library(factoextra)
library(data.table)
library(raster)
library(readr)
library(matrixStats)
library(tidyverse)


#load data.----
p1 <- readRDS(Product_1.path)
p2 <- readRDS(Product_2.path)

p1$PLT_CN <- as.numeric(p1$PLT_CN)
composite <- read.csv("D:/Zeus/ETH_zurich_MSc/ETHz_S2/Eco_Evo_ETH-D_assistantship/Leaf_Phenology_Project/Code/Niubi_transform/CovariatesExtracted_old/20210804_FIA_Merged_sampled_dataset_1000000.csv")

composite <- data.table(composite)

composite$PLT_CN <-as.numeric(composite$PLT_CN)
composite <- composite[composite$PLT_CN %in% p1$PLT_CN,]

plot.ref   <- composite[,c('PLT_CN','latitude','longitude')]


tic()
## PCA across all spatial environmental covariates
env.pca <- prcomp(composite, scale = T)
## get the variation captured by the leading 10 PCs
var.10 <- sum(summary(env.pca)$importance[2,1:10]) * 100
cat(paste0(round(var.10,2),'% of variation capture in by the first 10 principle components./n'))
toc()

# add the leading 10 PCs to the original dataset
p1 <- cbind(p1, as.data.frame(env.pca$x[,1:10]))



#run PCA on full composite. Takes ~20s.-----
tic()
env.pca <- prcomp(composite, scale = T)
var.10 <- sum(summary(env.pca)$importance[2,1:10]) * 100
cat(paste0(round(var.10,2),'% of variation capture in by the first 10 principle components./n'))
toc()

#grab first ten PCA values, append to PLT_CN codes w/ lat-lon. Also pull MAT and MAP.----
output <- cbind(plot.ref, 
                composite[,"CHELSA_Annual_Mean_Temperature"], 
                env.pca$x[,1:10])
#drop in mat independently, because it is useful for hysteresis simulation.
setnames(output,c('CHELSA_Annual_Mean_Temperature'),c('mat'))


#assign EPA ecoregions.----
ecoregions <- raster::shapefile(EPA_L2_ecoregions_raster.path)           #load shape file.
ecoregions <- spTransform(ecoregions, CRS('+proj=longlat +datum=WGS84')) #re-project raster to WGS84.
pts <- cbind(p1$longitude, p1$latitude)                          #bind up your points.
ecoreg.assign <- raster::extract(ecoregions, pts)                         #extract raster.
p1$ecoregion <- ecoreg.assign$NA_L2NAME

ecoreg.df <- as.data.frame(cbind(ecoregion=ecoreg.assign$NA_L2NAME, PLT_CN=p1$PLT_CN))
p2 <- merge(p2, ecoreg.df, all.x = T, by = 'PLT_CN')

#append environmental stats to p1 and p2, save.----
p1 <- merge(p1, output, all.x = T, by = 'PLT_CN')
p2 <- merge(p2, output, all.x = T, by = 'PLT_CN')
saveRDS(p1, Product_1.path)
saveRDS(p2, Product_2.path)
write.csv(comp.varnames, composite_variable_names.path)

