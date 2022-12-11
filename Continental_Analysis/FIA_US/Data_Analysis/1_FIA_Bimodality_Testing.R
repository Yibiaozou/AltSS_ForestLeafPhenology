rm(list=ls())

library(tidyverse)
library(tictoc)
# library(geosphere)
library(mousetrap)
# library(diptest)
library(moments)
# library(raster)
# library(BimodalIndex)
library(LaplacesDemon)
library(MASS)
library(ggpubr)
# library(rgdal)
library(data.table)
# library("rnaturalearth")
# library("rnaturalearthdata")
# library(rgeos)
library(matrixStats)
library(factoextra)
library(ggfortify)
# library(caret)
# library(disttree)
# library(ranger)
library(gamlss)
library(gamlss.add)

####----Continental Analysis: US----####
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')
AX<-readRDS(Product_1.path)

# For FIA we already remove managed plot, therefore here we only remove small plots with less than 10 trees
AX <- AX[AX$count>=10,]
AX <- AX %>% mutate(relEV_Density=ev.density/(ev.density+de.density))
hist(AX$relEV_Density,breaks=20)
hist(AX$relEV,breaks=20)


# the data
AX2<-AX[,c("latitude","longitude","county.id","ev.density","de.density","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10", "relEV")]
AX2<-AX2[complete.cases(AX2),]

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot
AX2$Cluster <- raster::extract(fishNet, AX2[,c("longitude", "latitude")])

###---GAMLSS zero adjusted Poisson (ZAP) distribution
# kv <- 5

## All parameters of the distribution (mean parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).

tic()
# model for evergreen stem density
b1zap <- gamlss(ev.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                family=ZAP(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))

toc()

tic()
# model for deciduous stem density
b2zap <- gamlss(de.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                family=ZAP(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()

saveRDS(b1zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_GAM_ZAP_EV_FSD.rds")
saveRDS(b2zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_GAM_ZAP_DE_FSD.rds")

iterations=1000
AzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b1zap)))
BzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b2zap)))

### randomly sample evergreen and deciduous stem densities for each plot 
# based on plot-specific environmental conditions using the GAMLSS models for 1000 times
for(i in 1:iterations){
  tic()
  AzapMx[i,] <- rZAP(rep(1,length(fitted(b1zap))),mu=fitted(b1zap))
  BzapMx[i,] <- rZAP(rep(1,length(fitted(b2zap))),mu=fitted(b2zap))
  toc()
}

saveRDS(AzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_zapMx_EV_FSD.rds")
saveRDS(BzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_zapMx_DE_FSD.rds")

###---GAMLSS zero-adjusted negative binomial (ZANBI) distribution
# kv <- 5

## All parameters of the distribution (mean parameter with log link, dispersion parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).
tic()
b1zanusi <- gamlss(ev.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b1zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_GAM_ZANBI_EV_FSD.rds")

tic()
b2zanusi <- gamlss(de.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b2zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_GAM_ZANBI_DE_FSD.rds")

iterations=1000
AbnMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b1zanusi)))
BbnMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b2zanusi)))

### randomly sample evergreen and deciduous stem densities for each plot 
# based on plot-specific environmental conditions using the GAMLSS models for 1000 times
for(i in 1:iterations){
  tic()
  AbnMx[i,] <- rZANBI(rep(1,length(fitted(b1zanusi))),mu=fitted(b1zanusi),sigma=fitted(b1zanusi,what="sigma"),nu=fitted(b1zanusi,what="nu"))
  BbnMx[i,] <- rZANBI(rep(1,length(fitted(b2zanusi))),mu=fitted(b2zanusi),sigma=fitted(b2zanusi,what="sigma"),nu=fitted(b2zanusi,what="nu"))
  toc()
}

saveRDS(AbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_bnMx_EV_FSD.rds")
saveRDS(BbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/FIA_bnMx_DE_FSD.rds")
