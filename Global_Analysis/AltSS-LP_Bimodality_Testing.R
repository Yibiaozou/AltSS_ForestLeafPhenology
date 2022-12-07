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




####----Global Analysis----####

memory.limit(5e5)

## load plot-level aggregated forest data
GFBI_Df2_Aggregated_Full_2 <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_Df2_aggregated_New.rds")

## load cluster-level aggregated forest data
GFBI_Df3_Scaled_10min <-readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_Df3_Scaled_10min_FullStandardized.rds")
## get names of spatial environmental covariates
covariateList <- colnames(GFBI_Df3_scaled_10min_RF)
covariateList <- covariateList[-which(covariateList%in%c("BI", "latitude", "longitude", "Resolve_Biome",
                                                         "SG_Soil_pH_H2O_0_100cm", "Nitrogen", "cnRatio",
                                                         "Organic_Carbon", "Cation","Cluster_10min"))]

# function to get density of one variable along another variable, for visualization
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# compute plot-level relative evergreen abundance based on stem density
GFBI_Df2_Aggregated_Full_2 <- GFBI_Df2_Aggregated_Full_2 %>% mutate(relEV_Density=EV_Density/(EV_Density+DE_Density))

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2[GFBI_Df2_Aggregated_Full_2$spp_maxprop<=0.75 & GFBI_Df2_Aggregated_Full_2$count>=10 & GFBI_Df2_Aggregated_Full_2$spp.count>=2,]

GFBI_Df2_Aggregated_Full_2_filtered <- GFBI_Df2_Aggregated_Full_2_filtered%>%drop_na()%>%mutate(density=get_density(relEV, relEV_Density, n=100))

# extract columns of spatial environmental covariates
composite <- as.data.frame(GFBI_Df2_Aggregated_Full_2_filtered[,covariateList])
composite <- data.table(composite)

tic()
## PCA across all spatial environmental covariates
env.pca <- prcomp(composite, scale = T)
## get the variation captured by the leading 10 PCs
var.10 <- sum(summary(env.pca)$importance[2,1:10]) * 100
cat(paste0(round(var.10,2),'% of variation capture in by the first 10 principle components./n'))
fviz_screeplot(env.pca, addlabels = TRUE, ylim = c(0, 35)) #factoextra package.
toc()

# add the leading 10 PCs to the original dataset
GFBI_Df2_Aggregated_Full_2_filtered <- cbind(GFBI_Df2_Aggregated_Full_2_filtered, as.data.frame(env.pca$x[,1:10]))


saveRDS(GFBI_Df2_Aggregated_Full_2_filtered, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_Df2_Aggregated_Full_2_FSD_Fig1.rds")

# set seed for reproducibility 
set.seed(2022)


### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 1-degree grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1 # set the resolution of the template as 1 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot
GFBI_Df2_Aggregated_Full_2_filtered$Cluster <- raster::extract(fishNet, GFBI_Df2_Aggregated_Full_2_filtered[,c("longitude", "latitude")])


##------GAMLSS zero adjusted Poisson (ZAP) distribution
kv <- 5 # To account for the possibility of nonlinear effects, all predictors were fit as independent spline terms with at most 5 ‘knots’

## All parameters of the distribution (mean parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).

tic()
# model for evergreen stem density
b1zap <- gamlss(EV_Density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv),discrete=T,method="fREML"),
                family=ZAP(),data=GFBI_Df2_Aggregated_Full_2_filtered, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b1zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZAP_EV_FSD.rds")


tic()
# model for deciduous stem density
b2zap <- gamlss(DE_Density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv),discrete=T,method="fREML"),
                family=ZAP(),data=GFBI_Df2_Aggregated_Full_2_filtered, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b2zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZAP_DE_FSD.rds")


b1zap <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZAP_EV_FSD.rds")
b2zap <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZAP_DE_FSD.rds")

iterations=1000
AzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b1zap)))
BzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b2zap)))

### randomly sample evergreen and deciduous stem densities for each plot 
# based on plot-specific environmental conditions using the GAMLSS models for 1000 times
for(i in 1:iterations){
  tic()
  AzapMx[i,] <- rZAP(rep(1,length(fitted(b1zap))),mu=fitted(b1zap),sigma=fitted(b1zap,what="sigma"))
  BzapMx[i,] <- rZAP(rep(1,length(fitted(b2zap))),mu=fitted(b2zap),sigma=fitted(b2zap,what="sigma"))
  toc()
}

saveRDS(AzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_zapMx_EV_FSD.rds")
saveRDS(BzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_zapMx_DE_FSD.rds")

####----GAMLSS zero-adjusted negative binomial (ZANBI) distribution
GFBI_Df2_Aggregated_Full_2_filtered <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_Df2_Aggregated_Full_2_FSD_Fig1.rds")

# kv <- 5# To account for the possibility of nonlinear effects, all predictors were fit as independent spline terms with at most 5 ‘knots’

## All parameters of the distribution (mean parameter with log link, dispersion parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).
tic()
# model for evergreen stem density
b1zanusi <- gamlss(EV_Density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   family=ZANBI(),data=GFBI_Df2_Aggregated_Full_2_filtered, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b1zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZANBI_EV_FSD.rds")

tic()
# model for deciduous stem density
b2zanusi <- gamlss(DE_Density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   family=ZANBI(),data=GFBI_Df2_Aggregated_Full_2_filtered, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b2zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_GAM_ZANBI_DE_FSD.rds")

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

saveRDS(AbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_bnMx_EV_FSD.rds")
saveRDS(BbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/GFBI_bnMx_DE_FSD.rds")



####----Continental Analysis: US----####
setwd("D:/Zeus/ETH_zurich_MSc/ETHz_S2/Eco_Evo_ETH-D_assistantship/Leaf_Phenology_Project/Code/YibiaoZ_202x_LeafHabit_ColinA/altSS_forest_mycorrhizas-master/")
source('paths.r')
AX<-readRDS(Product_1.1.path)

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
# extract spatial cluster for each plot, which later used for 
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

####----Continental Analysis: Europe----####
FunDivEUROPE_plots_EU_Full <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/FunDivEurope/FunDivEurope_Inventory_DE_SW_ES_2022-09/FunDivEUROPE_plots_EU_Full.rds")

# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
AX <- FunDivEUROPE_plots_EU_Full[FunDivEUROPE_plots_EU_Full$spp_maxprop<=0.75 & FunDivEUROPE_plots_EU_Full$spp.count>=2 & FunDivEUROPE_plots_EU_Full$stem.density>=10,]

# compute evergreen stem density and deciduous stem density for each plot
AX <- AX %>% mutate(ev.density=round(relEV*stem.density), de.density=round((1-relEV)*stem.density))

# compute relative evergreen abundance weighted by stem density
AX <- AX %>% mutate(relEV_Density=ev.density/(ev.density+de.density))

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot, which later used for 
AX$Cluster <- raster::extract(fishNet, AX[,c("longitude", "latitude")])


# the data
AX2<-AX[,c("latitude","longitude","Cluster","ev.density","de.density","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","relEV_Density")]
AX2<-AX2[complete.cases(AX2),]

###---GAMLSS zero adjusted Poisson (ZAP) distribution
tic()
b1zap <- gamlss(ev.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                family=ZAP(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
#+s(Cluster,bs="re")
toc()

tic()
b2zap <- gamlss(de.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                family=ZAP(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()

saveRDS(b1zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_GAM_ZAP_EV_FSD.rds")
saveRDS(b2zap, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_GAM_ZAP_DE_FSD.rds")

iterations=1000
AzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b1zap)))
BzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b2zap)))

for(i in 1:iterations){
  tic()
  AzapMx[i,] <- rZAP(rep(1,length(fitted(b1zap))),mu=fitted(b1zap),sigma=fitted(b1zap,what="sigma"))
  BzapMx[i,] <- rZAP(rep(1,length(fitted(b2zap))),mu=fitted(b2zap),sigma=fitted(b2zap,what="sigma"))
  toc()
}

saveRDS(AzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_zapMx_EV_FSD.rds")
saveRDS(BzapMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_zapMx_DE_FSD.rds")

###---GAMLSS zero-adjusted negative binomial (ZANBI) distribution

## All parameters of the distribution (mean parameter with log link, dispersion parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).

tic()
b1zanusi <- gamlss(ev.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b1zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_GAM_ZANBI_EV_FSD.rds")

tic()
b2zanusi <- gamlss(de.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b2zanusi, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_GAM_ZANBI_DE_FSD.rds")

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

saveRDS(AbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_bnMx_EV_FSD.rds")
saveRDS(BbnMx, file="D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/GFBI/EU_bnMx_DE_FSD.rds")

