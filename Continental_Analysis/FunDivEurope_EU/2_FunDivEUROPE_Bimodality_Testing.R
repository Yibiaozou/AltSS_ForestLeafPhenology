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

####----Continental Analysis: Europe----####
FunDivEUROPE_plots_EU_Full <- readRDS("**Continental_Analysis/FunDivEurope_EU/Data/FunDivEUROPE_plots_EU_Full.rds")

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

saveRDS(b1zap, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_GAM_ZAP_EV_FSD.rds")
saveRDS(b2zap, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_GAM_ZAP_DE_FSD.rds")

iterations=1000
AzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b1zap)))
BzapMx <- matrix(NA, nrow=iterations, ncol=length(fitted(b2zap)))

for(i in 1:iterations){
  tic()
  AzapMx[i,] <- rZAP(rep(1,length(fitted(b1zap))),mu=fitted(b1zap),sigma=fitted(b1zap,what="sigma"))
  BzapMx[i,] <- rZAP(rep(1,length(fitted(b2zap))),mu=fitted(b2zap),sigma=fitted(b2zap,what="sigma"))
  toc()
}

saveRDS(AzapMx, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_zapMx_EV_FSD.rds")
saveRDS(BzapMx, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_zapMx_DE_FSD.rds")

###---GAMLSS zero-adjusted negative binomial (ZANBI) distribution

## All parameters of the distribution (mean parameter with log link, dispersion parameter with log link, zero probability parameter with logit link) 
# were modelled as a function of the ten environmental principal components + spatial cluster (as random effect).

tic()
b1zanusi <- gamlss(ev.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b1zanusi, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_GAM_ZANBI_EV_FSD.rds")

tic()
b2zanusi <- gamlss(de.density~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"), discrete=T,method="fREML"),
                   nu.formula =~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re"),discrete=T,method="fREML"),
                   sigma.formula=~ba(~s(PC1,k=kv)+s(PC2,k=kv)+s(PC3,k=kv)+s(PC4,k=kv)+s(PC5,k=kv)+s(PC6,k=kv)+s(PC7,k=kv)+s(PC8,k=kv)+s(PC9,k=kv)+s(PC10,k=kv)+s(Cluster,bs="re") ,discrete=T,method="fREML"),
                   family=ZANBI(),data=AX2, control=gamlss.control(c.crit = 0.05,n.cyc = 20))
toc()
saveRDS(b2zanusi, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_GAM_ZANBI_DE_FSD.rds")

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

saveRDS(AbnMx, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_bnMx_EV_FSD.rds")
saveRDS(BbnMx, file="**Continental_Analysis/FunDivEurope_EU/Data/EU_bnMx_DE_FSD.rds")

