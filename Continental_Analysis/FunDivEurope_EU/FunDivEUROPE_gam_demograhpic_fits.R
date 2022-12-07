rm(list=ls())

library(data.table)
library(mgcv)
library(doParallel)
library(tidyverse)
library(LaplacesDemon)
library(raster)
library(tictoc)

# load the plot-level and tree-level FunDivEUROPE dataset
FunDivEUROPE_plots_EU_Full <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/FunDivEurope/FunDivEurope_Inventory_DE_SW_ES_2022-09/FunDivEUROPE_plots_EU_Full.rds")
FunDivEUROPE_trees_EU_Full <- readRDS("D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Data/FunDivEurope/FunDivEurope_Inventory_DE_SW_ES_2022-09/FunDivEUROPE_trees_EU_Full.rds")



# we accounted for the potential impact of human management on the 
# establishment of monocultures by only keeping plots 
# 1) with at least ten trees, 
# 2) with at least two species, 
# 3) in which the relative basal area of the tree species with the largest cumulative basal area in the plot was smaller than 75%
d1 <- FunDivEUROPE_plots_EU_Full[FunDivEUROPE_plots_EU_Full$spp_maxprop<=0.75 & FunDivEUROPE_plots_EU_Full$spp.count>=2 & FunDivEUROPE_plots_EU_Full$stem.density>=10,]

# subset the tree-level data to have the same plot-id as the plot-level dataset
d2 <- FunDivEUROPE_trees_EU_Full[FunDivEUROPE_trees_EU_Full$plotcode%in%d1$plotcode,]

# compute deciduous basal area of each plot
d1$BASAL.de <- d1$BA-d1$BA_EV
# change names of some columns
setnames(d1,'BA','BASAL.plot')
setnames(d1,'BA_EV','BASAL.ev')
setnames(d1,'BA_DE' ,'BASAL.de')
setnames(d2,'BA_EV','BASAL.ev')
setnames(d2,'BA_DE','BASAL.de')
setnames(d2,'dbh1' ,'PREVDIA.mm')
setnames(d2,'dbh2' ,'DIA.mm')

# add plot-level properties to tree-level data
d2 <- left_join(d2,d1[,c("plotcode","stem.density","relEV","BASAL.ev", "BASAL.plot")],by="plotcode")

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
d1$Cluster <- raster::extract(fishNet, d1[,c("longitude", "latitude")])
d2$Cluster <- raster::extract(fishNet, d2[,c("longitude", "latitude")])

output.path <- "D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Results/Analysis_Results/Demo_GAMs/demographic_fits_gam_separate_FunDivEurope_FSD.rds"

# set plot code and cluster as factor
d1$plotcode <- as.factor(d1$plotcode)
d2$plotcode <- as.factor(d2$plotcode)
d2$Cluster <- as.factor(d2$Cluster)
d1$Cluster <- as.factor(d1$Cluster)

# remove NA
d1 <- d1%>%drop_na()
d2 <- d2%>%drop_na()

# change units from mm to cm
d2$PREVDIA.cm <- d2$PREVDIA.mm/10
d2$DIA.cm <- d2$DIA.mm/10
d2$BASAL.plot <- d2$BASAL.plot/100
d1$BASAL.plot <- d1$BASAL.plot/100
d1$BASAL.ev <- d1$BASAL.ev/100

#Set global k.----
kk <- 5    #global max number of 'knots' to use in smoothing functions.
bs <- 'tp' #thin plate regression splines are default as well, but you can change things here if you like.
nt <- 8    #number of threads.

#Fit recruitment, mortality and growth models.----
#Environmental models with feedbacks.
#gam models.
cat('Fitting feedback models...\n');tic()
# evergreen recruitment
R.mod.ev <- bam(recruit.ev ~ s(relEV, k=kk, bs=bs) + s(BASAL.ev, k=kk, bs=bs) +  s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T)
toc()
tic()
# deciduous recruitment
R.mod.de <- bam(recruit.de ~ s(relEV, k=kk, bs=bs) + s(BASAL.de, k=kk, bs=bs) +  s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d1, family = 'poisson', nthreads=nt, discrete=T)
toc()
tic()
# evergreen mortality
M.mod.ev <- bam(mortality ~ s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$EV_Status == 1,], family = 'binomial', nthreads=nt, discrete=T, control=list(maxit = 500))
toc()
tic()
# deciduous mortality
M.mod.de <- bam(mortality ~ s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d2[d2$EV_Status == 0,], family = 'binomial', nthreads=nt, discrete=T, control=list(maxit = 500))
toc()
tic()
# evergreen growth
G.mod.ev <- bam(DIA.cm   ~ s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d3[d3$EV_Status == 1,], nthreads=nt, discrete=T)
toc()
tic()
# deciduous growth
G.mod.de <- bam(DIA.cm   ~ s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                  s(Cluster, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                data = d3[d3$EV_Status == 0,], nthreads=nt, discrete=T)
cat('Feedback models fit.\n');toc()

#wrap output and name.
y.feedback <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev)
names(y.feedback) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev')


#Grab plot environmental covariates for reference.----
#all plot-level environmental covariates. we will sample this later when simulating forests.
all.cov <- d1[,c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10', 'Soil_PH')]
#mean plot-level environmental covariates.
cov <- colMeans(all.cov, na.rm=T)
# cov <- all.cov[sample(1:nrow(all.cov), size=1),]
names(cov) <- colnames(all.cov)

#Save models and size categories.----
cat('Saving output...\n');tic()
output <- list(y.feedback, cov, all.cov)
names(output) <- c('y.feedback','env.cov','all.cov')

saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n');toc()