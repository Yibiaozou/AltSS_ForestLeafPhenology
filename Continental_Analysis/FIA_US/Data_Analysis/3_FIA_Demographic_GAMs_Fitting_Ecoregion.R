#Fitting gams w/ separate AM-EV models and PC scores across ecoregions.
rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')

library(data.table)
library(mgcv)
library(doParallel)


#set output path.----
output.path <- demographic_fits_gam_ecoregion.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

# For FIA we already remove managed plot, therefore here we only remove small plots with less than 10 trees
d1 <- d1[d1$count>=10,]

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
# d2$county.ID <- as.factor(paste0(d2$STATECD,'_',d2$COUNTYCD))
# d1$county.ID <- as.factor(paste0(d1$STATECD,'_',d1$COUNTYCD))
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.Evergreen','BASAL.ev')
setnames(d1,'BASAL.Deciduous' ,'BASAL.de')

### Get spatial cluster of each plot. we partitioned the global forest zones using a ‘fishing net’ with 10 arc-min (~20km) grid size.
fishNet <- raster() # create an empty raster template 
res(fishNet) <- 1/6 # set the resolution of the template as 1/6 degree. 
values(fishNet) <- 1:ncell(fishNet) # allocate values to the 'fishing net' raster
# extract spatial cluster for each plot
d1$Cluster <- raster::extract(fishNet, d1[,c("longitude", "latitude")])
d2$Cluster <- raster::extract(fishNet, d2[,c("longitude", "latitude")])



#grab ecoregions of interest.----
ecoregion.ref <- table(d1$ecoregion)
ecoregion.ref <- names(ecoregion.ref[ecoregion.ref > 500])

#set k and basis for spline.----
kk <- 5    #number of knots in smooth functions.
bs <- 'tp' #thin plate regression are default, but you can change this here.
nt <- 8    #number of cores to use while fitting.

#Fit growth, recruitment and mortality models by ecoregion.----
output <- list()
tic()
for(i in 1:length(ecoregion.ref)){
  #subset data to ecoregion of interest.
  d1.sub <- d1[d1$ecoregion == ecoregion.ref[i],]
  d2.sub <- d2[d2$ecoregion == ecoregion.ref[i],]
  
  #fit gam models - no random effect.
  R.mod.ev <- bam(recruit.ev ~  s(relEV, k=kk, bs=bs) + s(BASAL.ev, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  R.mod.de <- bam(recruit.de ~  s(relEV, k=kk, bs=bs) + s(BASAL.de, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  M.mod.ev <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 1,], family = 'binomial', nthreads=nt, discrete=T)
  M.mod.de <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 0,], family = 'binomial', nthreads=nt, discrete=T)
  G.mod.ev <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 1,], nthreads=nt, discrete=T)
  G.mod.de <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 0,], nthreads=nt, discrete=T)
  
  #Grab plot environmental covariates for reference.
  #all plot-level environmental covariates. we will sample this later when simulating forests.
  all.cov <- d1.sub[,c('BASAL.plot','BASAL.de','BASAL.ev','stem.density','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  #mean plot-level environmental covariates.
  cov <- colMeans(all.cov)
  cov <- c(mean(d2.sub$PREVDIA.cm, na.rm=T),cov)
  names(cov)[1] <- 'PREVDIA.cm'
  
  #store in list and report.
  region.out <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev, cov)
  names(region.out) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev','cov')
  
  #fit gam models with cluster level random effect.----
  R.mod.ev <- bam(recruit.ev ~  s(relEV, k=kk, bs=bs) + s(BASAL.ev, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  R.mod.de <- bam(recruit.de ~  s(relEV, k=kk, bs=bs) + s(BASAL.de, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  M.mod.ev <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d2.sub[d2.sub$ev == 1,], family = 'binomial', nthreads=nt, discrete=T)
  M.mod.de <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d2.sub[d2.sub$ev == 0,], family = 'binomial', nthreads=nt, discrete=T)
  G.mod.ev <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d2.sub[d2.sub$ev == 1,], nthreads=nt, discrete=T)
  G.mod.de <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs)
                  + s(Cluster, bs = 're'), 
                  data = d2.sub[d2.sub$ev == 0,], nthreads=nt, discrete=T)
  
  region.out.clusterRE <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev, all.cov, cov)
  names(region.out.clusterRE) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev','all.cov','mean.cov')
  
  #fit gam models with plot level random effect.----
  M.mod.ev <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 1,], family = 'binomial', nthreads=nt, discrete=T)
  M.mod.de <- bam(mortality  ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 0,], family = 'binomial', nthreads=nt, discrete=T)
  G.mod.ev <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 1,], nthreads=nt, discrete=T)
  G.mod.de <- bam(DIA.cm   ~  s(relEV, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$ev == 0,], nthreads=nt, discrete=T)
  
  region.out.plotRE <- list(G.mod.de, G.mod.ev, M.mod.de, M.mod.ev, R.mod.de, R.mod.ev, cov)
  names(region.out.plotRE) <- c('G.mod.de','G.mod.ev','M.mod.de','M.mod.ev','R.mod.de','R.mod.ev','cov')
  
  total.out <- list(region.out, region.out.clusterRE, region.out.plotRE)
  names(total.out) <- c('none.re','cluster.re','plot.re')
  #total.out <- list(region.out, region.out.clusterRE)
  #names(total.out) <- c('none.re','cluster.re')
  output[[i]] <- total.out
    msg <- paste0(ecoregion.ref[i],' fit. ',i,' of ',length(ecoregion.ref),' ecoregions complete.\n')
  cat(msg); toc()
}
names(output) <- ecoregion.ref

#Save models and size categories.----
saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n')


