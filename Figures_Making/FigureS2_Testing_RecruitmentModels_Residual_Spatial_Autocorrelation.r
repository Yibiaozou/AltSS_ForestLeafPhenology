#compare spatial autocorrelation in DE vs. EV recruitment among models.
rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')

library(gstat)
library(sp)
library(mgcv)

#set output path.----
#load data and models.----
d1 <- readRDS(Product_1.path)
d2 <- readRDS(Product_2.path)

## load subset data for the whole US
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
fit <- readRDS(demographic_fits_gam_separate.path)


#Get DE and EV recruitment residuals of an intercept only model, and of the gams.----
raw.de.resid <- residuals(glm(d1$recruit.de ~ 1, family = poisson))
raw.ev.resid <- residuals(glm(d1$recruit.ev ~ 1, family = poisson))
fit.de.resid <- residuals(fit$y.feedback$R.mod.de)
fit.ev.resid <- residuals(fit$y.feedback$R.mod.ev)
resid.list <- list(raw.de.resid, raw.ev.resid, fit.de.resid, fit.ev.resid)


#fit variograms.----
#define spatial coordinates.
lat <- d1$latitude
lon <- d1$longitude

#fit you variograms.
vario.list <- list()
vario.model.list <- list()
for(i in 1:length(resid.list)){
  test.dat <- data.frame(resid.list[[i]], lat, lon)
  colnames(test.dat)[1] <- 'resid'
  coordinates(test.dat) <- ~ lon + lat
  vario       <-     variogram(resid~ 1, data = test.dat)           #get semivariance values
  vario.model <- fit.variogram(vario, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)
  vario.list[[i]] <- vario
  vario.model.list[[i]] <- vario.model
}

## Figure S2
#plot results.----
lab <- c('Raw DE Residuals','Raw EV Residuals','Fitted DE Residuals','Fitted EV Residuals')
par(mfrow = c(2,2))
for(i in 1:length(vario.list)){
  #plot semivariance as a function of distance.
  plot(vario.list[[i]]$gamma ~ vario.list[[i]]$dist, 
       bty = 'l', 
       ylim = c(0, max(vario.list[[i]]$gamma)*1.05), 
       pch = 16, 
       ylab = 'semivariance', xlab = 'distance (km)',
       main = lab[i])
  #calculate spherical fit line.
  par.mat <- vario.model.list[[i]]
  fit <- par.mat[1,2] + par.mat[2,2] * ((1.5 * (vario.list[[i]]$dist / par.mat[2,3])) - (0.5*(vario.list[[i]]$dist^3 / par.mat[2,3]^3)))
  fit <- ifelse(vario.list[[i]]$dist > par.mat[2,3], 
                par.mat[1,2] + par.mat[2,2],
                fit)
  #plot spherical fit line.
  lines(smooth.spline(fit ~ vario.list[[i]]$dist))
}

