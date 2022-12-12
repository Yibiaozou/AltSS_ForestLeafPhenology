rm(list=ls())

storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/predict_gam_well.r')
library(mgcv)

#load ecoregion data.
d <- readRDS(demographic_fits_gam_ecoregion.path)
names(d) <- c('atlantic_highlands','mixed_wood_plains','mixed_wood_shield','appalachian_forests','southeastern_plains')
lab <- names(d)


#setup save paths.
out.dir <- '**/Continental_Analysis/FIA_US/Model_Ouput/'
cmd <- paste0('mkdir -p ',out.dir)
system(cmd)
output.paths <- paste0(out.dir,paste0(0,1:5),'_',lab,'.png')

for(i in 1:length(output.paths)){
  #grab ecoregion results of interest.----
         z <- d[[i]]$county.re
  path.out <- output.paths[i]
       env <- z$mean.cov
  #deal with county level random effect covariate, which will be ignored.
       check1 <- z$M.mod.de$model$county.ID
       check2 <- z$M.mod.ev$model$county.ID
       check3 <- z$R.mod.de$model$county.ID
       check4 <- z$R.mod.ev$model$county.ID
       check <- check1[check1 %in% check2]
       check <- check [check  %in% check3]
       check <- check [check  %in% check4]
       
  #setup covariates, generate predictions.----
  de.dat <- c(0,env)
  names(de.dat)[1] <- 'relEV'
  ev.dat <- de.dat
  ev.dat[1] <- 1
  dat <- data.frame(rbind(de.dat, ev.dat))
  dat$county.ID <- check[1]
  
  #make predictions on linear scale.
  de.mort <- predict_gam_well(z$M.mod.de, newdat=dat, ranef.lab='county.ID')
  ev.mort <- predict_gam_well(z$M.mod.ev, newdat=dat, ranef.lab='county.ID')
  de.recr <- predict_gam_well(z$R.mod.de, newdat=dat, ranef.lab='county.ID')
  ev.recr <- predict_gam_well(z$R.mod.ev, newdat=dat, ranef.lab='county.ID')
  N <- 300
  y.de.de.M <- rnorm(N, de.mort$fit[1], de.mort$se.fit[1])
  y.de.ev.M <- rnorm(N, de.mort$fit[2], de.mort$se.fit[2])
  y.ev.de.M <- rnorm(N, ev.mort$fit[1], ev.mort$se.fit[1])
  y.ev.ev.M <- rnorm(N, ev.mort$fit[2], ev.mort$se.fit[2])
  y.de.de.R <- rnorm(N, de.recr$fit[1], de.recr$se.fit[1])
  y.de.ev.R <- rnorm(N, de.recr$fit[2], de.recr$se.fit[2])
  y.ev.de.R <- rnorm(N, ev.recr$fit[1], ev.recr$se.fit[1])
  y.ev.ev.R <- rnorm(N, ev.recr$fit[2], ev.recr$se.fit[2])
  
  #back transform mortality and recruitment data.
  y.de.M <- boot::inv.logit(c(y.de.de.M,y.de.ev.M))
  y.ev.M <- boot::inv.logit(c(y.ev.de.M,y.ev.ev.M))
  y.de.R <-             exp(c(y.de.de.R,y.de.ev.R))
  y.ev.R <-             exp(c(y.ev.de.R,y.ev.ev.R))
  
  #convert mortality probability to survival probability.
  y.de.M <- 1 - y.de.M
  y.ev.M <- 1 - y.ev.M
  
  #setup x positions.
  jit <- 0.05
  x.pos <- c(0.5,1,1.5,2)
  x1 <- rnorm(N,x.pos[1],jit)
  x2 <- rnorm(N,x.pos[2],jit)
  x3 <- rnorm(N,x.pos[3],jit)
  x4 <- rnorm(N,x.pos[4],jit)
  x.de <- c(x1,x2)
  x.ev <- c(x3,x4)
  
  #Start plots.----
  png(path.out, width = 8, height = 3.5, units = 'in', res = 300)
  par(mfrow = c(1,2),
      mar = c(1,3.5,1,1))
  
  #Mortality plot.
  cols <- c('red','blue')
  col.1 <- c(rep(cols[1],N),rep(cols[2],N))
  limx <- c(min(c(x.de,x.ev)), max(c(x.de,x.ev)))
  limy <- c(min(c(y.de.M,y.ev.M))*0.98, max(c(y.de.M,y.ev.M))*1.02)
  #plot.
  plot(y.de.M ~ x.de, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
  points(y.ev.M ~ x.ev, pch = 16, cex = 0.3, col = col.1)
  abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
  #labels.
  mtext('Survival Probability', side = 2, cex = 1.0, line = 2.5)
  #legend('bottomleft',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
  #       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
  #       x.intersp = .75, xpd = T, 
  #       horiz = F)
  mtext('DE Forest',side = 3, line = -1, adj = 0.125, cex = 1)
  mtext('EV Forest',side = 3, line = -1, adj = 0.9, cex = 1)
  
  #Recruitment plot.
  limy <- c(0, max(c(y.de.R,y.ev.R))*1.05)
  #plot.
  plot(y.de.R ~ x.de, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
  points(y.ev.R ~ x.ev, pch = 16, cex = 0.3, col = col.1)
  abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
  #labels.
  lab <- expression(paste('Recruitment - Poisson ',lambda))
  mtext(lab, side = 2, cex = 1.0, line = 2.5)
  mtext('DE Forest',side = 3, line = -1, adj = 0.125, cex = 1)
  mtext('EV Forest',side = 3, line = -1, adj = 0.9, cex = 1)
  
  #end plots.-----
  dev.off()
  
}
