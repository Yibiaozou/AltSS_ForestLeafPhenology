#DE - EV hysteresis conceptual figure and simulation result.
#code below flips de vs. ev ramps, because i flipped the initial conditions in the actual simulation script.
rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')


#hysteresis data load and workup feedback results.----
#load data.
d <- readRDS(initial_condition_hysteresis_simulation.path)
#grab results after 2000 years, make these final plot tables to work with other code.
time.step <- 401 #400 time steps = 2000 years.
for(i in 1:(length(d$all.ev$nul)-1)){
  d$all.ev$nul    [[i]]$plot.table <- d$all.ev$nul    [[i]]$super.table[[time.step]]
  d$all.ev$alt.GRM[[i]]$plot.table <- d$all.ev$alt.GRM[[i]]$super.table[[time.step]]
  d$all.de$nul    [[i]]$plot.table <- d$all.de$nul    [[i]]$super.table[[time.step]]
  d$all.de$alt.GRM[[i]]$plot.table <- d$all.de$alt.GRM[[i]]$super.table[[time.step]]
}
n.tot <- nrow(d$all.ev$alt.GRM$l1$plot.table)
n.lev <- length(d$all.de$alt.GRM)

#get number of EV and DE plots in feedback simulations (>90%) and bootstrap CI.
#starting from EV state.
ev.ramp <- list()
for(i in 1:(length(d$all.ev$alt.GRM)-1)){
  plot.tab <- d$all.ev$alt.GRM[[i]]$plot.table
  N.ev     <- nrow(plot.tab[plot.tab$relEV > 0.9,])
  N.de     <- nrow(plot.tab[plot.tab$relEV < 0.1,])
  relEV    <- mean(plot.tab$relEV)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.ev    <- nrow(dat[dat$relEV > 0.9,])
    boot.de    <- nrow(dat[dat$relEV < 0.1,])
    boot.relEV <- mean(dat$relEV)
    boot.dat[[k]] <- c(boot.relEV,boot.ev,boot.de)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEV.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  ev.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  de.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  
  #return output.
  ev.ramp[[i]] <- c(relEV,N.ev,N.de,relEV.95,ev.95,de.95)
}
ev.ramp <- data.frame(do.call(rbind, ev.ramp))
colnames(ev.ramp) <- c('relEV','N.ev','N.de','relEV.lo95','relEV.hi95','ev.lo95','ev.hi95','de.lo95','de.hi95')

#starting from DE state.
de.ramp <- list()
for(i in 1:(length(d$all.de$alt.GRM)-1)){
  plot.tab <- d$all.de$alt.GRM[[i]]$plot.table
  N.ev     <- nrow(plot.tab[plot.tab$relEV > 0.9,])
  N.de     <- nrow(plot.tab[plot.tab$relEV < 0.1,])
  relEV    <- mean(plot.tab$relEV)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.ev    <- nrow(dat[dat$relEV > 0.9,])
    boot.de    <- nrow(dat[dat$relEV < 0.1,])
    boot.relEV <- mean(dat$relEV)
    boot.dat[[k]] <- c(boot.relEV,boot.ev,boot.de)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEV.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  ev.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  de.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  #return output.
  de.ramp[[i]] <- c(relEV,N.ev,N.de,relEV.95,ev.95,de.95)
}
de.ramp <- data.frame(do.call(rbind, de.ramp))
colnames(de.ramp) <- c('relEV','N.ev','N.de','relEV.lo95','relEV.hi95','ev.lo95','ev.hi95','de.lo95','de.hi95')

#Use average relative abundance EV trees instead.
up <- ev.ramp[,c('relEV','relEV.lo95','relEV.hi95')]
down <- de.ramp[,c('relEV','relEV.lo95','relEV.hi95')]
colnames(  up) <- c('n','lo95','hi95')
colnames(down) <- c('n','lo95','hi95')

#Assign MAT.
up  $MAT <- seq(-2, 23, length=13)
down$MAT <- seq(-2, 23, length=13)
up.feed <- up
down.feed <- down


#workup null results.----
#get number of EV and DE plots in feedback simulations (>90%) and bootstrap CI.
#starting from EV state.
ev.ramp <- list()
for(i in 1:(length(d$all.ev$nul)-1)){
  plot.tab <- d$all.ev$nul[[i]]$plot.table
  N.ev     <- nrow(plot.tab[plot.tab$relEV > 0.9,])
  N.de     <- nrow(plot.tab[plot.tab$relEV < 0.1,])
  relEV    <- mean(plot.tab$relEV)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.ev    <- nrow(dat[dat$relEV > 0.9,])
    boot.de    <- nrow(dat[dat$relEV < 0.1,])
    boot.relEV <- mean(dat$relEV)
    boot.dat[[k]] <- c(boot.relEV,boot.ev,boot.de)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEV.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  ev.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  de.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  
  #return output.
  ev.ramp[[i]] <- c(relEV,N.ev,N.de,relEV.95,ev.95,de.95)
}
ev.ramp <- data.frame(do.call(rbind, ev.ramp))
colnames(ev.ramp) <- c('relEV','N.ev','N.de','relEV.lo95','relEV.hi95','ev.lo95','ev.hi95','de.lo95','de.hi95')

#starting from DE state.
de.ramp <- list()
for(i in 1:(length(d$all.de$nul)-1)){
  plot.tab <- d$all.de$nul[[i]]$plot.table
  N.ev     <- nrow(plot.tab[plot.tab$relEV > 0.9,])
  N.de     <- nrow(plot.tab[plot.tab$relEV < 0.1,])
  relEV    <- mean(plot.tab$relEV)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.ev    <- nrow(dat[dat$relEV > 0.9,])
    boot.de    <- nrow(dat[dat$relEV < 0.1,])
    boot.relEV <- mean(dat$relEV)
    boot.dat[[k]] <- c(boot.relEV,boot.ev,boot.de)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEV.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  ev.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  de.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  #return output.
  de.ramp[[i]] <- c(relEV,N.ev,N.de,relEV.95,ev.95,de.95)
}
de.ramp <- data.frame(do.call(rbind, de.ramp))
colnames(de.ramp) <- c('relEV','N.ev','N.de','relEV.lo95','relEV.hi95','ev.lo95','ev.hi95','de.lo95','de.hi95')

#Use average relative abundance EV trees instead.
up.null <- ev.ramp[,c('relEV','relEV.lo95','relEV.hi95')]
down.null <- de.ramp[,c('relEV','relEV.lo95','relEV.hi95')]
colnames(  up.null) <- c('n','lo95','hi95')
colnames(down.null) <- c('n','lo95','hi95')

#Assign MAT.
up.null$MAT <- seq(-2, 23, length=13)
down.null$MAT <- seq(-2, 23, length=13)
# down.null <- down.null[-8,]

# Figure 3E
#Panel 1.  null simulation ramp up / ramp down.----
up <- up.null[1:13,]
down <- down.null[1:13,]
par(mfrow = c(1,2))
par(mar = c(4,4,1,2))
cols <- c('blue','red')
trans <- 0.3
#plot ramp up, transitioning away from EV dominated forests.
color <- cols[1]
max.y <- max(c(up.null$hi95,up.feed$hi95)) * 1.07
limy <- c(0,max.y)
#plot(up$n ~ up$MAT, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
plot(up$n ~ up$MAT, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = limy)
lines(smooth.spline(up$n ~ up$MAT, spar = .1), lwd = 3, col = color)
polygon(c(up$MAT, rev(up$MAT)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$MAT, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$MAT, rev(down$MAT)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#legend.
legend(x=-2.5,y= 0.2,legend=c("initially EV","initially DE"),lty=1,lwd=2, col = cols, box.lwd=0,
       bg="transparent",cex= 1, bty = "n")
# 
#outer labels.
mtext('Relative Abundance Evergreen Trees', side =2, cex=1.2, line = 2.5)
mtext('Null Simulation 3', side = 3, adj = 0.05, line = -1.5, cex = 1.2)
mtext('Mean Annual Temperature (°C)', side = 1, line = 2.75, cex = 1.2)
mtext('E', side = 3, line = -1.5, adj = 0.95, font = 2)
# mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 3.5,  cex = 0.7)

# Figure 3F
#Panel 2. hysteresis feedback simulation ramp up / ramp down.----
up <-   up.feed[1:13,]
down <- down.feed[1:13,]
# par(mar = c(4,4,1,2))
cols <- c('blue','red')
trans <- 0.3
#plot ramp up, transitioning away from EV dominated forests.
color <- cols[1]
plot(up$n ~ up$MAT, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = limy)
lines(smooth.spline(up$n ~ up$MAT, spar = .1), lwd = 3, col = color)
polygon(c(up$MAT, rev(up$MAT)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$MAT, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$MAT, rev(down$MAT)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#outer labels.
#mtext('Number EV Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('Feedback Simulation 3', side = 3, adj = 0.05, line = -1.5, cex = 1.2)
mtext('Mean Annual Temperature (°C)', side = 1, line = 2.75, cex = 1.2)
mtext('F', side = 3, line = -1.5, adj = 0.95, font = 2)
# mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 3.5,  cex = 0.7)

#end plot.----


