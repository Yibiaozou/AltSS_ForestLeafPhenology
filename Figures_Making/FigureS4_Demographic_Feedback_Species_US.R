rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/predict_gam_well.r')
library(mgcv)
library(data.table)
library(tidyverse)

#load species level recruitment and mortality data.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
spp.id <- 1:20
d$spp.key <- NULL
g.mu <- rep(0, length(spp.id))
g.sd <- rep(0, length(spp.id))
r.mu <- rep(0, length(spp.id))
r.sd <- rep(0, length(spp.id))
m.mu <- rep(0, length(spp.id))
m.sd <- rep(0, length(spp.id))
j <- 1
for(i in spp.id){
  g.mu[j] <- d[[i]]$post.contrast$g.mu
  g.sd[j] <- d[[i]]$post.contrast$g.sd
  m.mu[j] <- d[[i]]$post.contrast$m.mu
  m.sd[j] <- d[[i]]$post.contrast$m.sd
  r.mu[j] <- d[[i]]$post.contrast$r.mu
  r.sd[j] <- d[[i]]$post.contrast$r.sd
  j <- j + 1
}
spp.key$g.mu <- g.mu
spp.key$g.sd <- g.sd
spp.key$m.mu <- m.mu
spp.key$m.sd <- m.sd
spp.key$r.mu <- r.mu
spp.key$r.sd <- r.sd

#setup pch key.
spp.key$pch <- 16
spp.key$pch <- ifelse(spp.key$Leaf_Phenology == 'Evergreen' & spp.key$gymno == 1, 17, spp.key$pch)
spp.key$pch <- ifelse(spp.key$Leaf_Phenology ==  'Deciduous'                     ,  1, spp.key$pch)
spp.key$pch <- ifelse(spp.key$Leaf_Phenology ==  'Deciduous' & spp.key$gymno == 1,  2, spp.key$pch)


#setup x positions for first two panels.-----
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.de <- c(x1,x2)
x.ev <- c(x3,x4)


t.test(y.de.surv[1:N], y.de.surv[(N+1):(2*N)])
t.test(y.ev.surv[1:N], y.ev.surv[(N+1):(2*N)])
t.test(y.de.recr[1:N], y.de.recr[(N+1):(2*N)])
t.test(y.ev.recr[1:N], y.ev.recr[(N+1):(2*N)])
t.test(y.de.grow[1:N], y.de.grow[(N+1):(2*N)])
t.test(y.ev.grow[1:N], y.ev.grow[(N+1):(2*N)])

#Setup panels.----
par(mfrow = c(2,2),
    mar = c(1.5,6.5,5,1),
    oma = c(1,1, 1,1))

#plot survival effects.----
cols <- c('red','blue')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
limx <- c(min(c(x.de,x.ev)), max(c(x.de,x.ev)))
limy <- c(min(c(y.de.surv,y.ev.surv)), 1)
#plot.
plot(y.de.surv ~ x.de, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.ev.surv ~ x.ev, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
mtext('Survival Probability', side = 2, cex = 1.0, line = 2.5)
legend('topleft',legend = c('All EV trees','All DE trees'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('DE Forest',side = 1, line = 1, adj = 0.175, cex = 1)
mtext('EV Forest',side = 1, line = 1, adj = 0.825, cex = 1)
mtext('a.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)


#species level survival.-----
spp.key <- spp.key[order(spp.key$m.mu),]
x    <- spp.key$m.mu
x.sd <- spp.key$m.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))), na.rm=T)
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Survival Advantage\namong DE Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.7)
mtext('Survival Advantage\namong EV Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.7)
mtext('b.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)
#legend
legend(x = 1.5, y = 8, legend = c('DE broadleaf','DE needleleaf', 'EV broadleaf','EV needleleaf'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)

#plot recruitment effects.----
cols <- c('red','blue')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
labx <- c('DE in EV Forest','DE in DE Forest','EV in DE Forest','EV in EV Forest')
limx <- c(min(c(x.de,x.ev)), max(c(x.de,x.ev)))
limy <- c(0, max(c(y.de.recr,y.ev.recr)))
#plot.
plot(y.de.recr ~ x.de, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.ev.recr ~ x.ev, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
lab <- expression(paste('Recruitment - Poisson ',lambda))
mtext(lab, side = 2, cex = 1.0, line = 2.5)
mtext('DE Forest',side = 1, line = 1, adj = 0.175, cex = 1)
mtext('EV Forest',side = 1, line = 1, adj = 0.825, cex = 1)
mtext('b.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)

#species level recruitment.----
spp.key <- spp.key[order(spp.key$r.mu),]
x    <- spp.key$r.mu
x.sd <- spp.key$r.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))))
# limx <- c(-limx, limx)
limx <- c(-45, 45)
plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Recruitment Advantage\namong DE Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.7)
mtext('Recruitment Advantage\namong EV Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.7)
mtext('d.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)

