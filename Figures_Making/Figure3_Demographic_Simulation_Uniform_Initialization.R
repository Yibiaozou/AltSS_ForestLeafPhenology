#visualizing null vs. feedback demographic simulations.
#5-panel to show recruitment and mortality positive feedbacks, demo simulation altnerative stable states, and hysteresis.
rm(list=ls())
storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
library(mgcv)


#load demographic simulations.----
d <- readRDS(null_vs_feedback_simulation_output_RE_uniform.path)

y <- d$y.feedback.2$super.table[[401]] #2000 years in.
n <- d$n.feedback.2$super.table[[401]]


y <- rbind(y,y)
n <- rbind(n,n)
#set colors, calculate y limit.----
#colors.
cols <- c('#00acd9','#cfe83c') #pick colors
cols <- c('light green','light green')
trans <- 0.3 #set transparency

#number of breaks.
n.breaks <- 20

#generate colors.
rbPal <- colorRampPalette(c('red','blue'))
y$col <- rbPal(n.breaks)[as.numeric(cut(y$relEV,breaks = n.breaks))]
d.col <- y[order(y$relEV),]
# hist.cols <- unique(d.col$col)
hist.cols <- rbPal(n.breaks)
# hist.cols[n.breaks] <- "blue"

#y-limit.
y$cat <- cut(y$relEV, n.breaks)
count.y <- table(y$cat)
n$cat <- cut(n$relEV, n.breaks)
count.x <- table(n$cat)
limy <- c(0, max(c(count.x, count.y)))


#png save line.----
# png(output.path, width = 11, height = 5, units = 'in', res = 300)
par(mfrow = c(1,2))

# Figure 3A
#plot simulation without feedbacks.----
tab <- n
par( mar = c(4,4,1,2))
#plot line.
hist(tab$relEV, breaks = n.breaks, 
     xlim = c(0,1), ylim =  c(limy[1], limy[2]+40), 
     ylab = NA, xlab = NA, main = NA, 
     bty = 'l', col = hist.cols, lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance Evergreen Trees")), side = 1, line = 2.75, cex = 1.2)
mtext('Number of Forests', side = 2, line = 2.5, cex = 1.2)
# msg1 <- expression(paste('Simulation ',bolditalic('without')))
msg1 <- 'Null Simulation 1'
# msg2 <- 'con-LeafPhenology feedbacks'
mtext(msg1, side = 3, line = -1.5, adj = 0.05, cex = 1.2)
# mtext(msg2, side = 3, line = -4, adj = 0.40, cex = 1.2)
mtext('A', side = 3, line = -1.5, adj = 0.95, cex = 1, font = 2)

# Figure 3B
#plot simulation with feedbacks.----
tab <- y
#plot line.
hist(tab$relEV, breaks = n.breaks, 
     xlim = c(0,1), ylim = c(limy[1], limy[2]+40), #limy[2]+40
     ylab = NA, xlab = NA, main = NA, 
     bty = 'l', col = hist.cols, lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance Evergreen Trees")), side = 1, line = 2.75, cex = 1.2)
# mtext('Number of Forests', side = 2, line = 2.5, cex = 1)
# msg1 <- expression(paste('Simulation ',bolditalic('with')))
msg1 <- 'Feedback Simulation 1'
# msg2 <- 'con-LeafPhenology feedbacks'
mtext(msg1, side = 3, line = -1.5, adj = 0.05, cex = 1.2)
# mtext(msg2, side = 3, line = -4, adj = 0.40, cex = 1.2)
mtext('B', side = 3, line = -1.5, adj = 0.95, cex = 1, font = 2)

#end plot.----

