#Ttesting demographic simulation model. 
rm(list=ls())

storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/paths.r')
source('Project_Functions/tic_toc.r')
source('Project_Functions/gam.int_forest.sim.r')
source('Project_Functions/predict_gam_well.r')


library(randomForest)
library(mgcv)
library(rowr)
library(data.table)
library(tidyverse)

# library(rowr)

#set output path.----
output.path.1 <- null_vs_feedback_simulation_output_RE.path
output.path.2 <- null_vs_feedback_simulation_output_RE_uniform.path

#load gam model results.----

d <- readRDS(demographic_fits_gam_separate.path)

envCov  <- d$all.cov%>%drop_na()

##---Simulations with bimodal initialization----
#Just run the function.
tic()
null <- forest.sim(g.mod.de = d$n.feedback$G.mod.de, g.mod.ev = d$n.feedback$G.mod.ev,
                   r.mod.de = d$n.feedback$R.mod.de, r.mod.ev = d$n.feedback$R.mod.ev,
                   m.mod.de = d$n.feedback$M.mod.de, m.mod.ev = d$n.feedback$M.mod.ev,
                   LP.split = 'between_plot',
                   env.cov = envCov,
                   n.cores = 8,
                   n.plots = 1000, n.step = 400, mode="null")
cat('Null simulation complete.\n')
toc()


tic()
feed <- forest.sim(g.mod.de = d$y.feedback$G.mod.de, g.mod.ev = d$y.feedback$G.mod.ev,
                   r.mod.de = d$y.feedback$R.mod.de, r.mod.ev = d$y.feedback$R.mod.ev,
                   m.mod.de = d$y.feedback$M.mod.de, m.mod.ev = d$y.feedback$M.mod.ev,
                   LP.split = 'between_plot',
                   env.cov = envCov,
                   n.cores = 8,
                   n.plots = 1000, n.step = 400, mode="feedback")
cat('Feedback simulation complete.\n')
toc()

#save output.
out <- list(null,feed)
names(out) <- c('n.feedback','y.feedback')
saveRDS(out, output.path.1)


##---Simulations with uniform initialization----
null_2 <- forest.sim(g.mod.de = d$n.feedback$G.mod.de, g.mod.ev = d$n.feedback$G.mod.ev,
                     r.mod.de = d$n.feedback$R.mod.de, r.mod.ev = d$n.feedback$R.mod.ev,
                     m.mod.de = d$n.feedback$M.mod.de, m.mod.ev = d$n.feedback$M.mod.ev,
                     LP.split = 'uniform',
                     env.cov = envCov,
                     n.cores = 8,
                     n.plots = 1000, n.step = 400, mode="null")
cat('Feedback simulation complete.\n')
toc()

tic()
feed_2 <- forest.sim(g.mod.de = d$y.feedback$G.mod.de, g.mod.ev = d$y.feedback$G.mod.ev,
                   r.mod.de = d$y.feedback$R.mod.de, r.mod.ev = d$y.feedback$R.mod.ev,
                   m.mod.de = d$y.feedback$M.mod.de, m.mod.ev = d$y.feedback$M.mod.ev,
                   LP.split = 'uniform',
                   env.cov = envCov,
                   n.cores = 8,
                   n.plots = 1000, n.step = 400, mode="feedback")
cat('Feedback simulation complete.\n')
toc()



#save output.
out_2 <- list(null_2,feed_2)
names(out_2) <- c('n.feedback','y.feedback')
saveRDS(out_2, output.path.2)

