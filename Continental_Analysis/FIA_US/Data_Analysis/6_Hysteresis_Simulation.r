#running hysteresis simulations w/ stand replacing disturbance rate.
rm(list=ls())

storage.dir <- '**/Continental_Analysis/FIA_US/'
setwd(storage.dir)
source('Project_Functions/gam.int_forest.sim.r')

library(mgcv)
library(data.table)
library(tidyverse)
library(tictoc)


#set output path.----
output.path <- initial_condition_hysteresis_simulation.path

#load models and environmental covariates.----
fits <- readRDS(demographic_fits_gam_separate.path)
env.cov <- fits$all.cov
env.cov <- env.cov %>% drop_na()
#register parallel environment.----
n.cores <- 8

#Specify MAT up and down ramp ranges and number of plots.----
MAT.range <- seq(-2,23, length=13)
N.PLOTS <- 1000 #Must be even!
N.STEPS <- 400 #400 steps = 2000 years, longer to make sure it runs to something "stable".

#Run ramp up models.----
cat('Running all EV simulations...\n');tic()

#Run all DE simulatic() #start timer.
for(i in 1:length(MAT.range)){
  env.cov$mat <- MAT.range[i]
  # Null model.
  ev.nul[[i]] <- forest.sim(g.mod.de = fits$n.feedback$G.mod.de,
                            g.mod.ev = fits$n.feedback$G.mod.ev,
                            m.mod.de = fits$n.feedback$M.mod.de,
                            m.mod.ev = fits$n.feedback$M.mod.ev,
                            r.mod.de = fits$n.feedback$R.mod.de,
                            r.mod.ev = fits$n.feedback$R.mod.ev,
                            env.cov = env.cov,
                            LP.split = 'between_plot',
                            split_frac = 0.8,
                            silent = T,
                            n.plots = N.PLOTS,
                            n.cores = n.cores,
                            n.step = N.STEPS, mode='null')
  
    
  #Feedback model.
  ev.alt[[i]] <- forest.sim(g.mod.de = fits$y.feedback$G.mod.de,
                   g.mod.ev = fits$y.feedback$G.mod.ev,
                   m.mod.de = fits$y.feedback$M.mod.de,
                   m.mod.ev = fits$y.feedback$M.mod.ev,
                   r.mod.de = fits$y.feedback$R.mod.de, 
                   r.mod.ev = fits$y.feedback$R.mod.ev,
                   env.cov = env.cov, 
                   LP.split = 'between_plot', split_frac = 0.8, 
                   n.plots = N.PLOTS,
                   n.cores = n.cores,
                   n.step = N.STEPS, mode='feedback')

  msg <- paste0(i,' of ',length(PH.ramp.range),' all EV simulations complete. ')
  cat(msg);toc()
}

save(ev.nul, file="ev_nul_PH_FSD.RData")
save(ev.alt, file="ev_alt_PH_FSD.RData")


#reset env.cov levels.
env.cov <- fits$all.cov
de.nul <- list()
de.alt <- list()


cat('Running all DE simulations...\n');tic() #start timer.
for(i in 1:length(MAT.range)){
  env.cov$mat <- MAT.range[i]
  #Null model.
  de.nul[[i]] <- 
    forest.sim(g.mod.de = fits$n.feedback$G.mod.de,
               g.mod.ev = fits$n.feedback$G.mod.ev,
               m.mod.de = fits$n.feedback$M.mod.de,
               m.mod.ev = fits$n.feedback$M.mod.ev,
               r.mod.de = fits$n.feedback$R.mod.de, 
               r.mod.ev = fits$n.feedback$R.mod.ev,
               env.cov = env.cov, 
               LP.split = 'between_plot', split_frac = 0.2, 
               n.plots = N.PLOTS,
               n.cores = n.cores,
               n.step = N.STEPS, mode='null')
  
  #Feedback model.
  de.alt[[i]] <- 
    forest.sim(g.mod.de = fits$y.feedback$G.mod.de,
               g.mod.ev = fits$y.feedback$G.mod.ev,
               m.mod.de = fits$y.feedback$M.mod.de,
               m.mod.ev = fits$y.feedback$M.mod.ev,
               r.mod.de = fits$y.feedback$R.mod.de, 
               r.mod.ev = fits$y.feedback$R.mod.ev,
               env.cov = env.cov, 
               LP.split = 'between_plot', split_frac = 0.2, 
               n.plots = N.PLOTS,
               n.cores = n.cores,
               n.step = N.STEPS, mode='feedback')

  msg <- paste0(i,' of ',length(MAT.range),' all DE simulations complete. ')
  cat(msg);toc()
}

#wrap, name and return output.----
cat('Wrapping output and saving...\n')
all.ev <- list(ev.nul, ev.alt)
all.de <- list(de.nul, de.alt)
names(all.ev) <- c('nul','alt.GRM')
names(all.de) <- c('nul','alt.GRM')


lab <- paste0('l',MAT.range)
for(i in 1:length(all.ev)){names(all.ev[[i]]) <- lab}
for(i in 1:length(all.de)){names(all.de[[i]]) <- lab}

output <- list(all.ev, all.de)
names(output) <- c('all.ev','all.de')
saveRDS(output, output.path) 
cat('Script complete. ');toc()
