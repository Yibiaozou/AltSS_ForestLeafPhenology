library(tictoc)

forest.sim <- function(g.mod.de, g.mod.ev,
                       r.mod.de, r.mod.ev, 
                       m.mod.de, m.mod.ev, 
                       initial_density = 20, n.plots = 100, n.step = 20,
                       disturb_rate = 0.018,
                       step.switch = NA, switch.lev = NA, #if changing N level mid run.
                       env.cov = NA, n.cores = NA, silent = F,
                       LP.split = 'within_plot', split_frac = 0.5, mode='null'){ 

  # 
  #source other functions on which this depends.----
  source('Project_Functions/predict_gam_well.r')
  
  #Compatibility tests.----
  #check if mgcv forest is installed.
  if('mgcv' %in% rownames(installed.packages()) == F){
    stop('This function requires the mgcv package, please install it.\n')
  }
  if("mgcv" %in% (.packages()) == F){
    library(mgcv)
  }
  #Check inputs are comaptible.
  if(initial_density %% 2 != 0){
    stop('Initial stem density needs to be an even, integer value.\n')
  }
  if(n.plots         %% 2 != 0){
    stop('n.plots needs to be an even, integer value.\n')
  }
  if(is.data.frame(env.cov) == F){
    stop('env.cov must be a data frame, even if its a data frame with only one row.\n')
  }

  
  #build the initial tree and plot tables.----
  tree <- matrix(data = 12.7,nrow = initial_density, ncol = 1)
  colnames(tree) <- c('DIA.cm')
  tree <- data.frame(tree)
  if(LP.split == 'within_plot'){
    #Assign half trees ev, half de.
    tree$ev <- c(rep(0, initial_density/2), rep(1, initial_density/2))
    plot.list <- list()
    for(i in 1:n.plots){plot.list[[i]] <- tree}
  }
  if(LP.split == 'between_plot'){
    tree.1 <- tree
    tree.2 <- tree
    tree.1$ev <- 1
    tree.2$ev <- 0
    plot.list <- list()
    for(i in 1:n.plots){
      if(i <= n.plots * split_frac){
        plot.list[[i]] <- tree.1
        }
      if(i >  n.plots * split_frac){
        plot.list[[i]] <- tree.2
        }
    }
  }
  if(LP.split == 'uniform'){
    plot.list <- list()
    for(i in 1:n.plots){
      #draw relative abundance EV trees from uniform distribution.
      relEV <- runif(1,0,1)
      #calculate number of EV trees (rounded). Number of DE trees is just however trees remain post rounding.
      n.ev <- round(nrow(tree) * relEV)
      n.de <- nrow(tree) - n.ev
      #generate vector, drop in new tree table, 'tree.now'.
      ev <- c(rep(0, n.de), rep(1, n.ev))
      tree.now <- tree
      tree.now$ev <- ev
      plot.list[[i]] <- tree.now
    }
  }
  if(LP.split == 'all.ev'){
    plot.list <- list()
    for(i in 1:n.plots){
      ev             <- c(rep(1, nrow(tree)))
      tree.now       <- tree
      tree.now$ev    <- ev
      plot.list[[i]] <- tree.now
    }
  }
  if(LP.split == 'all.de'){
    plot.list <- list()
    for(i in 1:n.plots){
      ev             <- c(rep(0, nrow(tree)))
      tree.now       <- tree
      tree.now$ev    <- ev
      plot.list[[i]] <- tree.now
    }
  }
  
  #get plot table with plot level characteristics.
  plot.table <- list()
  for(i in 1:length(plot.list)){
    sum <- plot.list[[i]]
    density <- nrow(sum)
    plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
    plot.basal.ev <- sum(pi*((sum$ev*sum$DIA.cm)/2)^2)
    plot.basal.de <- plot.basal - plot.basal.ev
    ev.density <- sum(sum$ev)
    de.density <- length(sum$ev) - sum(sum$ev)
    relEV <- plot.basal.ev / plot.basal
    STDAGE <- 0
    return <- c(plot.basal, plot.basal.ev, plot.basal.de, density, de.density, ev.density, relEV, STDAGE)
    names(return) <- c('BASAL.plot','BASAL.ev','BASAL.de','stem.density','de.density','ev.density','relEV','STDAGE')
    #add the environmental covariates in (if you have any).
    if(sum(!is.na(env.cov)) > 0){
      #sample a row of the environmental covariate matrix.
      this.cov <- env.cov[sample(nrow(env.cov), 1),]
      return <- c(return, this.cov)
      return <- unlist(return)
    }
    plot.table[[i]] <- return
  }
  plot.table <- data.frame(do.call(rbind, plot.table))
  #track plot table through time in a list.
  super.table <- list(plot.table)
  #save a record of each plots environmental covariates.
  env.table <- plot.table[,colnames(plot.table) %in% colnames(env.cov)]
  
  #Begin simulation!----
  for(t in 1:n.step){
    tic()
    #1. Grow and kill your trees. Then recruit new trees.----
    new.plot.list <- list()
    
    for(j in 1:length(plot.list)){     # length(plot.list)
      #grab tree table for a given plot.
      if(nrow(plot.list[[j]])>0){
        cov <- plot.list[[j]]
        colnames(cov)[1] <- c('PREVDIA.cm')
        #merge plot-level covariates into tree table
        cov <- cbind(cov, plot.table[j,])
        cov <- as.data.frame(cov)
        if(is.na(cov$relEV[1])){
          cov$ev.density <- sum(cov$ev)
          cov$de.density <- nrow(cov) - sum(cov$ev)
          cov$BASAL.plot <- sum(pi*(cov$PREVDIA.cm/2)^2)
          cov$BASAL.ev <- sum(pi*(cov$PREVDIA.cm/2)^2*cov$ev)
          cov$BASAL.de <- cov$BASAL.plot - cov$BASAL.ev
          cov$relEV <-    cov$BASAL.ev/cov$BASAL.plot
          cov$stem.density <- nrow(cov)
        }
      }else{relEV <- mean(super.table[[t]]$relEV, na.rm=T)
      new.ev <- rbinom(n = initial_density,size = 1,prob = relEV)
      new.DIA.cm <- rep(12.7, initial_density)
      new.plot <- data.frame(new.DIA.cm, new.ev)
      colnames(new.plot) <- c("DIA.cm", "ev")
      cov <- new.plot
      colnames(cov)[1] <- c('PREVDIA.cm')

      plot.table[j,]$ev.density <- sum(cov$ev)
      plot.table[j,]$de.density <- nrow(cov) - sum(cov$ev)
      plot.table[j,]$BASAL.plot <- sum(pi*(cov$PREVDIA.cm/2)^2)
      plot.table[j,]$BASAL.ev <- sum(pi*(cov$PREVDIA.cm/2)^2*cov$ev)
      plot.table[j,]$BASAL.de <- plot.table[j,]$BASAL.plot[1] - plot.table[j,]$BASAL.ev[1]
      plot.table[j,]$relEV <-    plot.table[j,]$BASAL.ev/plot.table[j,]$BASAL.plot
      plot.table[j,]$stem.density <- nrow(cov)
      
      cov <- cbind(cov, plot.table[j,])
      }


      cov <- data.table(cov)
      #ORDER TREE TABLE BY ev STATUS.
      #So important otherwise you scramble ev status later down.
      cov <- cov[order(cov$ev),]
      
      #add Cluster ID that will be ignored to covariate table. Necessary to get around bam.predict() bug.
      check1 <- g.mod.de$model$Cluster
      check2 <- g.mod.ev$model$Cluster
      check3 <- r.mod.de$model$Cluster
      check4 <- r.mod.ev$model$Cluster
      check <- check1[check1 %in% check2]
      check <- check [check  %in% check3]
      check <- check [check  %in% check4]
      cov$Cluster <- check[1]
      
      #grow your trees.
      tree.new <- c()
      if(nrow(cov[cov$ev == 0,]) > 0){
        #tree.new.de <- predict(g.mod.de, newdata = cov[cov$ev == 0,], exclude = c("s(Cluster)","s(PLT_CN)"), newdata.guaranteed = T)
        tree.new.de <- predict_gam_well(g.mod.de, newdata=cov[cov$ev==0,], ranef.lab='Cluster')$fit
        tree.new    <- c(tree.new, tree.new.de)
      }
      if(nrow(cov[cov$ev == 1,]) > 0){
        #tree.new.ev <- predict(g.mod.ev, newdata = cov[cov$ev == 1,], exclude = c("s(Cluster)","s(PLT_CN)"), newdata.guaranteed = T)
        tree.new.ev <- predict_gam_well(g.mod.ev, newdata=cov[cov$ev==1,], ranef.lab='Cluster')$fit
        tree.new    <- c(tree.new, tree.new.ev)
      }
      tree.new <- data.frame(tree.new)
      
      #kill your trees.
      tree.dead <- c()
      if(nrow(cov[cov$ev == 0,]) > 0){
        #tree.dead.de <- predict(m.mod.de, newdata = cov[cov$ev == 0,], exclude = c("s(Cluster)","s(PLT_CN)"), newdata.guaranteed = T)
        tree.dead.de <- predict_gam_well(m.mod.de, newdata=cov[cov$ev==0,], ranef.lab='Cluster')$fit
        tree.dead    <- c(tree.dead, tree.dead.de)
      }
      if(nrow(cov[cov$ev == 1,]) > 0){
        #tree.dead.ev <- predict(m.mod.ev, newdata = cov[cov$ev == 1,], exclude = c("s(Cluster)","s(PLT_CN)"), newdata.guaranteed = T)
        tree.dead.ev <- predict_gam_well(m.mod.ev, newdata=cov[cov$ev==1,], ranef.lab='Cluster')$fit
        tree.dead    <- c(tree.dead, tree.dead.ev)
      }
      tree.dead    <- rbinom(length(tree.dead), 1, boot::inv.logit(tree.dead))   #logit, since model fit on logit scale.
      tree.new     <- data.frame(tree.new[!(tree.dead == 1),])                   #drop trees that died from tree table.
      tree.new$ev  <- cov$ev[!(tree.dead == 1)]                                  #insert ev status from covariate table.
      colnames(tree.new) <- c('DIA.cm','ev')
      
      #recruit new trees.
      r.dat <- plot.table[j,]
      r.dat$Cluster <- check[1]
      recruits.prob.de <- exp(predict_gam_well(r.mod.de, newdata=r.dat, ranef.lab='Cluster')$fit)
      recruits.prob.ev <- exp(predict_gam_well(r.mod.ev, newdata=r.dat, ranef.lab='Cluster')$fit)
      recruits.de      <- rpois(length(recruits.prob.de), recruits.prob.de)
      recruits.ev      <- rpois(length(recruits.prob.ev), recruits.prob.ev)
      #Hard limit on stem density.
      if(nrow(tree.new) > 120){
        recruits.de <- 0
        recruits.ev <- 0
      }
      
      #de recruits.
      if(!is.na(recruits.de) && recruits.de>0 && class(recruits.ev)=="integer"){
        new.recruit.de <- matrix(data = 12.7,nrow = recruits.de, ncol = 1)
        ev             <- matrix(data =    0,nrow = recruits.de, ncol = 1)
        colnames(new.recruit.de) <- 'DIA.cm'
        colnames(ev)             <- 'ev'
        new.recruit.de <- cbind(new.recruit.de, ev)
      }else{new.recruit.de <- NULL}

      #ev recruits.
      if(!is.na(recruits.ev) && recruits.ev>0 && class(recruits.ev)=="integer"){
        new.recruit.ev <- matrix(data = 12.7,nrow = recruits.ev, ncol = 1)
        ev             <- matrix(data =    1,nrow = recruits.ev, ncol = 1)
        colnames(new.recruit.ev) <- 'DIA.cm'
        colnames(ev)             <- 'ev'
        new.recruit.ev <- cbind(new.recruit.ev, ev)
      }else{new.recruit.ev <- NULL}

      
      #update your tree table.
      to_return <- rbind(tree.new, new.recruit.ev, new.recruit.de)
      
      #role stand replacing disturbance dice.
      annihilate <- rbinom(n = 1, size = 1, prob = disturb_rate)
      if(annihilate == 1){
        relEV <- sum(to_return[,2]) / nrow(to_return)
        new.ev <- rbinom(n = initial_density,size = 1,prob = relEV)
        new.DIA.cm <- rep(12.7, initial_density)
        new.plot <- data.frame(new.DIA.cm, new.ev)
        colnames(new.plot) <- colnames(to_return)
        to_return <- new.plot
      }
      
      #return result.
      # return(to_return)
      new.plot.list[[j]] <- to_return
    } #end parallel plot loop.
    plot.list <- new.plot.list
    
    #2. Update plot table.----
    plot.table <- list()
    for(i in 1:length(plot.list)){
      sum <- plot.list[[i]]
      density <- nrow(sum)
      ev.density <- sum(sum$ev)
      de.density <- length(sum$ev) - sum(sum$ev)
      plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
      plot.basal.ev <- sum(pi*((sum$ev*sum$DIA.cm)/2)^2)
      plot.basal.de <- plot.basal - plot.basal.ev
      relEV <- plot.basal.ev / plot.basal
      STDAGE <- t*5
      return <- c(plot.basal, plot.basal.ev, plot.basal.de, density, de.density, ev.density, relEV, STDAGE)
      names(return) <- c('BASAL.plot','BASAL.ev','BASAL.de','stem.density','de.density','ev.density','relEV','STDAGE')
      #add the static environmental covariates in (if you have any).
      if(sum(!is.na(env.cov)) > 0){
        this.env <- env.table[i,]
        return <- c(return, this.env)
        return <- unlist(return)
      }
      plot.table[[i]] <- return
    }
    plot.table <- data.frame(do.call(rbind, plot.table))
    #update super table.
    super.table[[t+1]] <- plot.table
    
    
    output_sub <- list(plot.table, super.table, env.table)
    names(output_sub) <- c('plot.table','super.table','env.table')

    
    toc()
    #3. report time step complete.----
    if(silent == F){
      current_time <- t*5
      talk <- paste0(current_time,' years of simulation complete.\n')
      cat(talk)
    }
  }
  #return simulation output.----
  output <- list(plot.table, super.table, env.table)
  names(output) <- c('plot.table','super.table','env.table')
  return(output)
}
