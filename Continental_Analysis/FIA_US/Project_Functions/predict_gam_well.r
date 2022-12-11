#' predict_gam_well
#' I built this because the gam.predict methods for bams fit with discretization is broken in mgcv.
#' This could stand to have one millions checks, warnings and error messages built into it. But, academia.
#'
#' @param fit          fitted gam object.
#' @param newdata      new observations for which you want to make predictions.
#' @param ranef.lab    name of the random effect you want to exclude from predictions.
#' @param N.posterior  number of posterior draws from which to calulate prediction mean uncertainty.
#'
#' @return             returns a list with the vector "fit"  predicted mean values for each observation,
#'                     and "se.fit" which is the standard deviation of the predicted mean value.
#' @export
#'
#' @examples
predict_gam_well <- function(fit, newdata, ranef.lab = NA, N.posterior = 1000){
  #function to draw from multivariate normal.
  rmvn <- function(n,mu,sig) { ## MVN random deviates
    L <- mroot(sig);m <- ncol(L);
    t(mu + L%*%matrix(rnorm(m*n),m,n)) 
  }
  
  #Get predictors transformed into scale of knots of the non-linear gam functions.
  exclude.lab <- paste0('s(',ranef.lab,')')
  preds <- predict(fit, newdata = newdata, type = 'lpmatrix', exclude = c(exclude.lab), newdata.guaranteed = T)
  #grab columns that are random effects.
  if(is.na(ranef.lab)){
    lab.drop <- NULL
  }else{lab.drop <- grep(ranef.lab, colnames(preds))}

  
  #grab parameters and parameter variance covariance matrix, excluding random effects.
  if(length(lab.drop != 0)){
    pred.sub <- preds[,-lab.drop]
    par.sub <- coef(fit)[-lab.drop]
    par.cov.sub <- fit$Vp[-lab.drop, -lab.drop]
  }
  if(length(lab.drop) == 0){
       pred.sub <- preds
        par.sub <- coef(fit)
    par.cov.sub <- fit$Vp
  }
  
  #draw from parameter multivariate distribtuion to get correlated preditors.
  posterior <- rmvn(N.posterior, par.sub, par.cov.sub)
  
  #multiply predictors by posterior drawn parameters, summarize mean and sd (se.fit) of predictons.
  pred.out <- list()
  for(i in 1:N.posterior){
    pred.out[[i]] <- t(pred.sub %*% posterior[i,])
  }
  pred.out <- do.call(rbind, pred.out)
  pred.mu  <- colMeans(pred.out)
  pred.sd  <- apply(pred.out, 2, sd)
  
  #wrap and return output.
  output <- list(pred.mu, pred.sd)
  names(output) <- c('fit','se.fit')
  return(output)
  
} #end function.
