summary.MFVAR <- function(m, probs=c(0.05, 0.50, 0.95))
{
  out <- list()
  
  out$K <- m$K
  out$p <- m$p
  
  ### Return the regressand as a time series at the specified percentiles
  out$monthly <- list()
  for(p in probs) {
    out$monthly[[paste("pctile_",as.character(p*100), sep="")]] <- 
      ts(matrix(apply(m$Y,2,quantile,probs=p),m$N,m$K),frequency=12,start=m$tsp[1])
    colnames(out$monthly[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  ### Provide the regressand at the quartlery frequency
  # Reshape matrix into an array
  post_Y_array <- array(t(m$Y),dim=c(m$N, m$K, m$draws)) 
  
  yy <- apply(post_Y_array, 3, 
              function(obs_array, N, K) {
                # Drop last quarter if incomplete, first quarter alwas complete (by assumption)
                monthly <- obs_array[1:(N-N%%3),]
                # Reshape and average over every 3 months to return quarterly series
                quarterly <- apply(array(monthly, dim=c(3, nrow(monthly)/3, K)), 2:3, mean)
              }, m$N, m$K)
  
  # Return the regressand as a time series at the specified percentiles
  out$quarterly <- list()
  for(p in probs) {
    out$quarterly[[paste("pctile_",as.character(p*100), sep="")]] <- 
      ts(matrix(apply(t(yy),2,quantile,probs=p),ncol=m$K),frequency=4,start=m$tsp[1])
    colnames(out$quarterly[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  # Return the coefficient matrix at the specified percentiles
  out$A <- list()
  for(p in probs) {
    out$A[[paste("pctile_",as.character(p*100), sep="")]] <- 
      matrix(apply(m$A,2,quantile,probs=p),ncol=m$K)
    colnames(out$A[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  # Print method to be added...

  
  class(out) <- "MFVAR.summary"
  return(out)
}