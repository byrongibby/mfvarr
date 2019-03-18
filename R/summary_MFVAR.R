summary.MFVAR <- function(m, probs=c(0.05, 0.50, 0.95))
{
  out <- list()
  
  out$K <- m$K
  out$p <- m$p
  M <- 1
    
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
                # Drop last quarter if incomplete, first quarter always complete (by assumption)
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
  
  # Return the coefficient matrices at the specified percentiles
  out$A <- list()
  for(p in probs) {
    out$A[[paste("pctile_",as.character(p*100), sep="")]] <- 
      matrix(apply(m$A,2,quantile,probs=p),ncol=m$K)
    colnames(out$A[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  out$S <- list()
  for(p in probs) {
    out$S[[paste("pctile_",as.character(p*100), sep="")]] <- 
      matrix(apply(m$S,2,quantile,probs=p),ncol=m$K)
    colnames(out$S[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
    rownames(out$S[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  # Reformulate estimated parameters in state-space form and return as a SSModel object
  Acomp <- rbind(c(1,rep(0,m$K*m$p)),t(out$A$pctile_50), 
                 cbind(rep(0,m$K*(m$p-1)),diag(m$K*(m$p-1)),matrix(0,m$K*(m$p-1),m$K)))
  Scomp <- matrix(0,m$K*m$p+M, m$K*m$p+M)
  Scomp[(M+1):(m$K+M),(M+1):(m$K+M)] <- out$S$pctile_50
  
  out$SSMvar <- SSModel(m$obs ~ -1 + SSMcustom(Z=m$ML_z, 
                                               T=Acomp, 
                                               R=diag(M+m$K*m$p), 
                                               Q=Scomp, 
                                               a1=m$y1), 
                        H=matrix(0,m$K,m$K))
  
  # Print method to be added...

  
  class(out) <- "MFVAR.summary"
  return(out)
}