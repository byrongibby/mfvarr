summary.MFVAR <- function(m, probs=c(0.05, 0.50, 0.95))
{
  out <- list()
  
  # Return the regressand as a time series at the specified percentiles
  out$monthly <- list()
  for(p in probs) {
    out$monthly[[paste("pctile_",as.character(p*100), sep="")]] <- 
      ts(matrix(apply(m$Y,2,quantile,probs=p),m$N,m$K),frequency=12,start=m$tsp[1])
    colnames(out$monthly[[paste("pctile_",as.character(p*100), sep="")]]) <- m$names
  }
  
  # Provide the regressand at the quartlery frequency
  out$quarterly <- lapply(out$monthly, mon2qtr)
  
  # Parameter estimates and print method to be added...
  
  class(out) <- "MFVAR.summary"
  return(out)
}