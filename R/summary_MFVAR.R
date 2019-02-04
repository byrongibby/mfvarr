summary.MFVAR <- function(m, probs=c(0.05, 0.50, 0.95))
{
  out <- list()
  
  # Return the regressand as a time series at the specified percentiles
  out$monthly <- list()
  for(p in probs) {
    out$monthly[[paste("_",as.character(p*100), "pct", sep="")]] <- 
      ts(matrix(apply(post$Y,2,quantile,probs=p),nrow(m$Y),K),frequency=12,start=m$tsp[1])
    colnames(out$monthly[[paste("_",as.character(p*100), "pct", sep="")]]) <- m$names
  }
  
  # Provide the regressand at the quartlery frequency
  out$quarterly <- lapply(out$monthly, mon2qtr)
  
  # Parameter estimates and print method to be added...
  
  class(out) <- "MFVAR.summary"
  return(out)
}