nowcast <- function(m)
{
  # Nowcast-quarter draws
  Y_array      <- array(t(m$Y), dim=c(nrow(m$obs),K,reps-burn))
  last_quarter <- apply(Y_array[(nrow(m$Y)-2):nrow(m$Y),,], c(2,3), mean)
  out$draws    <- split(last_quarter, slice.index(last_quarter, 1))
  # Name the densities after their corresponding variables.
  names(out$draws) <- m$names
  # Empirical cumulative densities
  out$ecdf <- lapply(out$draws, ecdf)
  # Percentiles
  out$pctiles <- lapply(out$draws, function(x) { 
    return(function(probs) {quantile(x, probs)}) })
  class(out) <- "nowcast"
  return(out)
}