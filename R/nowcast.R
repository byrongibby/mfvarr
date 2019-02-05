nowcast <- function(m)
{
  out <- list()
  # Nowcast-quarter draws
  Y_array      <- array(t(m$Y), dim=c(m$N,m$K,m$draws))
  last_quarter <- apply(Y_array[m$nowcast$index,,], c(2,3), mean)
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