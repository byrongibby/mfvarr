decomp_news <- function(m, variable="gdp") {
  # Return the decomposition of new information
  decomp <- decomp_smoother(m)
  
  # Index of variable of interest
  variable_index <- which(colnames(m$SSM$y)==variable)
  
  a <- decomp$a[m$nowcast_index, variable_index+1]
  alphahat <- decomp$alphahat[m$nowcast_index, variable_index+1]
  smoother_update <- decomp$smoother_update[variable_index+1,,m$nowcast_index]
  
  scale_factor <- m$pctiles[[variable_index]](0.5)/mean(alphahat)
  
  # Return news decomposition of variable of interest
  decomp_variable <- list()
  decomp_variable$a <- mean(a)*scale_factor
  decomp_variable$alphahat <- mean(alphahat)*scale_factor
  decomp_variable$smoother_update <- apply(smoother_update,1,mean)*scale_factor
  
  return(decomp_variable)
}
