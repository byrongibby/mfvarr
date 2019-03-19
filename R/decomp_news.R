decomp_news <- function(m, update_months, variable="gdp") {
  # Return the decomposition of new information
  decomp <- decomp_smoother(m)
  
  # Index of variable of interest
  variable_index <- which(colnames(m$SSMvar$y)==variable)
  
  # Return news decomposition of variable of interest
  decomp_variable <- list()
  decomp_variable$monthly <- decomp[variable_index+1,,update_months]
  decomp_variable$quarterly <- apply(decomp_variable$monthly,1,mean)
  
  return(decomp_variable)
}
