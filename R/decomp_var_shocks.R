decomp_shocks <- function(m, recursive=FALSE, theta0=NULL) {
  
  #------------------- COEFFICIENT MATRICES -------------------#
  
  # reduced form coefficients
  B    <- m$VAR$A
  A    <- t(B)[, -(1:m$VAR$M)]
  
  # regressand and regressor
  z <- as.matrix(m$VAR$y) 
  Z <- matrix(NA, nrow=nrow(z)-m$VAR$p, ncol=m$VAR$K*m$VAR$p) 
  # Create new data matrix with p lags for all the variables
  for(i in 1:m$VAR$p) Z[, (m$VAR$K*(i-1) + 1):(m$VAR$K*i) ] <- z[((m$VAR$p+1) - i):(nrow(z) - i), ]
  # Constant (dd) and number of deterministic variables
  dd <- as.matrix(rep(1,nrow(z)))
  # adjust data for lags
  Y <- z[(m$VAR$p+1):nrow(z),]
  Z <- cbind(dd[(m$VAR$p+1):nrow(z),],Z)
  
  # reduced form shocks
  Uhat <- Y - Z%*%B
  
  # sigma matrix
  Su <- tcrossprod(t(Uhat))/(nrow(Uhat) - m$VAR$M - m$VAR$K*m$VAR$p)
  
  #------------------- WOLD DECOMPOSITION -------------------#
  structural <- recursive || !is.null(theta0)
  
  # Instantaneous impact matrix, cholesky decomposition -> recursive VAR
  if(recursive || is.null(theta0)) theta0 <- t(chol(Su))
  
  consts <- array(NA, dim = c(m$VAR$K, m$VAR$K, nrow(Uhat)))
  shocks <- array(NA, dim = c(m$VAR$K, m$VAR$K, nrow(Uhat)))
  
  # impulse-response function
  irf <- impulse(A, theta0, m$VAR$K, m$VAR$p, nrow(Uhat))
  
  if(structural) {
    # Structural shocks
    What <- t(solve(theta0, t(Uhat)))
    
    # Structural deterministic variables
    B0star   <- solve(theta0, t(B)[, 1:m$VAR$M])
    B0star.x <- t(B0star%*%t(Z[, 1:m$VAR$M]))
    
    # Wold's decomposition
    # contributions to the various variables by the relevant 
    # historical structural shocks
    for(i in 1:m$VAR$K) {
      for(j in 1:m$VAR$K) {
        for(k in 1:nrow(What)) {
          consts[i,j,k] <- irf$theta[i,j,1:k] %*% B0star.x[k:1,j]
          shocks[i,j,k] <- irf$theta[i,j,1:k] %*% What[k:1,j]
        }
      }
    }
    
  } else {
    # Reduced form deterministic variables
    B0   <- t(B)[, 1:m$VAR$M]
    B0.x <- t(B0%*%t(Z[, 1:m$VAR$M]))
    
    # Wold's decomposition
    # contributions to the various variables by the relevant 
    # historical structural shocks
    for(i in 1:m$VAR$K) {
      for(j in 1:m$VAR$K) {
        for(k in 1:nrow(Uhat)) {
          consts[i,j,k] <- irf$phi[i,j,1:k] %*% B0.x[k:1,j]
          shocks[i,j,k] <- irf$phi[i,j,1:k] %*% Uhat[k:1,j]
        }
      }
    }
  }
  
  #------------------- OUTPUT -------------------#
  
  sdt <- tsp(m$VAR$y)[1]+m$VAR$p/12
  frq <- tsp(m$VAR$y)[3]
  
  # formatted shocks
  out <- list()
  
  shock.names <- colnames(m$VAR$y)
  
  for(i in 1:m$VAR$K) {
    out[[colnames(m$VAR$y)[i] ]] <- 
      list("const"       = ts(cbind(apply(t(consts[i, , ]), 1, sum), t(shocks[i, , ])),
                              freq=frq, start=sdt, names=c("Constant",shock.names)),
           "omit"        = ts(t(shocks[i, , ]), 
                              freq=frq, start=sdt, names=shock.names), 
           "decomp"      = ts(t(shocks[i,,]) + t(consts[i,,]), 
                              freq=frq, start=sdt, names=shock.names),
           "series"      = window(m$VAR$y[, i], start=sdt),
           "series.omit" = window(m$VAR$y[, i], start=sdt) - apply(t(consts[i,,]), 1, sum))
  } 
 
  class(out) <- "VAR.shocks"
  return(out)
}