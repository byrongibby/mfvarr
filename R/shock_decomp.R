shock.decomp <- function(m) {
  
  #------------------- INPUT VALIDATION -------------------#
  K <- m$K; p <- m$p; M <- 1
  
  #------------------- COEFFICIENT MATRICES -------------------#
  
  # reduced form coefficients
  B    <- m$A$pctile_50
  A    <- t(B)[, -(1:M)]
  
  # regressand and regressor
  z <- as.matrix(m$monthly$pctile_50) 
  Z <- matrix(NA, nrow=nrow(z)-p, ncol=K*p) 
  # Create new data matrix with p lags for all the variables
  for(i in 1:p) Z[, (K*(i-1) + 1):(K*i) ] <- z[((p+1) - i):(nrow(z) - i), ]
  # Constant (dd) and number of deterministic variables
  dd <- as.matrix(rep(1,nrow(z)))
  M  <- NCOL(dd) 
  # adjust data for lags
  Y <- z[(p+1):nrow(z),]
  Z <- cbind(dd[(p+1):nrow(z),],Z)
  
  # reduced form shocks
  Uhat <- t(Y - Z%*%B)
  
  # sigma matrix
  Su <- tcrossprod(Uhat)/(ncol(Uhat) - M - K*p)
  
  #------------------- WOLD DECOMPOSITION -------------------#
  
  # Instantaneous impact matrix, cholesky decomposition -> recursive VAR
  theta0 <- t(chol(Su))
  
  # structural shocks
  What <- t(solve(theta0, Uhat))
  
  # structural deterministic variables
  B0star   <- solve(theta0, t(B)[, 1:M])
  B0star.x <- t(B0star%*%t(Z[, 1:M]))
  
  consts <- array(NA, dim = c(K, K, nrow(What)))
  shocks <- array(NA, dim = c(K, K, nrow(What)))
  
  # impulse-response function
  irf <- impulse(A, theta0, K, p, nrow(What))
  
  # Wold's decomposition
  # contributions to the various variables by the relevant 
  # historical structural shocks
  for(i in 1:K) {
    for(j in 1:K) {
      for(k in 1:nrow(What)) {
        consts[i,j,k] <- irf$theta[i,j,1:k] %*% B0star.x[k:1,j]
        shocks[i,j,k] <- irf$theta[i,j,1:k] %*% What[k:1,j]
      }
    }
  }
  
  #------------------- OUTPUT -------------------#
  sdt <- tsp(m$monthly$pctile_50)[1]+p/12
  frq <- tsp(m$monthly$pctile_50)[3]
  
  # formatted shocks
  out <- list()
  
  shock.names <- colnames(m$monthly$pctile_50)
  
  for(i in 1:K) {
    out[[colnames(m$monthly$pctile_50)[i] ]] <- 
      list("const"       = ts(cbind(apply(t(consts[i, , ]), 1, sum), t(shocks[i, , ])),
                              freq=frq, start=sdt, names=c("Constant",shock.names)),
           "omit"        = ts(t(shocks[i, , ]), 
                              freq=frq, start=sdt, names=shock.names), 
           "decomp"      = ts(t(shocks[i,,]) + t(consts[i,,]), 
                              freq=frq, start=sdt, names=shock.names),
           "series"      = window(m$monthly$pctile_50[, i], start=sdt),
           "series.omit" = window(m$monthly$pctile_50[, i], start=sdt) - apply(t(consts[i,,]), 1, sum))
  } 
 
  class(out) <- "MFVAR.decomp"
  return(out)
}