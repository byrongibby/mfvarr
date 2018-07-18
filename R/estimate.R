estimate <-
function(m, reps, burn, save.draws=FALSE) {
  
  #------------------- SET UP -------------------#

  # VAR model parameters 
  M <- 1
  K <- ncol(m$obs)
  p <- m$lag
  YD <- m$dummy$YD
  ZD <- m$dummy$ZD
  S.draw <- m$filtmat$Scomp[(M+1):(K+M),(M+1):(K+M)]
  # State space model
  SSMvar <- SSModel(m$obs ~ -1 + SSMcustom(Z=m$filtmat$ML_z, 
                                           T=m$filtmat$Acomp, 
                                           R=diag(M+K*p), 
                                           Q=m$filtmat$Scomp, 
                                           a1=m$filtmat$y1,
                                           P1=m$filtmat$P1), 
                    H=matrix(0,K,K))
  ### List of function output
  post <- list()
  post$Y    <- matrix(NA, reps-burn, K*nrow(m$obs))
  post$A    <- matrix(NA, reps-burn, K*(M+K*p))
  post$S    <- matrix(NA, reps-burn, K*K)

  #------------------- MCMC LOOP -------------------#
  
  message("Starting the Gibbs sampler...")
  pb <- txtProgressBar(style=3)
  i <- 0
  
  while(reps > i) {
    setTxtProgressBar(pb, i/reps)
    # Draw missing observations.
    Z.draw <- simulateSSM(SSMvar, type="states")[,,1]
    # Append dummy data for VAR estimation.
    Y <- rbind(Z.draw[,(M+1):(K+1)], YD)
    Z <- rbind(Z.draw, ZD)
    # VAR posterior mean and variance.
    B. <- solve.qr(qr(Z),Y)
    b. <- as.vector(B.) 
    S. <- S.draw%x%solve(crossprod(Z))
    # Find nearest positive semi-definite matrix if S. doesn't factorise.
    cholS. <- chol(S.)
    # Look for stable VAR parameters, abandon if not found after 100 draws.
    for(try_i in 1:100) {
      # Draw VAR coefficients.
      b.draw <- b. + t(rnorm(K*(M+K*p))%*%cholS.)
      B.draw <- matrix(b.draw,M+K*p,K)
      # If stable VAR drawn...
      if(stable(B.draw,K,p,M,allow_ur=TRUE)) {
        # Draw VAR covariance coeffcients.
        S.draw <- riwish(nrow(Y), solve(crossprod(Y-Z%*%B.draw)))
        # Update VAR parameters in SSModel.
        SSMvar["T"][(M+1):(K+M),,1] <- t(B.draw)
        SSMvar["Q"][(M+1):(K+M),(M+1):(K+M),1] <- S.draw
        # Advance sampler.
        i = i + 1
        break
      }
    }
    # Save draws.
    if(i > burn) {
      post$Y[i - burn,] <- as.vector(Z.draw[,(M+1):(K+1)])
      post$A[i - burn,] <- b.draw
      post$S[i - burn,] <- as.vector(S.draw)
    }
  }
  
  message("\nsampling complete.")
  close(pb)

  #------------------- OUTPUT -------------------#

  out <- list()
  if(save.draws) out$post <- post
  # Return the median as the posterior point estimate.
  out$regressand <- matrix(apply(post$Y,2,median),nrow(m$obs),K)
  out$A          <- matrix(apply(post$A,2,median),K,M+K*p)
  out$S          <- matrix(apply(post$S,2,median),K,K)
  # Convert to time series.
  out$regressand <- ts(out$regressand,frequency=12,start=tsp(m$obs)[1])
  colnames(out$regressand) <- colnames(m$obs)
  return(out)
}
