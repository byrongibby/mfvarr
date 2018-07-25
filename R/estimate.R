estimate <-
function(m, reps, burn) {
  
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
                                           a1=m$filtmat$y1), 
                    H=matrix(0,K,K))
  ### List of function output
  post <- list()
  post$p <- p
  post$K <- K 
  post$Y <- matrix(NA, reps-burn, K*nrow(m$obs))
  post$Z <- matrix(NA, reps-burn, (K*p+M)*nrow(m$obs))
  post$A <- matrix(NA, reps-burn, K*(M+K*p))
  post$S <- matrix(NA, reps-burn, K*K)

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
      post$Z[i - burn,] <- as.vector(Z.draw[,])
      post$B[i - burn,] <- b.draw
      post$S[i - burn,] <- as.vector(S.draw)
    }
  }
  
  message("\nsampling complete.")
  close(pb)

  #------------------- OUTPUT -------------------#

  out <- list()
  out$post <- post
  out$regressand <- NULL 
  # Return the median as the posterior point estimate.
  out$regressand$median <- matrix(apply(post$Y,2,quantile,probs=0.50),nrow(m$obs),K)
  out$regressand$upper  <- matrix(apply(post$Y,2,quantile,probs=0.95),nrow(m$obs),K)
  out$regressand$lower  <- matrix(apply(post$Y,2,quantile,probs=0.05),nrow(m$obs),K)
  # Convert to time series.
  out$regressand$median <- ts(out$regressand$median,frequency=12,start=tsp(m$obs)[1])
  out$regressand$upper  <- ts(out$regressand$upper,frequency=12,start=tsp(m$obs)[1])
  out$regressand$lower  <- ts(out$regressand$lower,frequency=12,start=tsp(m$obs)[1])
  # Name the variables accordingly
  colnames(out$regressand$median) <- colnames(m$obs)
  colnames(out$regressand$upper) <- colnames(m$obs)
  colnames(out$regressand$lower) <- colnames(m$obs)
  return(out)
}
