predict.MFVAR <- function(m, h, pctiles=c(0.05,0.50,0.95))
{
  # Model parameters
  K <- m$post$K; M <- 1; p <- m$post$p

  # Array for storing forecasts
  Yhat_mat <- matrix(NA,nrow(m$post$Y),ncol(m$post$Y)+h*K) 

  for(i in 1:nrow(m$post$Y)) {
    Y <- matrix(m$post$Y[i,], length(m$post$Y[i,])/K, K)
    Z <- matrix(m$post$Z[i,], length(m$post$Y[i,])/K, M+K*p)
    B <- matrix(m$post$B[i,], M+K*p, K)
    
    # Reformulate as VAR(1) - companion form 
    Acomp <- rbind(t(B[-(1:M),]), cbind(diag(K*(p-1)), matrix(0,K*(p-1),K)))
    Bcomp <- rbind(t(B[1:M,,drop=F]), matrix(0,K*(p-M),M))
    J     <- cbind(diag(K), matrix(0,K,K*(p-1)))

    # Forecast of Y and Z using estimated coefficients 
    Yhat <- rbind(m$par$Y[1:nrow(Y),], matrix(NA,h,K))
    Zhat <- rbind(m$par$Z[1:nrow(Y),], cbind(matrix(1,h,M),matrix(NA,h,K*p)))
    
    # Array for products of Acomp
    A_prd <- array(NA, dim=c(dim(Acomp),h+1))
    
    # Calculate products of Acomp
    A_prd[,,1] <- diag(K*p)
    for(j in 2:(h+1)) A_prd[,,j] <- A_prd[,,j-1] %*% Acomp

    # Forecast function
    for(j in 1:h) {
      dummy_sum <- matrix(0,K,1)
      for(k in 1:j) dummy_sum <- dummy_sum + 
        J%*%A_prd[,,k]%*%Bcomp%*%Zhat[nrow(Y)+(j+1)-k,1:M]
      Yhat[nrow(Y)+j,] <- J%*%A_prd[,,(j+1)]%*%Zhat[nrow(Y),-(1:M)] + dummy_sum
    }

    # Store forecasts
    Yhat_array[i,] <- as.vector(Yhat) 
  } 
  
  # ith percentile forecast
  out <- list()
  # Return the median as the posterior point estimate.
  out$median <- matrix(apply(post$Y,2,quantile,probs=0.50),nrow(m$obs),K)
  out$upper  <- matrix(apply(post$Y,2,quantile,probs=0.95),nrow(m$obs),K)
  out$lower  <- matrix(apply(post$Y,2,quantile,probs=0.05),nrow(m$obs),K)
  # Convert to time series.
  out$median <- ts(out$median,frequency=12,start=tsp(m$obs)[1])
  out$upper  <- ts(out$upper,frequency=12,start=tsp(m$obs)[1])
  out$lower  <- ts(out$lower,frequency=12,start=tsp(m$obs)[1])
  # Name the variables accordingly
  colnames(out$regressand$median) <- colnames(m$obs)
  colnames(out$regressand$upper) <- colnames(m$obs)
  colnames(out$regressand$lower) <- colnames(m$obs) 
  
  return(out)
}