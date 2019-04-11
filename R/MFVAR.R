MFVAR <- function(monthly, quarterly, p=3, prior="default", mcmc="default") 
{
  # Restricted options for MFVAR
  type="const"
  dummy=NULL
 
  if(p < 3) p=3
  
  if(any(grepl("default",prior))) {
    prior = list("prior"="minnesota", 
                 "hyperparams"=list("tau"    = 0.1,
                                    "lambda" = 0.0, 
                                    "rho"    = 2.0, 
                                    "kappa"  = 1e5))
  }

  if(any(grepl("default",mcmc))) {
    mcmc = list("mcmc$reps" = 5000, 
                "mcmc$burn" = 4000, 
                "seed" = 1337)
  }

  #------------------- FORMAT DATA -------------------#

  # Breadth of monthly and quarterly data 
  Kq <- NCOL(quarterly)
  Km <- NCOL(monthly)
  # Start month of quarterly data (2 months prior to initial quarter)
  month_start <- c(floor(tsp(quarterly)[1]), round(tsp(quarterly)[1]%%1*12+1))
  # Create quarterly and monthly observation time series
  yq <- ts(matrix(rep(quarterly,each=3),NROW(quarterly)*3,Kq), 
           frequency = 12, 
           start = month_start)
  ym <- monthly
  # Data for estimating the initial VAR parameters
  y.init <- na.omit(cbind(ym, yq))
  colnames(y.init) <- c(colnames(monthly), colnames(quarterly))
  # Insert NA for unobserved data
  for(i in 1:nrow(yq)) if(i%%3 != 0) yq[i,] <- NA
  # Combine monthly and quarterly time series
  y <- cbind(ym, yq)
  colnames(y) <- c(colnames(monthly), colnames(quarterly))

  #------------------- 'TRAINING DATA' MEAN AND VARIANCE FOR PRIOR/DUMMY DATA -------------------#

  # number of variables
  K <- ncol(y.init)
  # Determine cut-off date for training data (10%)
  # NB! this also ensures that the data for estimation always starts on the
  # first month of a quarter (rest of code assumes this!)
  dates   <- time(y.init)
  cut.off <- c(floor(dates[ceiling(length(dates)*0.1)])-1,12)
  # Training data
  YY <- window(y.init, end=cut.off)
  # Sample means
  mu <- apply(YY,2,mean)
  # Variances of error terms from AR(1) regressions
  sigma <- rep(NA,K)
  for(i in 1:K) {
    yy  <- YY[2:nrow(YY),i]
    xx  <- cbind(1, YY[1:(nrow(YY)-1),i])
    b0 <- solve(crossprod(xx), crossprod(xx,yy))
    sigma[i] <- sqrt(crossprod(yy-xx%*%b0) / (NROW(yy)-2))
  }
  # Remove training data from dataset
  y      <- window(y, start=cut.off+c(0,1))
  y.init <- window(y.init, start=cut.off+c(0,1))
  # Create a prior for starting state (Kalman Filter) from training data
  y1 <- as.matrix(c(1,as.vector(t(YY[(nrow(YY)-p+1):nrow(YY),]))))

  #------------------- NOWCAST INFORMATION -------------------#
 
  # Length of data 
  N <- nrow(y)
   
  # The nowcast will be the first incomplete quarter 
  nowcast_quarter <- tsp(na.omit(quarterly))[2]+1/4
  # Nowcast quarter position vector
  nowcast_index <- which(round(nowcast_quarter*12)==round(time(y)*12))+0:2
  # Ensure the data includes the nowcast quarter 
  if(N < nowcast_index[3]) {
    # append NAs to end of data while preserving ts class
    y <- ts(rbind(y, matrix(NA,nowcast_index[3]-N,ncol(y))), freq=12, start=tsp(y)[1])
    # Update length of data 
    N <- nrow(y)
  }
  #------------------- CREATE OBSERVATION ARRAY -------------------#

  # Typical quarterly, monthly, and combined observation matrix
  L_qz <- matrix(rep(cbind(matrix(0,Kq,Km),diag(1/3,Kq)),3),Kq,Km*3+Kq*3)
  L_mz <- cbind(diag(Km),matrix(0,Km,Km*2+Kq*3))
  L_z  <- cbind(matrix(0,Kq+Km,1), rbind(L_mz,L_qz))
  # Allowance for VAR of order more than 3
  if(p > 3) L_z <- cbind(L_z,matrix(0,Kq+Km,(Kq+Km)*(p-3)))
  # Time varying observation array (allowing for unobserved data)
  ML_z <- array(L_z, dim=c(dim(L_z),N)) 
  # Create (column) indices for "unobserved observations"
  na.index.ym <- which(is.na(apply(y[,1:Km,drop=FALSE],1,prod)))
  na.index.yq <- which(is.na(apply(y[,(Km+1):(Km+Kq),drop=FALSE],1,prod)))
  # Replace all NAs with zero
  y[is.na(y)] <- 0
  # Replace transform on quarterly variables (L_qz) with zero where data unobserved
  for(i in na.index.yq) {
    if(i%%3 == 0)
      ML_z[Km+which(!(abs(y[i,(Km+1):(Km+Kq)])>0)),,i] <- 0
    else
      ML_z[(Km+1):(Km+Kq),,i] <- 0
  }
  # Replace transform on monthly variables (L_mz) with zero where data unobserved
  if(length(na.index.ym)>0) 
    for(i in na.index.ym) 
      ML_z[which(!(abs(y[i,1:Km])>0)),,i] <- 0

  #-------------- OLS ESTIMATES OF INITIAL STATE (VAR) MATRICES --------------#
   
  # Construct Z for Y = ZB + U
  z <- as.matrix(y.init) 
  Z <- matrix(NA, nrow=nrow(z)-p, ncol=K*p) 
  # Create new data matrix with p lags for all the variables
  for(i in 1:p) Z[, (K*(i-1) + 1):(K*i) ] <- z[((p+1) - i):(nrow(z) - i), ]
  # Constant (dd) and number of deterministic variables
  dd <- as.matrix(rep(1,nrow(z)))
  M  <- NCOL(dd) 
  # adjust data for lags
  Y <- z[(p+1):nrow(z),]
  Z <- cbind(dd[(p+1):nrow(z),],Z)
  # OLS estimate
  B  <- solve.qr(qr(Z),Y)
  # VAR (companion) matrices
  Acomp <- rbind(c(1,rep(0,K*p)),t(B), 
    cbind(rep(0,K*(p-1)),diag(K*(p-1)),matrix(0,K*(p-1),K)))
  Scomp <- matrix(0,K*p+M, K*p+M)
  Scomp[(M+1):(K+M),(M+1):(K+M)] <- crossprod(Y-Z%*%B)/(nrow(Y)-M-K*p)

  #------------------- DUMMY DATA FOR BVAR -------------------#
  
  # artificial data
  YD <- NULL
  ZD <- NULL
  # artificial data for Minnesota prior
  if(any(grepl("minnesota", prior$prior))) {
    m1 <- diag(sigma)
    m2 <- matrix(0,K*(p-1),K)
    m3 <- matrix(0,M,K)
    m4 <- matrix(0,K*p,M)
    m5 <- t(m3)
    m6 <- diag(M)
    m7 <- t(m4)
    m8 <- diag(1:p)
    m9 <- matrix(0,K,K*p)
    YD <- rbind(prior$hyperparams$lambda*m1*(1/prior$hyperparams$tau), m2, m1, m3)
    ZD <- cbind(rbind(m4, m5, m6/prior$hyperparams$kappa),
                rbind((m8^prior$hyperparams$rho)%x%m1*(1/prior$hyperparams$tau), m9, m7))
  }
  # artificial data for sum of coefficients prior
  if(any(grepl("sumcoefficients", prior$prior))) {
    m1 <- diag(mu)
    m2 <- matrix(0,K,M)
    m3 <- matrix(rep(1,p),1,p)
    YD <- rbind(YD, prior$hyperparams$lambda*m1/prior$hyperparams$gamma)
    ZD <- rbind(ZD, cbind(m2,m3%x%(prior$hyperparams$lambda*m1/prior$hyperparams$gamma)))
  }
  # artificial data for common stochastic trends prior
  if(any(grepl("stochastictrends", prior$prior))) {
    m1 <- c(1,rep(0,M-1))
    YD <- rbind(YD, prior$hyperparams$delta*mu)
    ZD <- rbind(ZD, prior$hyperparams$delta*c(m1, rep(mu,p)))
  }

  #------------------- SET UP FOR MCMC ROUTINE -------------------#
  
  set.seed(mcmc$seed)
  
  var_names <- c("constant", colnames(y))
  for(i in 2:p)
    var_names <- c(var_names, paste(colnames(y), ".l", i-1, sep=""))

  # Initial covariance matrix
  S.draw <- Scomp[(M+1):(K+M),(M+1):(K+M)]
  # State space model
  SSM <- SSModel(y ~ -1 + SSMcustom(Z=ML_z, 
                                    T=Acomp, 
                                    R=diag(M+K*p), 
                                    Q=Scomp, 
                                    a1=y1,
                                    state_names=var_names),
                 H=matrix(0,K,K))
  # Matrices to store draws of system matrices
  Y_draws <- matrix(NA, mcmc$reps-mcmc$burn, K*N)
  A_draws <- matrix(NA, mcmc$reps-mcmc$burn, K*(M+K*p))
  S_draws <- matrix(NA, mcmc$reps-mcmc$burn, K*K)
  
  #------------------- MCMC LOOP -------------------#
  
  message("Starting the Gibbs sampler...")
  pb <- txtProgressBar(style=3)
  i <- 0
  
  while(mcmc$reps > i) {
    setTxtProgressBar(pb, i/mcmc$reps)
    # Draw missing observations
    Z.draw <- simulateSSM(SSM, type="states")[,,1]
    # Append dummy data for VAR estimation
    Y <- rbind(Z.draw[,(M+1):(K+1)], YD)
    Z <- rbind(Z.draw, ZD)
    # VAR posterior mean and variance
    B. <- solve.qr(qr(Z),Y)
    b. <- as.vector(B.) 
    S. <- S.draw%x%solve(crossprod(Z))
    # Find nearest positive semi-definite matrix if S. doesn't factorise
    cholS. <- chol(S.)
    # Look for stable VAR parameters, abandon if not found after 100 draws
    for(try_i in 1:100) {
      # Draw VAR coefficients.
      b.draw <- b. + t(rnorm(K*(M+K*p))%*%cholS.)
      B.draw <- matrix(b.draw,M+K*p,K)
      # If stable VAR drawn...
      if(stable(B.draw,K,p,M,allow_ur=TRUE)) {
        # Draw VAR covariance coefficients
        S.draw <- riwish(nrow(Y), solve(crossprod(Y-Z%*%B.draw)))
        # Update VAR parameters in SSModel
        SSM["T"][(M+1):(K+M),,1] <- t(B.draw)
        SSM["Q"][(M+1):(K+M),(M+1):(K+M),1] <- S.draw
        # Advance sampler
        i = i + 1
        break
      }
    }
    # Save draws
    if(i > mcmc$burn) {
      Y_draws[i - mcmc$burn,] <- as.vector(Z.draw[,(M+1):(K+1)])
      A_draws[i - mcmc$burn,] <- b.draw
      S_draws[i - mcmc$burn,] <- as.vector(S.draw)
    }
  }
  
  message("\nsampling complete.")
  close(pb)
  
  #------------------- OUTPUT -------------------#
 
  # Probabilities for quantile function 
  probs <- c(0.05, 0.50, 0.95)
  # Return the regressand as a time series at the specified percentiles
  y_ <- list()
  for(i in probs) {
    y_[[paste("p",as.character(i*100), sep="")]] <- 
      ts(matrix(apply(Y_draws,2,quantile,probs=i),N,K),frequency=12,start=tsp(y)[1])
    colnames(y_[[paste("p",as.character(i*100), sep="")]]) <- colnames(y)
  }
  # Return the coefficient matrices at the specified percentiles
  A_ <- list()
  for(i in probs) {
    A_[[paste("p",as.character(i*100), sep="")]] <- 
      matrix(apply(A_draws,2,quantile,probs=i),ncol=K)
    colnames(A_[[paste("p",as.character(i*100), sep="")]]) <- colnames(y)
  }
  S_ <- list()
  for(i in probs) {
    S_[[paste("p",as.character(i*100), sep="")]] <- 
      matrix(apply(S_draws,2,quantile,probs=i),ncol=K)
    colnames(S_[[paste("p",as.character(i*100), sep="")]]) <- colnames(y)
    rownames(S_[[paste("p",as.character(i*100), sep="")]]) <- colnames(y)
  }
  
  # Update the state space model with the median posterior matrices
  SSM["T"][(M+1):(K+M),,1] <- t(A_$p50)
  SSM["Q"][(M+1):(K+M),(M+1):(K+M),1] <- S_$p50
 
  # Create a new state-space model where it is assumed that the 
  # quarterly variables are observed on a monthly basis as y_$p50
  SSM_fixed <- SSM
  SSM_fixed["Z"][1:K,(M+1):(M+K),1:(nowcast_index[1]-1)] <- diag(K)
  SSM_fixed["y"] <- y_$p50
  
  out <- list()
  out$nowcast_index <- nowcast_index
  out$pctiles <- list("y"=y_, "A"=A_, "S"=S_)
  out$SSM <- SSM
  out$SSM_fixed <- SSM_fixed 
  
  class(out) <- "MFVAR"
  
  return(out)
}
