mean.BVAR <-
function(m) {
  cat("Unconditional posterior mean: \n\n")
  # model parameters
  K <- m$par$K; M <- m$par$M; p <- m$par$p
  
  # coefficient matrix
  b    <- apply(m$post$b, 2, mean)
  B    <- matrix(b,M+K*p,K)
  
  # create A(L) and A(1) matrix
  AL <- array(as.vector(t(B[-(1:M),])), dim=c(K,K,p))
  A1 <- diag(K)
  for(i in 1:p) A1 <- A1 - AL[,,i]
  
  # initial estimate of C0 - unconditional mean - (OLS)
  # WARNING: estimates may differ slightly from sample moments
  #          in the case of steady-state priors
  C0 <- t(solve(A1,t(B[(1:M),,drop=F])))
  
  # add constant mean to dummy period means
  if(M > 1) {
    for(i in 2:nrow(C0))
      C0[i,] <- C0[i,] + C0[1,]
  }
  
  # row and column labels
  varnames  <- colnames(m$regressand)
  dvarnames <- colnames(m$regressor[,1:m$par$M])
  colnames(C0) <- varnames
  rownames(C0) <- dvarnames
  
  round(C0,2)
}
