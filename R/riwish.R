riwish <-
function(v, S) {
  K <- NROW(S)
  Z <- apply(matrix(rnorm(K*v),K,v),2,function(x) as.vector(x%*%chol(S)))
  
  siginv <- matrix(0,K,K)
  for(i in 1:v) siginv <- siginv + tcrossprod(Z[,i])
  siginv <- solve(siginv)
  
  siginv
}
