#Decompose the matrix representing the difference between the predicted and 
#smoothed state vectors into the components of the state variable, i.e. P_t*r_{t-1}
#on page 91 of Time Series Analysis by State Space Methods (Durbin and Koopman 2012).
#
#The logic is that the difference represents all the information contained in the 
#current observation as well as future observations, i.e. the contribution of new 
#information to the  estimation of the state vector. This is used to represent the 
#news content of successive vintages in the nowcasting framework.

decomp_smoother <- function(m) {
  # Return filtered and smoothed states
  kfs <- KFS(m$SSMvar, filtering="state", smoothing="state", simplify = FALSE)  
  # Return multivariate innovations (some quantities in KFS are not accurate for the multivariate model)
  mvi <- mvInnovations(kfs)

  # Calculate the inverse of a matrix dropping all rows and columns containing a zero on the diagonal
  inv <- function(X, index) {
    if(nrow(X[index,index]) != 0)
      X[index,index] <- solve(X[index,index])
    return(X) }
  
  # Dimensions of the observation matrix
  dimZ <- dim(m$SSMvar$Z)
  
  # Calculate F inverse and L matrices (pages 85 and 87 respectively)
  Finv <- array(NA,dim=c(dimZ[1],dimZ[1],dimZ[3]))
  L <- array(NA,dim=c(dimZ[2],dimZ[2],dimZ[3]))
  Finv[,,1] <- mvi$F[,,1]
  L[,,1] <- m$SSMvar$T[,,1]%*%(diag(dimZ[2]) - kfs$P[,,1]%*%t(m$SSMvar$Z[,,1])%*%Finv[,,1]%*%m$SSMvar$Z[,,1])
  for(i in 2:dimZ[3]) {
    Finv[,,i] <- inv(mvi$F[,,i], which(apply(m$SSMvar$Z[,,i],1,sum)!=0))
    L[,,i] <- m$SSMvar$T[,,1]%*%(diag(dimZ[2]) - kfs$P[,,i]%*%t(m$SSMvar$Z[,,i])%*%Finv[,,i]%*%m$SSMvar$Z[,,i])
  }
  
  # Calculate r_t and P_t*r_{t-1} - the rows of r are not "collapsed", hence it is an array and not a vector
  r_array <- array(NA, dim=c(dimZ[2],dimZ[1],dimZ[3]))
  smoother_update <- array(NA, dim=c(dimZ[2],dimZ[1],dimZ[3]))
  r_array[,,dimZ[3]] <- matrix(0, dimZ[2], dimZ[1])
  for(j in (dimZ[3]-1):1) {
    r_array[,,j] <- t(m$SSMvar$Z[,,j+1])%*%Finv[,,j+1]*matrix(mvi$v[j+1,],dimZ[2],dimZ[1],byrow=TRUE) + 
      t(L[,,j+1])%*%r_array[,,j+1]
    smoother_update[,,j+1] <- kfs$P[,,j+1]%*%r_array[,,j] 
  }
  smoother_update[,,1] <- kfs$P[,,1]%*%(t(m$SSMvar$Z[,,1])%*%Finv[,,1]*matrix(mvi$v[1,],dimZ[2],dimZ[1],byrow=TRUE) + 
                                          t(L[,,1])%*%r_array[,,1])
  return(smoother_update)
}