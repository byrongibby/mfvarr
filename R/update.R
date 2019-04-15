update <- function(monthly, quarterly, m, nowcast=FALSE)
{
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
  # Data for estimating the starting state of the Kalman Filter 
  y.init <- na.omit(cbind(ym, yq))
  colnames(y.init) <- c(colnames(monthly), colnames(quarterly))
  # Insert NA for unobserved data
  for(i in 1:nrow(yq)) if(i%%3 != 0) yq[i,] <- NA
  # Combine monthly and quarterly time series
  y <- cbind(ym, yq)
  colnames(y) <- c(colnames(monthly), colnames(quarterly))

  #------------------- FURTHER DATA FORMATTING AND INITIALISATION FOR FILTER -------------------#

  # number of variables
  K <- ncol(y.init)
  # Determine cut-off date for training data (10%)
  # NB! this also ensures that the data for estimation always starts on the
  # first month of a quarter (rest of code assumes this!)
  dates   <- time(y.init)
  cut.off <- c(floor(dates[ceiling(length(dates)*0.1)])-1,12)
  # Training data
  YY <- window(y.init, end=cut.off)
  # Remove training data from dataset
  y      <- window(y, start=cut.off+c(0,1))
  y.init <- window(y.init, start=cut.off+c(0,1))
  # Create a prior for starting state (Kalman Filter) from training data
  y1 <- as.matrix(c(1,as.vector(t(YY[(nrow(YY)-p+1):nrow(YY),]))))

  #------------------- NOWCAST INFORMATION AND FORMATTING -------------------#
 
  if(nowcast) {
    # The nowcast will be the first incomplete quarter 
    nowcast_quarter <- tsp(na.omit(quarterly))[2]+1/4
    # Nowcast quarter position vector
    nowcast_index <- round((nowcast_quarter-tsp(y)[1])*12)+1:3
    # Ensure the data includes the full nowcast quarter 
    if(nrow(y) < nowcast_index[3]) {
      # append NAs to end of data while preserving ts class
      y <- ts(rbind(y, matrix(NA,nowcast_index[3]-nrow(y),ncol(y))), freq=12, start=tsp(y)[1])
    }
  }
  
  #------------------- CREATE OBSERVATION ARRAY -------------------#

  # Length of data 
  N <- nrow(y)

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


  #------------------- UPDATE DATA AND OBSERVATION MATRICES -------------------#

  # Update the state-space model with the new data and observation equation
  m$SSM["y1"] <- y1
  m$SSM["y"] <- y
  m$SSM["Z"] <- ML_z

  # Don't update SSM_fixed if the nowcast quarter has possibly changed
  if(!nowcast) {
    # Update the state-space model with fixed monthly observations
    m$SSM_fixed["y"][m$nowcast_index,] <- y[m$nowcast_index,]
    m$SSM_fixed["Z"][,,m$nowcast_index] <- ML_z[,,m$nowcast_index]
  } else {
    warning("SSM_fixed not updated becasuse nowcast quarter has changed.")
  }

  out <- m 
  
  return(out)
}
