stable <-
function(coef, nvar, nlag, ndet) {
  if(ndet > 0) { m1 <- array(t(coef[-(1:ndet),]), dim = c(nvar, nvar, nlag))
  } else { m1 <- array(t(coef), dim = c(nvar, nvar, nlag)) }
  m2 <- NULL
  for(i in 1:nlag) m2 <- cbind(m2, m1[,,i])
  m3 <- rbind(m2,cbind(diag(nvar*(nlag-1)), matrix(0,nvar*(nlag-1),nvar)) ) 
  out <- max(abs(eigen(m3)$values)) <= 1
  out
}
