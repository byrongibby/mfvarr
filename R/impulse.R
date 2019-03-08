#creates the sequence of matrices, phi (theta), that constitutes a part of the canonical
#representation laid out in Lutkepohl (2005, Ch2, pg. 15-18 & pg. 58)

impulse <- function(Amat, P, K, p, nsteps)
{
  if(p > 1)
  {
    A <- rbind(Amat, cbind(diag(K*(p-1)), matrix(0, K*(p-1), K) ) )
    J <- cbind(diag(K), matrix(0, K, (p-1)*K) ) 
    
    Aprod <- diag(K*p)
    phi   <- J%*%Aprod%*%t(J) 
    theta <- J%*%Aprod%*%t(J)%*%P
    
    for(h in 1:(nsteps - 1))
    {
      Aprod <- Aprod%*%A
      phi   <- cbind(phi,   J%*%Aprod%*%t(J)     ) 
      theta <- cbind(theta, J%*%Aprod%*%t(J)%*%P )
    }
    
    phi   <- array(phi,   dim = c(K, K, nsteps))
    theta <- array(theta, dim = c(K, K, nsteps))
    
    list("phi" = phi, "theta" = theta)
    
  } else
  {
    Aprod <- diag(K) 
    phi   <- Aprod
    theta <- P
    
    for(h in 1:(nsteps - 1))
    {
      Aprod <- Aprod%*%A
      phi   <- cbind(phi,   Aprod ) 
      theta <- cbind(theta, Aprod%*%P )
    }
    
    phi   <- array(phi,   dim = c(K, K, nsteps))
    theta <- array(theta, dim = c(K, K, nsteps))
    
    list("phi" = phi, "theta" = theta) 
  }
}