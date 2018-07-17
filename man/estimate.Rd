\name{estimate}
\alias{estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate(m, reps, burn, save.draws = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{burn}{
%%     ~~Describe \code{burn} here~~
}
  \item{save.draws}{
%%     ~~Describe \code{save.draws} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m, reps, burn, save.draws = FALSE) 
{
    M <- 1
    K <- ncol(m$obs)
    p <- m$lag
    YD <- m$dummy$YD
    ZD <- m$dummy$ZD
    S.draw <- m$filtmat$Scomp[(M + 1):(K + M), (M + 1):(K + M)]
    SSMvar <- SSModel(m$obs ~ -1 + SSMcustom(Z = m$filtmat$ML_z, 
        T = m$filtmat$Acomp, R = diag(M + K * p), Q = m$filtmat$Scomp, 
        a1 = m$filtmat$y1), H = matrix(0, K, K))
    post <- list()
    post$Y <- matrix(NA, reps - burn, K * nrow(m$obs))
    post$A <- matrix(NA, reps - burn, K * (M + K * p))
    post$S <- matrix(NA, reps - burn, K * K)
    message("Starting the Gibbs sampler...")
    pb <- txtProgressBar(style = 3)
    i <- 0
    while (reps > i) {
        setTxtProgressBar(pb, i/reps)
        Z.draw <- simulateSSM(SSMvar, type = "states")[, , 1]
        Y <- rbind(Z.draw[, (M + 1):(K + 1)], YD)
        Z <- rbind(Z.draw, ZD)
        B. <- solve(crossprod(Z), crossprod(Z, Y))
        b. <- as.vector(B.)
        S. <- S.draw \%x\% solve(crossprod(Z))
        try(cholS. <- chol(S.), cholS. <- chol(nearestSPD(S.)))
        for (try_i in 1:100) {
            b.draw <- b. + t(rnorm(K * (M + K * p)) \%*\% cholS.)
            B.draw <- matrix(b.draw, M + K * p, K)
            if (stable(B.draw, K, p, M)) {
                S.draw <- riwish(nrow(Y), solve(crossprod(Y - 
                  Z \%*\% B.draw)))
                SSMvar["T"][(M + 1):(K + M), , 1] <- t(B.draw)
                SSMvar["Q"][(M + 1):(K + M), (M + 1):(K + M), 
                  1] <- S.draw
                i = i + 1
                break
            }
        }
        if (i > burn) {
            post$Y[i - burn, ] <- as.vector(Z.draw[, (M + 1):(K + 
                1)])
            post$A[i - burn, ] <- b.draw
            post$S[i - burn, ] <- as.vector(S.draw)
        }
    }
    message("\nsampling complete.")
    close(pb)
    out <- list()
    if (save.draws) 
        out$post <- post
    out$regressand <- matrix(apply(post$Y, 2, median), nrow(m$obs), 
        K)
    out$A <- matrix(apply(post$A, 2, median), K, M + K * p)
    out$S <- matrix(apply(post$S, 2, median), K, K)
    out$regressand <- ts(out$regressand, frequency = 12, start = tsp(m$obs)[1])
    colnames(out$regressand) <- colnames(m$obs)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
