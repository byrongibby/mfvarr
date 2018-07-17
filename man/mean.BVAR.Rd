\name{mean.BVAR}
\alias{mean.BVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mean.BVAR(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (m) 
{
    cat("Unconditional posterior mean: \n\n")
    K <- m$par$K
    M <- m$par$M
    p <- m$par$p
    b <- apply(m$post$b, 2, mean)
    B <- matrix(b, M + K * p, K)
    AL <- array(as.vector(t(B[-(1:M), ])), dim = c(K, K, p))
    A1 <- diag(K)
    for (i in 1:p) A1 <- A1 - AL[, , i]
    C0 <- t(solve(A1, t(B[(1:M), , drop = F])))
    if (M > 1) {
        for (i in 2:nrow(C0)) C0[i, ] <- C0[i, ] + C0[1, ]
    }
    varnames <- colnames(m$regressand)
    dvarnames <- colnames(m$regressor[, 1:m$par$M])
    colnames(C0) <- varnames
    rownames(C0) <- dvarnames
    round(C0, 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
