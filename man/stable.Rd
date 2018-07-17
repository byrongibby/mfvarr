\name{stable}
\alias{stable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stable(coef, nvar, nlag, ndet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{
%%     ~~Describe \code{coef} here~~
}
  \item{nvar}{
%%     ~~Describe \code{nvar} here~~
}
  \item{nlag}{
%%     ~~Describe \code{nlag} here~~
}
  \item{ndet}{
%%     ~~Describe \code{ndet} here~~
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
function (coef, nvar, nlag, ndet) 
{
    if (ndet > 0) {
        m1 <- array(t(coef[-(1:ndet), ]), dim = c(nvar, nvar, 
            nlag))
    }
    else {
        m1 <- array(t(coef), dim = c(nvar, nvar, nlag))
    }
    m2 <- NULL
    for (i in 1:nlag) m2 <- cbind(m2, m1[, , i])
    m3 <- rbind(m2, cbind(diag(nvar * (nlag - 1)), matrix(0, 
        nvar * (nlag - 1), nvar)))
    out <- max(abs(eigen(m3)$values)) <= 1
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
