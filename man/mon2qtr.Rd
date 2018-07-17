\name{mon2qtr}
\alias{mon2qtr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mon2qtr(monthly, aggregation = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monthly}{
%%     ~~Describe \code{monthly} here~~
}
  \item{aggregation}{
%%     ~~Describe \code{aggregation} here~~
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
function (monthly, aggregation = "mean") 
{
    month.start = round(tsp(monthly)[1]\%\%1 * 12) + 1
    year.start = round(tsp(monthly)[1])
    month.end = round(tsp(monthly)[2]\%\%1 * 12) + 1
    year.end = round(tsp(monthly)[2])
    quarter.start <- NULL
    if (any(month.start == c(2, 5, 8, 11))) {
        message("Dropping two months from start.")
        quarter.start <- month.start + 2
    }
    else if (any(month.start == c(3, 6, 9, 12))) {
        message("Dropping one month from start.")
        quarter.start <- month.start + 1
    }
    else {
        quarter.start <- month.start
    }
    quarter.end <- NULL
    if (any(month.end == c(1, 4, 7, 10))) {
        message("Dropping one month from end.")
        quarter.end <- month.end - 1
    }
    else if (any(month.end == c(2, 5, 8, 11))) {
        message("Dropping two months from end.")
        quarter.end <- month.end - 2
    }
    else {
        quarter.end <- month.end
    }
    monthly <- window(monthly, start = c(year.start, quarter.start), 
        end = c(year.end, quarter.end))
    agg.transform <- function(x, agg = "mean") {
        y <- switch(agg, mean = apply(x, 2, mean), sum = apply(x, 
            2, sum), last = apply(x, 2, function(z) z[length(z)]))
        return(y)
    }
    quarterly <- NULL
    for (i in 1:NCOL(monthly)) {
        monthly.reshape <- matrix(monthly[, i], 3, NROW(monthly)/3)
        quarterly <- cbind(quarterly, agg.transform(monthly.reshape))
    }
    quarterly <- ts(quarterly, freq = 4, start = c(year.start, 
        (quarter.start - 1)/3 + 1))
    colnames(quarterly) <- colnames(monthly)
    return(quarterly)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
