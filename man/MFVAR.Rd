\name{MFVAR}
\alias{MFVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MFVAR(monthly, quarterly, p = 3, prior = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monthly}{
%%     ~~Describe \code{monthly} here~~
}
  \item{quarterly}{
%%     ~~Describe \code{quarterly} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
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
function (monthly, quarterly, p = 3, prior = "default") 
{
    type = "const"
    dummy = NULL
    nq <- NCOL(quarterly)
    nm <- NCOL(monthly)
    month.start <- c(floor(tsp(quarterly)[1]), tsp(quarterly)[1]\%\%1 * 
        12 + 1)
    xq <- ts(matrix(rep(quarterly, each = 3), NROW(quarterly) * 
        3, nq), frequency = 12, start = month.start)
    xm <- monthly
    y.init <- na.omit(cbind(xm, xq))
    colnames(y.init) <- c(colnames(monthly), colnames(quarterly))
    for (i in 1:nrow(xq)) if (i\%\%3 != 0) 
        xq[i, ] <- NA
    x <- cbind(xm, xq)
    colnames(x) <- c(colnames(monthly), colnames(quarterly))
    K <- ncol(y.init)
    dates <- time(y.init)
    cut.off <- c(floor(dates[ceiling(length(dates) * 0.1)]) - 
        1, 12)
    YY <- window(y.init, end = cut.off)
    mu <- apply(YY, 2, mean)
    sigma <- rep(NA, K)
    for (i in 1:K) {
        yy <- YY[2:nrow(YY), i]
        xx <- cbind(1, YY[1:(nrow(YY) - 1), i])
        b0 <- solve(crossprod(xx), crossprod(xx, yy))
        sigma[i] <- sqrt(crossprod(yy - xx \%*\% b0)/(NROW(yy) - 
            2))
    }
    x <- window(x, start = cut.off + c(0, 1))
    y.init <- window(y.init, start = cut.off + c(0, 1))
    y1 <- as.matrix(c(1, as.vector(t(YY[(nrow(YY) - p + 1):nrow(YY), 
        ]))))
    L_qz <- matrix(rep(cbind(matrix(0, nq, nm), diag(1/3, nq)), 
        3), nq, nm * 3 + nq * 3)
    L_mz <- cbind(diag(nm), matrix(0, nm, nm * 2 + nq * 3))
    L_z <- cbind(matrix(0, nq + nm, 1), rbind(L_mz, L_qz))
    if (p > 3) 
        L_z <- cbind(L_z, matrix(0, nq + nm, (nq + nm) * (p - 
            3)))
    ML_z <- array(L_z, dim = c(dim(L_z), nrow(x)))
    na.index.xm <- which(is.na(apply(x[, 1:nm, drop = FALSE], 
        1, prod)))
    na.index.xq <- which(is.na(apply(x[, (nm + 1):(nm + nq), 
        drop = FALSE], 1, prod)))
    x[is.na(x)] <- 0
    for (i in na.index.xq) {
        if (i\%\%3 == 0) 
            ML_z[nm + which(!(abs(x[i, (nm + 1):(nm + nq)]) > 
                0)), , i] <- 0
        else ML_z[(nm + 1):(nm + nq), , i] <- 0
    }
    if (length(na.index.xm) > 0) 
        for (i in na.index.xm) ML_z[1:nm, (1 + 1):(nm + 1), i] <- diag(abs(x[i, 
            1:nm]) > 0) * 1
    z <- as.matrix(y.init)
    Z <- matrix(NA, nrow = nrow(z) - p, ncol = K * p)
    for (i in 1:p) Z[, (K * (i - 1) + 1):(K * i)] <- z[((p + 
        1) - i):(nrow(z) - i), ]
    dd <- as.matrix(rep(1, nrow(z)))
    M <- NCOL(dd)
    Y <- z[(p + 1):nrow(z), ]
    Z <- cbind(dd[(p + 1):nrow(z), ], Z)
    B <- solve(crossprod(Z), crossprod(Z, Y))
    Acomp <- rbind(c(1, rep(0, K * p)), t(B), cbind(rep(0, K * 
        (p - 1)), diag(K * (p - 1)), matrix(0, K * (p - 1), K)))
    S <- crossprod(Y - Z \%*\% B)/(nrow(Y) - M - K * p)
    Scomp <- matrix(0, K * p + M, K * p + M)
    Scomp[(M + 1):(K + M), (M + 1):(K + M)] <- S
    YD <- NULL
    ZD <- NULL
    if (any(prior$prior == "mn")) {
        m1 <- diag(sigma)
        m2 <- matrix(0, K * (p - 1), K)
        m3 <- if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else matrix(0, M, K)
        m4 <- if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else matrix(0, K * p, M)
        m5 <- if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else t(m3)
        m6 <- if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else diag(M)
        m7 <- if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else t(m4)
        m8 <- diag(1:p)
        m9 <- matrix(0, K, K * p)
        YD <- rbind(prior$hp$lambda * m1 * (1/prior$hp$tau), 
            m2, m1, m3)
        ZD <- cbind(if (M == 0 || any(prior$prior == "ss")) 
            NULL
        else rbind(m4, m5, m6/prior$hp$kappa), rbind((m8^prior$hp$rho) \%x\% 
            m1 * (1/prior$hp$tau), m9, m7))
    }
    if (any(prior$prior == "sc")) {
        m1 <- diag(mu)
        m2 <- if (M == 0) 
            NULL
        else matrix(0, K, M)
        m3 <- matrix(rep(1, p), 1, p)
        YD <- rbind(YD, prior$hp$lambda * m1/prior$hp$gamma)
        ZD <- rbind(ZD, cbind(m2, m3 \%x\% (prior$hp$lambda * m1/prior$hp$gamma)))
    }
    if (any(prior$prior == "st")) {
        YD <- rbind(YD, prior$hp$delta * mu)
        m1 <- if (type == "const") 
            c(1, rep(0, M - 1))
        else rep(0, M)
        ZD <- rbind(ZD, prior$hp$delta * c(m1, rep(mu, p)))
    }
    out <- list(obs = x, lag = p, filtmat = list(ML_z = ML_z, 
        Acomp = Acomp, Scomp = Scomp, y1 = y1), dummy = list(YD = YD, 
        ZD = ZD))
    class(out) <- "MFVAR model"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
