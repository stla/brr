#' @name Beta2Dist 
#' @rdname Beta2Dist
#' @title Beta distribution of the second kind
#' 
#' @details The Beta distribution of the second kind with shape parameters 
#' \eqn{c>0} and \eqn{d>0} and scale parameter \eqn{k>0} is the distribution of 
#' \eqn{k*(U/(1-U))} where \eqn{U} is a random variable following the Beta distribution 
#' with shape parameters  \eqn{c} and \eqn{d}. 
#' \cr
#' It is also the distribution of 
#' 
#' @param x,q vector of quantiles 
#' @param p vector of probabilities
#' @param c,d non-negative shape parameters
#' @param scale non-negative scale parameter
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{FDist}}
#' 
#' @return \code{dbeta2} gives the density, \code{pbeta2} the distribution function, \code{qbeta2} the quantile function, and \code{rbeta2} generates random observations.
#' 
#' @note \code{Beta2Dist} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dbeta2(x, 3, 10, scale=2), from=0, to=3)
#' u <- rbeta(1e5, 3, 10)
#' lines(density(2*u/(1-u)), col="blue", lty="dashed")
#' 
NULL

#' @rdname Beta2Dist
#' @export 
dbeta2 <- function(x, c, d, scale, log=FALSE, ...){
  k <- d/c/scale
  switch(as.character(log),
         "FALSE" = k*df(k*x, df1=2*c, df2=2*d, log=FALSE, ...), 
         "TRUE" = log(k)+df(k*x, df1=2*c, df2=2*d, log=TRUE, ...)
  )
}
#' 
#' @rdname Beta2Dist 
#' @export 
pbeta2 <- function(q, c, d, scale, ...){
  pf(d/c/scale*q, df1=2*c, df2=2*d, ...)
}
#'
#' @rdname Beta2Dist 
#' @export 
rbeta2 <- function(nsims,c, d, scale){
  scale*c/d*rf(nsims, df1=2*c, df2=2*d)
}


#' @name Prior_mu 
#' @rdname Prior_mu
#' @title Prior distribution on the rate in the control group
#' 
#' @details The prior distribution on the rate \eqn{\mu} is the Gamma distribution 
#' with shape parameter \eqn{a} and rate parameter \eqn{b}
#' 
#' @param x,q vector of quantiles 
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{GammaDist}}
#' 
#' @return \code{dprior_mu} gives the density, \code{pprior_mu} the distribution function, \code{qprior_mu} the quantile function, and \code{rprior_mu} samples from the distribution.
#' 
#' @note \code{Prior_mu} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dprior_mu(x, 2, 2), from=0, to=3)
#' 
NULL

#' @rdname Prior_mu
#' @export 
dprior_mu<-function(x, a, b, ...){
  dgamma(x, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
pprior_mu<-function(q, a, b, ...){
  pgamma(q, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
qprior_mu<-function(p, a, b, ...){
  qgamma(p, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
rprior_mu<-function(n, a, b, ...){
  rgamma(nsims, a, b, ...)
}


#' @name Prior_phi 
#' @rdname Prior_phi
#' @title Prior distribution on the relative risk and the vaccine efficacy
#' 
#' @details The prior distribution on the relative risk \eqn{\phi} is the Beta2 distribution 
#' with shape parameters \eqn{c} and \eqn{d} and scale parameter \eqn{(T+b)/S}.
#' 
#' @param x,q vector of quantiles 
#' @param p vector of probabilities
#' @param b non-negative rate parameter
#' @param c,d non-negative shape parameters 
#' @param S,T sample sizes in control group and treated group
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{Beta2Dist}}
#' 
#' @return \code{dprior_phi} gives the density, \code{pprior_phi} the distribution function, \code{qprior_phi} the quantile function, and \code{rprior_phi} samples from the distribution.
#' 
#' @note \code{Prior_phi} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dprior_phi(x, 2, 2, 2, 10, 10), from=0, to=7)
#' 
NULL

#' @rdname Prior_phi
#' @export 
dprior_phi<-function(phi, b, c, d, S, T, ...){
  scale <- (T+b)/S
  dbeta2(phi,c,d,scale, ...)
}
#' @rdname Prior_phi
#' @export 
dprior_VE<-function(VE, b, c, d, S, T, ...){
  dprior_phi(1-VE, b, c, d, S, T, ...)
}
#'
#' @rdname Prior_phi
#' @export 
pprior_phi<-function(q, b, c, d, S, T, ...){
  scale <- (T+b)/S
  pbeta2(q,c,d,scale, ...)
}
#' @rdname Prior_phi
#' @export 
pprior_VE<-function(q, b, c, d, S, T, ...){
  1-pprior_phi(1-q, b, c, d, S, T, ...)
}

#'
#' @rdname Prior_phi
#' @export 
qprior_phi <- function(p, b, c, d, S, T, ...){
  scale <- (T+b)/S
  qbeta2(p, c, d, scale, ...)
}
#' @rdname Prior_phi
#' @export 
qprior_VE <- function(p, b, c, d, S, T, ...){
  1-qprior_phi(1-p, b, c, d, S, T, ...)
}
#'
#' @rdname Prior_phi
#' @export 
rprior_phi<-function(n, b, c, d, S, T){
  rbeta2(n, c, d, scale=(T+b)/S)
}
