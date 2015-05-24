#' @name Beta2Dist 
#' @rdname Beta2Dist
#' @title Beta distribution of the second kind
#' @description Density, distribution function, quantile function and random 
#' generation for the Beta distribution of the second kind with shape parameters 
#' \code{c} and \code{d} and scale parameter \code{scale}. 
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
#'
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


#' @name GammaInverseBetaDist 
#' @rdname GammaInverseBetaDist
#' @title Gamma-Inverse Beta distribution
#' @description Density and random  generation for the Gamma-Inverse Beta distribution 
#' with shape parameters \code{a}, \code{c}, \code{d} and scale parameter \code{rho}. 
#' @details This is the mixture distribution obtained by sampling a value \eqn{b} from a Beta distribution with parameters \eqn{\alpha}, \eqn{\beta} 
#' and then sampling a Gamma distribution with shape \eqn{a} and rate \eqn{\rho/b}.
#' 
#' @param x vector of quantiles
#' @param a non-negative shape parameter of the Gamma distribution
#' @param alpha,beta non-negative shape parameters of the mixing Beta distribution
#' @param n number of observations to be simulated
#' 
#' @return \code{dGIB} gives the density, and \code{rGIB} samples from the distribution.
#' 
#' @note \code{GammaInverseBetaDist } is a generic name for the functions documented. 
#' 
#' @importFrom gsl lnpoch lngamma hyperg_U
#' @examples
#' curve(dGIB(x,3,4,2,2.5), from=0, to=3)
#' sims <- rgamma(100000, 3, 2.5/rbeta(100000,4,2))
#' lines(density(sims, from=0), col="red")
#' lines(density(rGIB(100000, 3, 4, 2, 2.5), from=0), col="green")
#' 
NULL
#'
#' @rdname GammaInverseBetaDist
#' @export
dGIB <- function(x,a,alpha,beta,rho){
  exp(lnpoch(alpha,beta)-lngamma(a))*rho^a*x^(a-1)*exp(-rho*x)*hyperg_U(beta,a-alpha+1,rho*x)
}
#'
#' @rdname GammaInverseBetaDist
#' @export
rGIB <- function(n,a,alpha,beta,rho){
  rbeta(n,alpha,beta)*rgamma(n, a, rho)
}
