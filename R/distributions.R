#' @name Beta2 
#' @rdname Beta2
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
#' @note \code{Beta2} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dbeta2(x, 3, 10, scale=2), from=0, to=3)
#' u <- rbeta(1e5, 3, 10)
#' lines(density(2*u/(1-u)), col="blue", lty="dashed")
#' 
NULL

#' @rdname Beta2
#' @export 
dbeta2 <- function(x, c, d, scale, log=FALSE, ...){
  k <- d/c/scale
  switch(as.character(log),
         "FALSE" = k*df(k*x, df1=2*c, df2=2*d, log=FALSE, ...), 
         "TRUE" = log(k)+df(k*x, df1=2*c, df2=2*d, log=TRUE, ...)
  )
}
#' 
#' @rdname Beta2 
#' @export 
pbeta2 <- function(q, c, d, scale, ...){
  pf(d/c/scale*q, df1=2*c, df2=2*d, ...)
}
#'
#' @rdname Beta2 
#' @export 
rbeta2 <- function(nsims,c, d, scale){
  scale*c/d*rf(nsims, df1=2*c, df2=2*d)
}