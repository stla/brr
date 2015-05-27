#' @name Posterior_lambda 
#' @rdname Posterior_lambda
#' @title Posterior distribution on the incidence rate in the treated group
#' @description Density and random generation for the posterior distribution on 
#' the rate in the treated group. The distribution function and the quantile function 
#' are not available.
#' @details The pdf of the posterior distribution of the incidence rate \eqn{\lambda} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' 
#' 
#' @param lambda vector of quantiles 
#' @param a non-negative shape parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S sample size in treated group
#' @param x,y counts in the treated group and control group 
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{dGIB}}
#' 
#' @return \code{dpost_lambda} gives the density, and \code{rpost_lambda} samples from the distribution.
#' 
#' @note \code{Posterior_lambda} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dpost_lambda(x, 2, 2, 2, 20, 1, 10), from=0, to=0.4)
#' 
NULL
#'
#' @rdname Posterior_lambda
#' @export 
dpost_lambda <- function(lambda, a, c, d, S, x, y, ...){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- a+d+y
  dGIB(lambda, a.post, c.post, d.post, S, ...)
}
#'
#' @rdname Posterior_lambda
#' @export 
rpost_lambda <- function(n, a, c, d, S, x, y){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- a+d+y
  rGIB(n, a.post, c.post, d.post, S)
}


#' @name Posterior_mu 
#' @rdname Posterior_mu
#' @title Posterior distribution on the rate in the control group
#' @description Density and random generation for the posterior distribution on 
#' the rate in the control group. The distribution function and the quantile function 
#' are not available.
#' @details The pdf of the posterior distribution of the incidence rate \eqn{\mu} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' 
#' @param mu vector of quantiles 
#' @param a,b non-negative shape and rate parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param T sample size in control group 
#' @param x,y counts in the treated group and control group 
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{dGIB}}
#' 
#' @return \code{dpost_mu} gives the density, and \code{rpost_mu} samples from the distribution.
#' 
#' @note \code{Posterior_mu} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dpost_mu(x, 2, 2, 2, 2, 10, 3, 8), from=0, to=2)
#' lines(density(rpost_mu(1e6, 2, 2, 2, 2, 10, 3, 8)), col="red", lty="dashed")
#' 
NULL
#'
#' @rdname Posterior_mu
#' @export 
dpost_mu <- function(mu, a, b, c, d, T, x, y, ...){
  a.post <- a+x+y
  b.post <- b+T
  c.post <- c+x
  d.post <- a+d+y
  dGIB(mu, a.post, d.post, c.post, b.post, ...)
}
#'
#' @rdname Posterior_mu
#' @export 
rpost_mu <- function(n, a, b, c, d, T, x, y){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- a+d+y
  psi <- rbeta2(n, c.post, d.post, 1)
  return(rgamma(n, a.post, 1+psi)/(T+b))
}

