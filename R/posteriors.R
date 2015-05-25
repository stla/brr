#' @name Posterior_lambda 
#' @rdname Posterior_lambda
#' @title Posterior distribution on the incidence rate in the treated group
#' @description Density and random  generation for the prior distribution on 
#' the rate in the treated group. The distribution function and the quantile function 
#' are not available.
#' @details The pdf of the posterior distribution of the incidence rate \eqn{\lambda} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' 
#' 
#' @param lambda vector of quantiles 
#' @param a,b non-negative shape and rate parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S,T sample sizes in control group and treated group
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
