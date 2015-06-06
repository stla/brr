
#' @name Post_x 
#' @rdname Post_x
#' @title Posterior predictive distribution of the count in the treated group
#' @description Density, distribution function, quantile function and random 
#' generation for the posterior predictive distribution of the count in the treated group.
#' @details The posterior predictive distribution of the count in the treated group is a  
#' \code{\link[=GammaInverseBetaDist]{Gamma-Inverse Beta distribution}}.
#' 
#' @param x2,q vector of non-negative \strong{integer} quantiles 
#' @param x1,y1 counts (integer) in the treated group and control group of the observed experiment
#' @param p vector of probabilities
#' @param a non-negative shape parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S1,S2 sample sizes of the control group in the observed experiment and the
#' predicted experiment
#' @param n number of observations to be simulated
#' 
#' @return \code{dpost_x} gives the density, \code{ppost_x} the distribution function, \code{qpost_x} the quantile function, and \code{rpost_x} samples from the distribution.
#' 
#' @note \code{Post_x} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dpost_x(0:10, 10, 2, 3, 4, 5, 3, 4, 10))
#' qpost_x(0.5, 2, 3, 4, 5, 10)
#' ppost_x(5, 2, 3, 4, 5, 10)
#' 
NULL
#'
#' @rdname Post_x
#' @export 
dpost_x <- function(x2, S2, a=0.5, c=0.5, d=0, x1, y1, S1){
  a.post <- a+x1+y1
  c.post <- c+x1
  d.post <- d+a+y1
  return( setNames(dPGIB(x2,a.post,d.post,c.post,S2/S1), paste("x2=",x2,sep="")) )
}
#'
#' @rdname Post_x
#' @export 
ppost_x <- function(q, S2, a=0.5, c=0.5, d=0, x1, y1, S1){
  a.post <- a+x1+y1
  c.post <- c+x1
  d.post <- d+a+y1
  return( setNames(pPGIB(q,a.post,d.post,c.post,S2/S1), paste("x2\u2264",q,sep="")) )
}
#'
#' @rdname Post_x
#' @export 
qpost_x <- function(p, S2, a=0.5, c=0.5, d=0, x1, y1, S1){
  a.post <- a+x1+y1
  c.post <- c+x1
  d.post <- d+a+y1
  return( qPGIB(p,a.post,d.post,c.post,S2/S1) ) 
}
#'
#' @rdname Post_x
#' @export 
rpost_x <- function(n, S2, a=0.5, c=0.5, d=0, x1, y1, S1){
  a.post <- a+x1+y1
  c.post <- c+x1
  d.post <- d+a+y1
  return( rPGIB(n,a.post,d.post,c.post,S2/S1) )
}
