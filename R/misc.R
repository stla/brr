#' Pipe operator
#' 
#' Pipe operator from \code{\link{magrittr}}
#' 
#' @importFrom magrittr "%>%"
#' @export
'%>%' <- magrittr::`%>%`

#' Unicode encoding of Greek letters
#' 
greek_utf8 <- function(letter){
  return(switch(letter, mu="\u03BC", phi="\u03d5", lambda="\u03BB"))
}

#' Gauss hypergeometric function
#' 
#' @importFrom gsl hyperg_2F1
Gauss2F1 <- function(a, b, c, x){
  if(x>=0 & x<1){
    hyperg_2F1(a,b,c,x)
  }else{
    #exp(log(hyperg_2F1(c-a,b,c,1-1/(1-x))) - b*log1p(-x))
    hyperg_2F1(c-a,b,c,1-1/(1-x)) / (1-x)^b
  }
}

#' Inverse cdf of a discrete distribution
#' 
#' @param pmf a probability mass function
#' @param p probability
#' @param ... arguments passed to \code{pmf}
#'
#' @examples
#' icdf(dpois, 0.5, lambda=10)
#' qpois(0.5, 10)
#'
#' @export 
icdf <- function(pmf, p, ...){
  q <- 0
  prob <- pmf(0, ...)
  while(prob < p){
    q <- q+1
    prob <- prob + pmf(q, ...)
  }
  return(q)
}

#' Moment of a discrete distribution
#' 
#' @param pmf a probability mass function
#' @param k order
#' @param accuracy accuracy
#' @param ... arguments passed to \code{pmf}
#' 
#' @examples
#' dd_moment(dpois, lambda=5)
#' dd_moment(dpois, lambda=5.5795791557050280)
#' dd_moment(dpois, k=2, lambda=5)
#' @export
dd_moment <- function(pmf, k=1, accuracy=.Machine$double.eps, ...){
 m0 <- 0
 x <- 1
 px <- pmf(x, ...)
 m1 <- m0+x^k*px
 while(px==0 || (m1-m0)>accuracy){
   x <- x + 1
   px <- pmf(x, ...)
   m0 <- m1
   m1 <- m0+x^k*px
 }
 return(m1)
}