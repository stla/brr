#' Pipe operator
#' 
#' Pipe operator from \code{\link{magrittr}}
#' 
#' @importFrom magrittr "%>%"
#' @export
'%>%' <- magrittr::`%>%`

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
#' iquantiles(dpois, 0.5, lambda=10)
#' qpois(0.5, 10)
#'
#' @export 
iquantiles <- function(pmf, p, ...){
  q <- 0
  prob <- pmf(0, ...)
  while(prob < p){
    q <- q+1
    prob <- prob + pmf(q, ...)
  }
  return(q)
}