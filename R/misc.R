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
Gauss2F1 <- function(a,b,c,x){
  if(x>=0 & x<1){
    hyperg_2F1(a,b,c,x)
  }else{
    hyperg_2F1(c-a,b,c,1-1/(1-x))/(1-x)^b
  }
}