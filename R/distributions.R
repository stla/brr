#' ## cumulative distribution of B2(c, d, scale)
#'
#'@export
pbeta2<-function(t, c, d, scale, ...){
  pf(d/c/scale*t, df1=2*c, df2=2*d, ...)
}