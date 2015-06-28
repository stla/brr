#' @name Inference
#' @rdname Inference
#' 
#' @title Inference summaries
#' @description Credibility intervals, estimates
#' @note \code{Inference} is a generic name for the functions documented. 
#' 
#' @param conf confidence level
#' @param intervals a character vector, the intervals to be returned
#' @param ... arguments passed to \link{IntrinsicInference}
#' @return A list of confidence intervals
#' @seealso \code{\link{confint.brr}}
#' 
#' @examples 
#' brr_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0)
#' brr_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0, intervals=c("left","equi-tailed"))
#' brr_estimates(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0)
NULL

#' @rdname Inference
#' @importFrom TeachingDemos hpd
#' @export
brr_intervals <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, conf=.95, intervals="equi-tailed*", ...){
  post.icdf <- function(q){
    qpost_phi(q, a, b, c, d, S, T, x, y)
  }
  hpd2 <- function(x, y, S, T, a, b, c, d, conf){
    if(c+x<1){
      bounds <- c(0, post.icdf(conf))
    }else{
      bounds <- hpd(post.icdf, conf=conf)
    }
    return(bounds)
  }
  bounds <- sapply(intervals, function(interval) setNames(
    switch(interval, 
           left=c(0, post.icdf(conf)), 
           right=c(post.icdf(1-conf), Inf),
           "right*"=c(sign(x)*post.icdf(1-conf), Inf),
           "equi-tailed"=post.icdf(c((1-conf)/2, (1+conf)/2)),
           "equi-tailed*"=c(sign(x)*post.icdf((1-conf)/2),post.icdf((1+conf)/2)),
           hpd=hpd2(x, y, S, T, a, b, c, d, conf), 
           intrinsic=intrinsic_bounds(x, y, S, T, a, b, c, d, conf, ...),
           intrinsic2=if(a==0.5 && b==0) c(NA,NA) else intrinsic2_bounds(x, y, S, T, a, b, c, d, conf, ...)
    ), c("lwr", "upr")
  ), simplify=FALSE)
  if(is.null(bounds)) stop("invalid interval name")
  return(bounds)
}
#'
#' @rdname Inference
#' @export
brr_estimates <- function(x, y, S, T, a=0.5, b=0, c=0.5, d=0, parameter="phi", ...){
  if(!parameter %in% c("phi","VE")) stop("parameter must be 'phi' or 'VE'")
  estimates <- spost_phi(a, b, c, d, S, T, x, y)[c("mode","mean","Q2")]
  names(estimates)[3] <- "median"
  intrinsic <- intrinsic_estimate(x, y, S, T, a, b, c, d, ...)
  out <- c(estimates, intrinsic=intrinsic)
  if(parameter=="VE") out <- lapply(out, function(x) 1-x)
  return(out)
}
