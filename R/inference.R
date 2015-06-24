#' @name Inference
#' @rdname Inference
#' 
#' @title Inference summaries
#' @description Credibility intervals, estimates
#' @note \code{Inference} is a generic name for the functions documented. 
#' 
#' @param a
#' @return xxx
#' 
#' @examples 
#' credibility_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0)
#' credibility_intervals(x=4, y=5, S=10, T=10, a=0.5, b=0, c=0.5, d=0, intervals=c("left","equi"))
NULL


#' @rdname Inference
#' @export
credibility_intervals <- function(x, y, S, T, a, b, c, d, conf=.95, intervals="equi.star"){
  post.icdf <- function(q){
    qpost_phi(q, a, b, c, d, S, T, x, y)
  }
  hpd2 <- function(x, y, S, T, a, b, c, d, conf){
    post.icdf <- function(q){
      qpost.phi(q, a, b, c, d, S, T, x, y)
    }
    if(c+x<1){
      bounds<-c(0, post.icdf(conf))
    }else{
      bounds <- hpd(post.icdf, conf=conf)
    }
    return(bounds)
  }
  bounds <- sapply(intervals, function(interval) switch(interval, 
                   left=c(0, post.icdf(conf)), 
                   right=c(post.icdf(1-conf), Inf),
                   right.star=c(sign(x)*post.icdf(1-conf), Inf),
                   equi=post.icdf(c((1-conf)/2, (1+conf)/2)),
                   equi.star=c(sign(x)*post.icdf((1-conf)/2),post.icdf((1+conf)/2)),
                   hpd=hpd2(x, y, S, T, a, b, c, d, conf), 
                   intrinsic=bounds.intrinsic(x, y, S, T, a, b, c, d, conf)
  ), simplify=FALSE)
  if(is.null(bounds)) stop("invalid interval name")
  return(bounds)
}
