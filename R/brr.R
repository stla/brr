#' Set up
#'
#'@export
brr <- function(...){
  parameters <- sapply(list(...), identity, simplify=FALSE) #names(list(...))
  params <- c("a","b","c","d","S","T","x","y")
  if((!is.null(names(parameters))) && ! all(names(parameters) %in% params)){
    not <- names(parameters)[which(! names(parameters) %in% params)]
    if(length(not)==1){
      return( stop(sprintf("Invalid parameter '%s'.", not)) )
    }else{
      return( stop(sprintf("Invalid parameters '%s'.", paste(not, collapse=", "))) )
    }
  }
  out <- function(...){
    if(length(list(...))==0){
      return(parameters)
    }else{
#       P <- list(...)
#       return(function(...){
#         if(length(list(...))==0){
#           return(c(parameters, do.call(brr,P)()))
#         }else{
#           return(c(parameters, brr(...)()))
#         }
#       })
    #return(do.call(brr, c(list(...), parameters))) # pb error message : "Erreur dans function()" au lieu de "Erreur dans brr"
    parameters <- c(list(...), parameters)
    return( eval(parse(
      text=sprintf("brr(%s)", paste(sprintf("%s=%s", names(parameters), sapply(parameters, function(x) as.character(x))), collapse=",")))) )
    }
  }
  class(out) <- "brr"
  return(out)
}

#' plot brr
#' 
#' @examples
#' model <- brr(a=2, b=3)
#' plot(model)
#' @export
plot.brr <- function(brr){ # marche car plot a déjà méthode S3 ; test.brr marche pas : il faudrait définir test() avec UseMethod
  params <- brr()
  for(i in seq_along(params)) assign(names(params)[i], params[[i]])
  # prior mu 
  bounds <- qprior_mu(c(1e-4, 1-1e-4), a=a, b=b)
  mu <- seq(bounds[1], bounds[2], length.out=100)
  mu %>% {plot(., dprior_mu(., a=a, b=b), 
            type="l", axes=FALSE, 
            xlab=expression(mu), ylab=NA)}
  axis(1)
  readline(prompt="Press [enter] to continue")
  return(params$a+params$b)
  # faire output ggplots qui s'affichent ou liste de ggplot
}

#' Type of the prior
#' 
#' 
prior <- function(params){
  #params <- brr()
  if(all(c("a","b","c","d","S","T") %in% names(params))){
    return("informative prior")
  }else{
    if(!all(c("a","b") %in% names(params))){
      return("non-informative prior")
    }else{
      if(all(c("a","b") %in% names(params))){
        return("semi-informative prior")
      }
    }
  }
  # warning si c, d mais pas S et T
}

#' Summary brr
#' 
#' @examples
#' model <- brr()
#' summary(model)
#' model <- brr(x=3, y=4)
#' summary(model)
#' model <- brr(a=2, b=4, T=10)
#' summary(model)
#' model <- brr(a=2, b=4, c=3, d=5, S=10, T=10)
#' summary(model)
#' model <- model(x=5, y=10)
#' summary(model)
#' @export
summary.brr <- function(brr){
  params <- brr()
  type <- prior(params)
  cat("----------\n")
  cat(type)
  cat("\n\n")
  cat("*Prior distribution on µ*\n")
  if(all(c("a","b") %in% names(params))){
    cat(with(params, sprintf("  Gamma(a=%s,b=%s)", a, b)))
    summary_gamma(params$a, params$b, type="pandoc", style="rmarkdown")
  }else{
    cat("  Non-informative prior\n")
  }
  cat("\n")
  cat("*Prior distribution on phi*\n")
  if(all(c("c","d","b","S","T") %in% names(params))){
    cat(with(params, sprintf("  Beta2(c=%s,d=%s,scale=%s)", c, d, (T+b)/S)))
    with(params, summary_prior_phi(b, c, d, S, T, type="pandoc", style="rmarkdown"))
  }else{
    if(type=="non-informative prior" || type=="semi-informative prior"){
      cat("  Non-informative prior")
    }else{
      cat("  c, d, b, S and T must be supplied")
    }
    cat("\n")
  }
  cat("\n")
  cat("*Sample sizes*\n")
  cat(sprintf("  S (treated group): %s", ifelse("S" %in% names(params), params$S, "not supplied")))
  cat("\n")
  cat(sprintf("  T (control group): %s", ifelse("T" %in% names(params), params$T, "not supplied")))
  cat("\n\n")
  cat("*Observed counts*\n")
  cat(sprintf("  x (treated group): %s", ifelse("x" %in% names(params), params$x, "not supplied")))
  cat("\n")
  cat(sprintf("  y (control group): %s", ifelse("y" %in% names(params), params$y, "not supplied")))
  cat("\n\n")
  cat("*Posterior distribution on phi*\n")
  if(all(c("a","b","c","d","S","T","x","y") %in% names(params))){
    cat(with(params, sprintf("  Beta2(%s,%s,scale=%s)", c+x, d+a+y, (T+b)/S)))
    with(params, summary_post_phi(a, b, c, d, S, T, x, y, type="pandoc", style="rmarkdown"))
  }else{
      cat("  a, b, c, d, S, T, x and y must be supplied")
  }
  cat("\n")
  return(invisible())
}