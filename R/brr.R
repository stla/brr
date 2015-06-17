#' Set up
#'
#' @examples
#' model <- Brr(a=2, b=3)
#' model()
#' # add parameters
#' model <- model(c=4, d=5)
#' model() 
#' # replace parameters
#' model <- model(a=10, b=11)
#' model()
#'@export
Brr <- function(...){
  parameters <- sapply(list(...), identity, simplify=FALSE) #names(list(...))
  params <- c("a","b","c","d","S","T","x","y","Snew","Tnew")
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
#           return(c(parameters, do.call(Brr,P)()))
#         }else{
#           return(c(parameters, Brr(...)()))
#         }
#       })
    #return(do.call(Brr, c(list(...), parameters))) # pb error message : "Erreur dans function()" au lieu de "Erreur dans Brr"
    parameters <- c(list(...), parameters[!names(parameters) %in% names(list(...))])
    return( eval(parse(
      text=sprintf("Brr(%s)", paste(sprintf("%s=%s", names(parameters), sapply(parameters, function(x) if(is.null(x)) "NULL" else as.character(x))), collapse=",")))) )
    }
  }
  class(out) <- "brr"
  return(out)
}

#' plot brr
#' 
#' @examples
#' model <- Brr(a=2, b=3)
#' plot(model)
#' model <- model(c=4, d=6, S=10, T=10)
#' plot(model)
#' plot(model, prior(phi))
#' @export
plot.brr <- function(brr, what="summary"){ # marche car plot a déjà méthode S3 ; test.brr marche pas : il faudrait définir test() avec UseMethod
  params <- brr()
  type <- prior(params)
  for(i in seq_along(params)) assign(names(params)[i], params[[i]])
  # specific plot 
  if(substitute(what)!="summary"){
    f <- eval(parse(text=paste(as.character(substitute(what)), collapse="_")))
    return(f)
  }
  
  # summary 
  if(substitute(what)=="summary"){
    # prior mu 
    if(type != "non-informative"){
      bounds <- qprior_mu(c(1e-4, 1-1e-4), a=a, b=b)
      mu <- seq(bounds[1], bounds[2], length.out=100)
      mu %>% {
        plot(., dprior_mu(., a=a, b=b), 
             type="l", axes=FALSE, 
             xlab=expression(mu), ylab=NA, 
             main=expression(paste("Prior distribution of ", mu)) )
      }
      axis(1)
      readline(prompt="Press [enter] to continue")
    }
    # prior phi
    if(type == "informative"){
      bounds <- qprior_phi(c(1e-4, 1-1e-4), b=b, c=c, d=d, S=S, T=T)
      phi <- seq(bounds[1], bounds[2], length.out=100)
      phi %>% {
        plot(., dprior_phi(., b=b, c=c, d=d, S=S, T=T), 
             type="l", axes=FALSE, 
             xlab=expression(phi), ylab=NA, 
             main=expression(paste("Prior distribution of ", phi)) )
      }
      axis(1)
      readline(prompt="Press [enter] to continue")
    }
    # posteriors
    if(!all(c("x","y","S","T") %in% names(params))) return(invisible())
    missings <- NULL
    if(!"a" %in% names(params)){
      a <- 0.5
      missings <- c(missings, "a")
    }
    if(!"b" %in% names(params)){
      b <- 0
      missings <- c(missings, "b")
    }
    if(!"c" %in% names(params)){
      c <- 0.5
      missings <- c(missings, "c")
    }
    if(!"a" %in% names(params)){
      d <- 0
      missings <- c(missings, "d")
    }
    return(missings)
    # faire output ggplots qui s'affichent ou liste de ggplot
  }
}

#' Type of the prior
#' 
#' 
prior <- function(params){
  #params <- brr()
  # remove NULL components
  params <- as.list(unlist(params))
  if(all(c("a","b","c","d","S","T") %in% names(params))){
    return("informative")
  }else{
    if(!all(c("a","b") %in% names(params))){
      return("non-informative")
    }else{
      if(all(c("a","b") %in% names(params))){
        return("semi-informative")
      }
    }
  }
  # warning si c, d mais pas S et T
}

#' Summary brr
#' 
#' @examples
#' model <- Brr()
#' summary(model)
#' model <- Brr(x=3, y=4)
#' summary(model)
#' model <- Brr(a=2, b=4, T=10)
#' summary(model)
#' model <- model(a=2, b=4, c=3, d=5, S=10, T=10)
#' summary(model)
#' model <- model(x=5, y=10)
#' summary(model)
#' @export
summary.brr <- function(brr){
  params <- brr()
  type <- prior(params)
  # remove NULL components
  params <- as.list(unlist(params))
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
    if(type=="non-informative" || type=="semi-informative"){
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

#' Generic function
#' 
brr_generic <- function(fun, model, parameter, ...){
  # virer les NULL mais mettre les valeurs pour posterior => mettre ça en attribut de prior()
  if(class(model)!="brr") stop("First argument is not a brr class object (see the given example and ?Brr)")
  fun <- sprintf("%s_%s", fun, parameter)
  if(! fun %in% ls(pos = "package:brr")) stop(sprintf("%s does not exist in brr package.", fun))
  fun <- eval(parse(text=fun))
  args <- formalArgs(fun) %>% subset(!. %in% "...") %>% .[-1]
  params <- model()
  if(!all(args %in% names(params))) stop(sprintf("Missing parameters. You must supply %s.", 
                                                 paste(args, collapse=", ")))
  return( do.call(fun, c(list(...), params[names(params) %in% args])) )
}

#' @name PriorAndPosterior
#' @rdname PriorAndPosterior
#' @title Prior and posterior distributions
#' @description Generic functions for prior and posterior distributions
#' 
#' @param model an object of class \code{\link[=Brr]{brr}}
#' @param parameter a character string among \code{mu}, \code{phi}, \code{lambda}, \code{x}, \code{y}
#' @param ... the first argument of the function called 
#' 
#' @examples
#' model <- Brr(a=2, b=4)
#' dprior(model, "mu", 1:3)
#' # the same:
#' dprior_mu(mu=1:3, a=2, b=4)
#' \dontrun{
#' dprior(model, "lambda", 1:3)}
#' model <- model(c=4, d=5, S=10, T=10)
#' dprior(model, "lambda", 1:3)
#' model <- model(x=5, y=10)
#' ppost(model, "phi", 1)
NULL 
#' 
#' @rdname PriorAndPosterior
#' @export
dprior <- function(model, parameter, ...){
  return( brr_generic("dprior", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
pprior <- function(model, parameter, ...){
  return( brr_generic("pprior", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
qprior <- function(model, parameter, ...){
  return( brr_generic("qprior", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
rprior <- function(model, parameter, ...){
  return( brr_generic("rprior", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
sprior <- function(model, parameter, ...){
  return( brr_generic("sprior", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
dpost <- function(model, parameter, ...){
  return( brr_generic("dpost", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
ppost <- function(model, parameter, ...){
  return( brr_generic("ppost", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
qpost <- function(model, parameter, ...){
  return( brr_generic("qpost", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
rpost <- function(model, parameter, ...){
  return( brr_generic("rpost", model, parameter, ...) )
}
#' @rdname PriorAndPosterior
#' @export
spost <- function(model, parameter, ...){
  return( brr_generic("spost", model, parameter, ...) )
}
