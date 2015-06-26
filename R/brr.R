prior <- function(params){
  #   # remove NULL components
  #   params <- as.list(unlist(params))
  #   if(all(c("a","b","c","d","S","T") %in% names(params))){
  #     return("informative")
  #   }else{
  #     if(!all(c("a","b") %in% names(params))){
  #       return("non-informative")
  #     }else{
  #       if(all(c("a","b") %in% names(params))){
  #         return("semi-informative")
  #       }
  #     }
  #   }
  params <- as.list(unlist(params))
  for(i in seq_along(params)) assign(names(params)[i], params[[i]])
  if((all(c("a","b","c","d") %in% names(params))) && all(c(a,b,c,d)>0)){
    return("informative")
  }
  info_mu <- all(c("a","b") %in% names(params)) && all(c(a,b)>0)
  noninfo_mu <- (is.null(params$a) && is.null(params$b)) || (all(c("a","b") %in% names(params)) && a==0.5 && b==0)
  noninfo_phi <- (is.null(params$c) && is.null(params$d)) || (all(c("c","d") %in% names(params)) && c==0.5 && d==0)
  if(noninfo_mu && noninfo_phi) return("non-informative")
  if(noninfo_phi && info_mu) return("semi-informative")
  stop("Invalid combination of parameters a,b,c,d.")
}

#' @importFrom stringr str_detect
brr_generic <- function(fun, model, parameter, ...){
  # virer les NULL mais mettre les valeurs pour posterior => mettre ça en attribut de prior()
  if(class(model)!="brr") stop("First argument is not a brr class object (see the given example and ?Brr)")
  fun_ <- sprintf("%s_%s", fun, parameter)
  if(! fun_ %in% ls(pos = "package:brr")) stop(sprintf("%s does not exist in brr package.", fun_))
  posterior <- stringr::str_detect(fun, "post")
  args <- formalArgs(fun_) %>% subset(!. %in% "...")
  if(!fun %in% c("sprior", "spost")) args <- args[-1] 
  fun <- eval(parse(text=fun_))
  params <- model()
  if(!posterior && !all(args %in% names(params))) stop(sprintf("Missing parameters. You must supply %s.", 
                                                               paste(args, collapse=", ")))
  if(posterior){
    type <- prior(params)
    if(type=="semi-informative" || type=="non-informative"){
      if(is.null(params$c)) params$c <- 0.5
      if(is.null(params$d)) params$d <- 0
    }
    if(type=="non-informative"){
      params$a <- 0.5; params$b <- 0
    }
    if(!all(args %in% names(params))) stop(sprintf("Missing parameters. You must supply %s (or at least %s, if you want to use the non-informative prior).", 
                                                   paste(args, collapse=", "),
                                                   paste(args[!args %in% c("a","b","c","d")], collapse=", ")))
  }
  return( do.call(fun, c(list(...), params[names(params) %in% args])) )
}


#' @name Brr
#' @rdname Brr
#' @title Brr object
#' @description Set up the Bayesian model and the observations
#' 
#' @param ... prior parameters \code{a}, \code{b}, \code{c}, \code{d}, 
#' samples sizes \code{S}, \code{T}, observed counts \code{x}, \code{y},  
#' future sample sizes \code{Snew}, \code{Tnew}, to be set as in a list (see examples) 
#' @param phi0 the value of interest of the rate ratio
#' @param hypothesis \code{"greater"} to return \eqn{Pr(\phi>\phi_0)}, 
#' \code{"lower"} to return \eqn{Pr(\phi<\phi_0)}  
#' @param table.style the style of the table to print 
#' (passed to  \code{\link[=pander]{pandoc.table.return}})
#' 
#' @return \code{Brr} returns an object of class \code{brr}, \code{summary.brr} 
#' returns a list but prints its contents through \code{print.summary.brr}
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
NULL

#' @rdname Brr
#' @export
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
#'
#' @rdname Brr
#' @export 
summary.brr <- function(brr, phi0=1, hypothesis="greater"){
  out <- list()
  class(out) <- "summary.brr"
  params <- brr()
  type <- prior(params)
  out$type <- type
  # remove NULL components
  params <- as.list(unlist(params)) # pas necessaire je pense - si pour print
  out$params <- params
  if(type!="non-informative"){
    out$prior_mu <- sprior_mu(params$a, params$b)
  }else{
    params$a <- 0.5; params$b <- 0
    out$prior_mu <- "non-informative prior"
  }
  if(all(c("c","d","b","S","T") %in% names(params))){
    out$prior_phi <- with(params, sprior_phi(b, c, d, S, T))
  }else{
    if(type=="non-informative" || type=="semi-informative"){
      params$c <- 0.5; params$d <- 0
      out$prior_phi <- "non-informative prior"
    }
  }
  if(all(c("a","b","c","d","S","T","x","y") %in% names(params))){
    out$post_phi <- with(params, spost_phi(a, b, c, d, S, T, x, y))
    out$Pr <- with(params, ppost_phi(phi0, a, b, c, d, S, T, x, y, 
                                     lower.tail=hypothesis=="greater"))
    out$hypothesis <- hypothesis
    out$phi0 <- phi0
  }else{
    out$Pr <- NULL
  }
return(out)
}
#' 
#' @rdname Brr
#' @importFrom pander pandoc.table.return
#' @export
print.summary.brr <- function(summary, table.style="grid"){
  cat("----------\n")
  cat(type <- summary$type, "prior")
  cat("\n\n")
  cat(sprintf("*Prior distribution on %s*:", greek_utf8("mu")))
  if(type!="non-informative"){
    cat(with(summary$params, sprintf("  Gamma(a=%s,b=%s)", a, b)))
    cat(pandoc.table.return(data.frame(summary$prior_mu), style=table.style))
    #cat(pander(data.frame(summary$prior_mu), style=table.style))
  }else{
    cat("\n  Non-informative prior")
  }
  #cat("\n")
  cat(sprintf("*Prior distribution on %s*:", greek_utf8("phi")))
  if(all(c("c","d","b","S","T") %in% names(summary$params))){
    cat(with(summary$params, sprintf("  Beta2(c=%s,d=%s,scale=%s)", c, d, (T+b)/S)))
    cat(pandoc.table.return(data.frame(summary$prior_phi), style=table.style))
  }else{
    if(type=="non-informative" || type=="semi-informative"){
      cat("  Non-informative prior")
    }else{
      cat("  c, d, b, S and T must be supplied")
    }
    cat("\n\n")
  }
  #cat("\n")
  params <- summary$params
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
  cat(sprintf("*Posterior distribution on %s*:", greek_utf8("phi")))
  if(all(c("a","b","c","d","S","T","x","y") %in% names(params))){
    cat(with(params, sprintf("  Beta2(%s,%s,scale=%s)", c+x, d+a+y, (T+b)/S)))
    cat(pandoc.table.return(data.frame(summary$post_phi), style=table.style))
    #cat("\n")
    cat(sprintf(" Pr('relative risk is %s than %s') = %s",
                summary$hypothesis, summary$phi0, 
                summary$Pr)
    )
  }else{
    cat("\n  a, b, c, d, S, T, x and y must be supplied")
  }
  cat("\n")
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
#' model <- Brr()
#' \dontrun{
#' ppost(model, "phi", 1)}
#' model <- model(x=5, y=10, S=3, T=10)
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



#' plot brr
#' 
#' @examples
#' model <- Brr(a=2, b=3)
#' plot(model)
#' plot(model, dprior(mu))
#' plot(model, pprior(mu))
#' model <- model(c=4, d=6, S=10, T=10)
#' plot(model)
#' plot(model, dprior(phi))
#' plot(model, dprior(x))
#' 
#' @importFrom stringr str_sub
#' @export
plot.brr <- function(brr, what="summary"){ # marche car plot a déjà méthode S3 ; test.brr marche pas : il faudrait définir test() avec UseMethod
  params <- brr()
  type <- prior(params)
  for(i in seq_along(params)) assign(names(params)[i], params[[i]])
  # specific plot 
  if(substitute(what)[1] != "summary"){
    what <- as.character(substitute(what)) # e.g. dprior(phi) => "dprior" "phi"
    f <- what[1]; param <- what[2]
    fun <- paste(what, collapse="_")
    if(! fun %in% ls(pos = "package:brr")) stop(sprintf("%s does not exist in brr package.", fun))
    if(f %in% c("dprior", "pprior", "dpost", "ppost")){
      qfun <- paste0("q", stringr::str_sub(f, 2))
      bounds <- try(eval(parse(text=qfun))(brr, param, c(1e-4, 1-1e-3)), silent=TRUE)
      if(class(bounds)=="try-error"){
        summ <- eval(parse(text=paste0("s", stringr::str_sub(f, 2))))(brr, param)
        bounds <- c(max(0, summ$mean-3*summ$sd), summ$mean+3*summ$sd)
      }
      if(!param %in% c("x","y")){
        seq(bounds[1], bounds[2], length.out=301) %>% {
          plot(., eval(parse(text=f))(brr, param, .), 
               lwd=2, 
               type="l", 
               axes=FALSE,
               ylab=ifelse(f %in% c("dprior","dpost"), NA, sprintf("P(%s \u2264 . )", greek_utf8(param))),
               xlab=ifelse(f %in% c("dprior","dpost"), parse(text=param), NA)
          )
          if(f %in% c("dprior","dpost")) axis(1) else { axis(1); axis(2) }
        }
        return(invisible())
      } else {
        barplot(eval(parse(text=f))(brr, param, bounds[1]:bounds[2]))
        return(invisible())
      }
    } else if(f %in% c("qprior", "qpost")){
      seq(0, 1, length.out=100) %>% {
        plot(., eval(parse(text=f))(brr, param, .), 
             type="l", 
             axes=FALSE, ylab="p", xlab="q")
        axis(1); axis(2)
      }
      return(invisible())
    } else{
      stop("Unvalid 'what' argument.")
    }
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

#' @name inference.brr
#' @rdname inference_brr
#' @title Credibility intervals and estimates
#' @description Get credibility intervals and estimates from a \code{brr} object
#' @details \code{confint.brr} is a wrapper to \code{\link{brr_intervals}}
#' 
#' @param model a \code{\link[Brr]{brr}} object
#' @param conf confidence level
#' @param intervals a character vector, the intervals to be returned
#' @param style the style of the table to print 
#' (passed to  \code{\link[=pander]{pandoc.table.return}})
#' @param ... other aguments passed to \code{\link{brr_intervals}} (currently nothing)
#' 
#' @return \code{confint.brr} returns a list of confidence intervals
#' 
#' @examples 
#' model <- Brr(x=10, y=10, S=100, T=100)
#' confint(model)
NULL

#' @rdname inference_brr
#' @export
confint.brr <- function(model, conf=0.95, intervals="all", ...){
  params <- model()
  type <- prior(params)
  if(type=="semi-informative" || type=="non-informative"){
    params$a <- 0.5; params$b <- 0
    if(type=="non-informative"){
      params$c <- 0.5; params$d <- 0
    }
  }
  args <- formalArgs("brr_intervals") %>% subset(!. %in% "...")
  if(intervals=="all") intervals <- c("equi", "equi.star", "hpd", "intrinsic")
  confints <- do.call(brr_intervals, c(list(...), params[names(params) %in% args], list(conf=conf, intervals=intervals)))
  class(confints) <- "confint.brr"
  attr(confints, "level") <- conf
  return(confints)
}
#'
#' @rdname inference_brr
#' @importFrom pander pandoc.table.return
#' @export
print.confint.brr <- function(x, style="grid"){
  cat(sprintf("%s-credibility intervals about %s", paste0(100*attr(x,"level"),"%"), greek_utf8("phi")))
  cat("\n\n")
  for(i in seq_along(x)){
    cat(names(x)[i], ": ")
    #attr(x[[i]], "caption") <- names(x)[i]
    cat(pandoc.table.return(x[[i]], style=style))
    cat("\n")
  }
}