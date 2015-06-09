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
  return(params$a+params$b)
  # faire output ggplots qui s'affichent ou liste de ggplot
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
#' @export
summary.brr <- function(brr){
  params <- brr()
  cat("----------\n")
  cat("*Prior distribution on µ*\n")
  if(all(c("a","b") %in% names(params))){
    cat(with(params, sprintf("Gamma(a=%s,b=%s)\n", a, b)))
  }else{
    cat("a and b not supplied")
  }
  return(invisible())
}