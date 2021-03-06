% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/priors.R
\name{Prior_phi}
\alias{Prior_phi}
\alias{dprior_VE}
\alias{dprior_phi}
\alias{pprior_VE}
\alias{pprior_phi}
\alias{qprior_VE}
\alias{qprior_phi}
\alias{rprior_phi}
\alias{sprior_phi}
\title{Prior distribution on the relative risk and the vaccine efficacy}
\usage{
dprior_phi(phi, b, c, d, S, T, ...)

dprior_VE(VE, b, c, d, S, T, ...)

pprior_phi(q, b, c, d, S, T, ...)

pprior_VE(q, b, c, d, S, T, ...)

qprior_phi(p, b, c, d, S, T, ...)

qprior_VE(p, b, c, d, S, T, ...)

rprior_phi(n, b, c, d, S, T)

sprior_phi(b, c, d, S, T, ...)
}
\arguments{
\item{phi,VE,q}{vector of quantiles}

\item{b}{non-negative rate parameter}

\item{c,d}{non-negative shape parameters}

\item{S,T}{sample sizes in control group and treated group}

\item{...}{other arguments passed to \code{\link{Beta2Dist}}}

\item{p}{vector of probabilities}

\item{n}{number of observations to be simulated}
}
\value{
\code{dprior_phi} gives the density, \code{pprior_phi} the distribution
function, \code{qprior_phi} the quantile function, \code{rprior_phi} samples from
the distribution, and \code{sprior_phi} gives a summary of the distribution.
}
\description{
Density, distribution function, quantile function and random
generation for the prior distribution on relative risk or the vaccine efficacy.










}
\details{
The prior distribution on the relative risk \eqn{\phi} is the Beta2 distribution
with shape parameters \eqn{c} and \eqn{d} and scale parameter \eqn{(T+b)/S}.
}
\note{
\code{Prior_phi} is a generic name for the functions documented.
}
\examples{
curve(dprior_phi(x, 2, 2, 2, 10, 10), from=0, to=7)
sprior_phi(2, 2, 2, 10, 10, output="pandoc")
}

