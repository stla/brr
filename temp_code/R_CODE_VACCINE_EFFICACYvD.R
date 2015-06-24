### version 04 nov 2012




## required packages 
library(TeachingDemos) ## package for function hpd()
library(SuppDists)  ## for function ghyper()
library(gsl)  ## for function hyperg_2F1()



####################################################
##### Frequentist confidence intervals #############
####################################################

bounds.SK <- function(conf, x, y, S, T){
	alpha <- 1-conf
	z <- qnorm(alpha, 0, 1, lower.tail=FALSE)
	D <- x+y+1-0.25*z^2
	if(!(x==0 & y==0) & D>0){
		t1 <- sqrt((x+0.5)*(y+0.5))
		t2 <- 0.5*z*sqrt(D)
		d <- y+0.5-0.25*z^2
		LB <- ((t1-t2)/d)^2
		UB <- ((t1+t2)/d)^2		
	}else{
		LB<-0 
		UB <- Inf
	}
	bounds <- T/S*c(LB,UB)
return(bounds)
}


bounds.binomial <- function(conf, x, y, S, T){
	alpha <- 1-conf
	n <- x+y
	p.L <- function(alpha) {
        if (x == 0) 
            0
        else qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(alpha) {
        if (x == n) 
            1
        else qbeta(1 - alpha, x + 1, n - x)
    }
	L <- p.L(alpha)
	U <- p.U(alpha)
	bounds <- T/S*c(L,U)/(1-c(L,U))
return(bounds)
}

Confidence.Interval <- function(x, y, S, T, conf, interval="two.sided", type="binomial"){
	if(interval=="two.sided"){conf<-conf+(1-conf)/2}
	Bounds <- switch(type, 
		binomial=bounds.binomial,
		SK=bounds.SK)
	if(is.null(Bounds)) stop("invalid type name")
	Bounds <- Bounds(conf,x,y,S,T)
	bounds <- switch(interval, 
			left=c(0, Bounds[2]), 
			right=c(Bounds[1], Inf),
			two.sided=Bounds
		)
		if(is.null(bounds)) stop("invalid interval name")
return(bounds)
}



##########################################################
##########################################################
###### B2(c, d, scale) distribution ('Beta of second kind') ########
##  = Beta prime when scale=1      ########################
##########################################################
##########################################################

## cumulative distribution of B2(c, d, scale)
pbeta2<-function(t, c, d, scale, ...){
	pf(d/c/scale*t, df1=2*c, df2=2*d, ...)
}

## quantile function of B2(c, d, scale)
qbeta2 <- function(p, c, d, scale, lower.tail=TRUE){
	scale*c/d*qf(p, df1=2*c, df2=2*d, lower.tail=lower.tail)
}

## density function of B2(c, d, scale)
dbeta2<-function(x, c, d, scale, log=FALSE, ...){
	k <- d/c/scale
	switch(as.character(log),
		"FALSE" = k*df(k*x, df1=2*c, df2=2*d, log=log, ...), 
		"TRUE" = log(k)+df(k*x, df1=2*c, df2=2*d, log=log, ...)
	)
}

## simulation of B2(c, d, scale)
rbeta2<-function(nsims,c, d, scale){
	scale*c/d*rf(nsims, df1=2*c, df2=2*d)
}

## mode of B2(c, d, scale)
beta2.mode <- function(c, d, scale){
	m <- 0
	if(c>1){
		m <- scale*(c-1)/(d+1)
	}
m
}

### mean of B2(c, d, scale)
beta2.mean <- function(c, d, scale){
	mean <- Inf
	if(d>1){
		mean <- scale*c/(d-1)
	}
mean
}

### variance of B2(c, d, scale)
beta2.var <- function(c, d, scale){
	var <- Inf
	if(d>2){
		var <- scale^2*c*(c+d-1)/((d-1)^2*(d-2))
	}
var
}

#####
beta2.summary <- function(c,d,scale){
	list(mode=beta2.mode(c,d,scale),  
		mean=beta2.mean(c,d,scale),  
		var=beta2.var(c,d,scale),  
		Q1=qbeta2(0.25,c,d,scale),
		Q2=qbeta2(0.5,c,d,scale),
		Q3=qbeta2(0.75,c,d,scale)
	)
}



###################################################################
#######  prior distribution of phi when mu ~ Gamma(a,b)          ##
#######  posterior : replace  c  and  d  by  c.post  and  d.post ## 
####### Remark : the prior do not depend on a                    ##
###################################################################


## summary of prior
## ? virer => inclure ?a dans la fonction summary(brrprior) uniquement
summary.prior.phi <- function(b, c, d, S, T){
	beta2.summary(c,d,scale=(T+b)/S)
}


#
## prior density function of phi and VE
dprior.phi<-function(phi, b, c, d, S, T, ...){
	scale <- (T+b)/S
	dbeta2(phi,c,d,scale, ...)
}
dprior.VE<-function(VE, b, c, d, S, T, ...){
	dprior.phi(1-VE, b, c, d, S, T, ...)
}

## prior cumulative distribution of phi and VE
pprior.phi<-function(q, b, c, d, S, T, ...){
	scale <- (T+b)/S
	pbeta2(q,c,d,scale, ...)
}
pprior.VE<-function(q, b, c, d, S, T, ...){
	1-pprior.phi(1-q, b, c, d, S, T, ...)
}

## prior quantile function of phi and VE
qprior.phi<-function(p, b, c, d, S, T, lower.tail=TRUE){
	scale <- (T+b)/S
	qbeta2(p,c,d,scale,lower.tail=lower.tail)
}
qprior.VE <- function(p, b, c, d, S, T, lower.tail=TRUE){
	1-qprior.phi(1-p, b, c, d, S, T, lower.tail)
}

## simulation of prior of phi
rprior.phi<-function(nsims, b, c, d, S, T){
	scale <- (T+b)/S
	rbeta2(nsims,c,d,scale)
}

## prior mode of phi 
mode.prior.phi <- function(b, c, d, S, T){
	scale <- (T+b)/S
	beta2.mode(c,d,scale)
}

## prior mean of phi
mean.prior.phi <- function(b, c, d, S, T){
	scale <- (T+b)/S
	beta2.mean(c,d,scale)
}

## prior variance of phi
var.prior.phi <- function(b, c, d, S, T){
	scale <- (T+b)/S
	beta2.var(c,d,scale)
}


####################################################
########  prior distribution on mu  ############
####################################################

## prior density function of mu 
dprior.mu<-function(mu, a, b, ...){
	dgamma(mu, a, b, ...)
}


## prior cumulative distribution of mu 
pprior.mu<-function(q, a, b, ...){
	pgamma(q, a, b, ...)
}


## prior quantile function of mu and VE
qprior.mu<-function(p, a, b, lower.tail=TRUE){
	qgamma(p, a, b, lower.tail=lower.tail)
}


## simulation of prior of mu
rprior.mu<-function(nsims, a, b){
	rgamma(nsims, a, b, ...)
}


####################################################
########  prior distribution on lambda  ############
####################################################

dprior.lambda <- function(lambda, a, b, c, d, S, T){ # =0 en 0 quand c>1 
	ifelse(c>1 & lambda<.Machine$double.eps, 0, 
	(b*S/(b+T))^a*exp(lnpoch(a,d)-lnbeta(c,d))*lambda^(a-1)*hyperg_U(a+d,a-c+1,b*S/(b+T)*lambda))
}

rprior.lambda <- function(nsims,a,b,c,d,S,T){
	mu.sims <- rgamma(nsims,a,b)
	phi.sims <- rprior.phi(nsims, b, c, d, S, T)
mu.sims*phi.sims
}

#####################################################
##### posterior distributions on lambda and mu  #####
#####################################################

# Gamma-Inverse Beta distribution  
dGIB <- function(x,alpha,beta,gamma,rho){
	exp(lnpoch(beta,gamma)-lngamma(alpha))*rho^alpha*x^(alpha-1)*exp(-rho*x)*hyperg_U(gamma,alpha-beta+1,rho*x)
}


# posterior distribution on mu
dpost.mu <- function(mu, a, b, c, d, T, x, y){
	a.post <- a+x+y
	b.post <- b+T
	c.post <- c+x
	d.post <- a+d+y
dGIB(mu,a.post,d.post,c.post,b.post)
}

# simulates from posterior distribution on mu 
rpost.mu <- function(nsims, a, b, c, d, T, x, y){
	a.post <- a+x+y
	c.post <- c+x
	d.post <- a+d+y
	psi <- rbeta2(nsims,c.post,d.post,1)
rgamma(nsims, a.post, 1+psi)/(T+b)
}


# posterior distribution on lambda
dpost.lambda <- function(lambda, a, c, d, S, x, y){
	a.post <- a+x+y
	c.post <- c+x
	d.post <- a+d+y
dGIB(lambda,a.post,c.post,d.post,S)
}


# simulates from posterior distribution on lambda
rpost.lambda <- function(nsims, a, c, d, S, x, y){
	a.post <- a+x+y
	c.post <- c+x
	d.post <- a+d+y
	psi <- rbeta2(nsims,c.post,d.post,1)
	mu0 <- rgamma(nsims, a.post, 1+psi)
mu0*psi/S
}


###################################################################
#######  posterior distribution of phi when mu ~ Gamma(a,b)      ##
#######                                                          ## 
###################################################################


## posterior cumulative distribution of  phi  and  VE
ppost.phi<-function(q, a, b, c, d, S, T, x, y, ...){
	c.post<-c+x ; d.post<-a+d+y
	pprior.phi(q, b, c.post, d.post, S, T, ...)
}
ppost.VE<-function(q, a, b, c, d, S, T, x, y, ...){
	1-ppost.phi(1-q, a, b, c, d, S, T, x, y, ...)
}

## posterior quantile function of  phi  and  VE
qpost.phi<-function(p, a, b, c, d, S, T, x, y, lower.tail=TRUE){
	c.post<-c+x ; d.post<-a+d+y
	qprior.phi(p, b, c.post, d.post, S, T, lower.tail=lower.tail)
}
qpost.VE<-function(p, a, b, c, d, S, T, x, y, lower.tail=TRUE){
	1-qpost.phi(1-p, a, b, c, d, S, T, x, y, lower.tail=lower.tail)
}


## posterior density function of  phi  and  VE
dpost.phi<-function(phi, a, b, c, d, S, T, x, y, ...){
	c.post<-c+x ; d.post<-a+d+y
	dprior.phi(phi, b, c.post, d.post, S, T, ...)
}
dpost.VE<-function(VE, a, b, c, d, S, T, x, y, ...){
	dpost.phi(1-VE, a, b, c, d, S, T, x, y, ...)
}

## simulation of posterior of  phi  
rpost.phi<-function(nsims, a, b, c, d, S, T, x, y){
	c.post<-c+x ; d.post<-a+d+y
	rprior.phi(nsims, b, c.post, d.post, S, T)
}


## posterior mode of  phi  and  VE 
mode.post.phi <- function(a, b, c, d, S, T, x, y){
	c.post<-c+x ; d.post<-a+d+y
	mode.prior.phi(b, c.post, d.post, S, T) 
}
mode.post.VE <- function(a, b, c, d, S, T, x, y){
	1-mode.post.phi(a, b, c, d, S, T, x, y)
}

## posterior mean of  phi  and  VE 
mean.post.phi <- function(a, b, c, d, S, T, x, y){
	c.post<-c+x ; d.post<-a+d+y
	mean.prior.phi(b, c.post, d.post, S, T) 
}
mean.post.VE <- function(a, b, c, d, S, T, x, y){
	1-mean.post.phi(a, b, c, d, S, T, x, y)
}

## posterior variance of  phi  and  VE 
var.post.phi <- function(a, b, c, d, S, T, x, y){
	c.post<-c+x ; d.post<-a+d+y
	var.prior.phi(b, c.post, d.post, S, T) 
}
var.post.VE <- var.post.phi

## estimators of phi and VE
estimates.phi <- function(a, b, c, d, S, T, x, y){
	mean <- mean.post.phi(a, b, c, d, S, T, x, y)
	median <- qpost.phi(0.5,a, b, c, d, S, T, x, y)
	mode <- mode.post.phi(a, b, c, d, S, T, x, y)
	intrinsic <- intrinsic.estimator(x, y, S, T, a, b, c, d)$estimator
list(mean=mean,median=median,mode=mode,intrinsic=intrinsic)
}
estimates.VE <- function(a, b, c, d, S, T, x, y){
	estims.phi <- estimates.phi(a, b, c, d, S, T, x, y)
lapply(estims.phi, function(x){1-x})
}
	



############################################################################
##############  Mixtures priors & posteriors on  phi  ######################
############################################################################
dmprior.phi <- function(phi, b, w, c1, d1, c2, d2, S, T){
	if(any(c(c1,d1,c2,d2)<=0)){stop("c1,d1,c2,d2 must be >0")}
	w*dprior.phi(phi,b,c1,d1,S,T)+(1-w)*dprior.phi(phi,b,c2,d2,S,T)
}

pmprior.phi <- function(q, b, w, c1, d1, c2, d2, S, T){
	if(any(c(c1,d1,c2,d2)<=0)){stop("c1,d1,c2,d2 must be >0")}
	w*pprior.phi(q,b,c1,d1,S,T)+(1-w)*pprior.phi(q,b,c2,d2,S,T)
}
dmpost.phi <-  function(phi, a, b, w, c1, d1, c2, d2, S, T, x, y){
	if(any(c(c1,d1,c2,d2)<=0)){stop("c1,d1,c2,d2 must be >0")}
	log.wpost1 <- log(w)+lbeta(c1+x,a+d1+y)-lbeta(c1,d1)
	log.wpost2 <- log(1-w)+lbeta(c2+x,a+d2+y)-lbeta(c2,d2)
	C <- exp(log.wpost1)+exp(log.wpost2)
	log.f1 <- log.wpost1+dpost.phi(phi,a,b,c1,d1,S,T,x,y, log=TRUE)
	log.f2 <- log.wpost2+dpost.phi(phi,a,b,c2,d2,S,T,x,y, log=TRUE)
	out <- 1/C*(exp(log.f1)+exp(log.f2))
	mostattributes(out) <- list(posterior.w=exp(log.wpost1))
out
}
pmpost.phi <-  function(q, a, b, w, c1, d1, c2, d2, S, T, x, y){
	if(any(c(c1,d1,c2,d2)<=0)){stop("c1,d1,c2,d2 must be >0")}
	log.wpost1 <- log(w)+lbeta(c1+x,a+d1+y)-lbeta(c1,d1)
	log.wpost2 <- log(1-w)+lbeta(c2+x,a+d2+y)-lbeta(c2,d2)
	C <- exp(log.wpost1)+exp(log.wpost2)
	log.F1 <- log.wpost1+ppost.phi(q,a,b,c1,d1,S,T,x,y, log=TRUE)
	log.F2 <- log.wpost2+ppost.phi(q,a,b,c2,d2,S,T,x,y, log=TRUE)
	out <- 1/C*(exp(log.F1)+exp(log.F2))
	mostattributes(out) <- list(posterior.w=exp(log.wpost1))
out
}
























###########################################################
###########################################################
###########################################################
###########################################################
######### Prior and Posterior Predictive Probabilities ####
###########################################################
###########################################################
###########################################################
###########################################################

### Related function: Gauss hypergeometric ###
Gauss2F1 <- function(a,b,c,x){
	if(x>=0 & x<1){
		hyperg_2F1(a,b,c,x)
	}else{
			hyperg_2F1(c-a,b,c,1-1/(1-x))/(1-x)^b
		}
}

###########################################################
############### Related distributions #####################
###########################################################

###  Bivariate Negative Binomial   ####
# this is the distribution of  (x,y)  mid  phi
# obtained by integrating over the prior  mu ~ Gamma(a,b)
bivariate.nbinom <- function(phi,x,y,S,T,a,b){
	dnbinom(y, a, b/(b+T))*dnbinom(x,y+a,(T+b)/(T+b+phi*S))
}

### Beta Negative Binomial ###
## if X|psi ~ PG(a,psi) and psi ~ Beta'(c,d)
## then  X ~ BetaNegBin(a,c,d)
dbeta.nbinom <- function(x,a,c,d,log=FALSE){
	dghyper(x,-d,-a,c-1,log=log)
}
pbeta.nbinom <- function(q,a,c,d){
	pghyper(q,-d,-a,c-1)
}
qbeta.nbinom <- function(p,a,c,d){
	qghyper(p,-d,-a,c-1)
}
rbeta.nbinom <- function(nsims,a,c,d){
	rghyper(nsims,-d,-a,c-1)
}
summary.beta.nbinom <- function(a,c,d){
	sghyper(-d,-a,c-1)
}

### Poisson-Gamma-Beta2 distribution ###
PGB2 <- function(x,a,c,d,tau){
	tau^(a)*dbeta.nbinom(x,a,c,d)*Gauss2F1(a+c,a+x,a+x+c+d,1-tau)
}
rPGB2 <- function(nsims, a, c, d, tau){
	psi <- rbeta2(nsims, c, d, tau)
rnbinom(nsims, a, psi/(1+psi))
}

## Hyperscaled Poisson-Gamma-Inverse Beta ##
PGIB <- function(x,a,c,d,tau){  # cas particulier quand tau=1 ?? je crois que non
	out <- sapply(x, function(x){
		ccpoch <- exp(lnpoch(a,x)+lnpoch(d,x)-lnpoch(c+d,x)-lnfact(x))
		ccpoch*tau^x*
		 Gauss2F1(d+x, a+x, c+d+x, -tau)
		})
return(out)
}

rPGIB <- function(nsims,a,c,d,tau){
	psi <- rbeta2(nsims,c,d,1)
	rnbinom(nsims, a, (1+psi)/(tau+1+psi))
}

###########################################################
######### Prior and Posterior Predictive Probabilities ####
###########################################################


## prior predictive probability to observe y placebo cases
dprior.y <- function(y, T, a, b){
	out <- sapply(y, function(y){
			prob <- dnbinom(y, a, b/(b+T))
		return(prob)
	})
	names(out) <- paste("y=",y,sep="")
return(out)
}
#
pprior.y <- function(y, T, a, b){
	out <- sapply(y, function(y){
			prob <- pnbinom(y, a, b/(b+T))
		return(prob)
	})
	names(out) <- paste("y=",y,sep="")
return(out)
}
#
qprior.y <- function(q, T, a, b){
	out <- sapply(q, function(q){
			p <- qnbinom(q, a, b/(b+T))
		return(p)
	})
return(out)
}

## prior predictive probability to observe x vaccine cases
dprior.x <- function(x, T, a, b, c, d){
	tau <- b/(b+T)
        out <- sapply(x, function(x){
                        prob <- PGB2(x,a,d,c,tau)
                return(prob)
        })
        names(out) <- paste("x=",x,sep="")
return(out)
}
## simulation of prior predictive of x
rprior.x <- function(nsims, T, a, b, c, d){
	tau <- b/(b+T)
rPGB2(nsims,a,d,c,tau)
}
## cumulative prior predictive of x
pprior.x <- function(x, T, a, b, c, d){
	sum(dprior.pred.x(0:x, T, a, b, c, d))
}

#### Prior Predictive of x mid y ################
dprior.x.mid.y <- function(x, y, a, b, c, d){
        out <- sapply(x, function(x){
                        prob <-dbeta.nbinom(x,a+y,d,c)
                return(prob)
        })
        names(out) <- paste("x=",x,sep="")
return(out)
}
#
pprior.x.mid.y <- function(q, y, a, b, c, d){
        out <- sapply(q, function(q){
                        prob <- pbeta.nbinom(q,a+y,d,c)
                return(prob)
        })
        names(out) <- paste("q=",q,sep="")
return(out)
}
#
qprior.x.mid.y <- function(p, y, a, b, c, d){
        out <- sapply(p, function(p){
                        prob <- qbeta.nbinom(p,a+y,d,c)
                return(prob)
        })
        names(out) <- paste(100*p,"%",sep="")
return(out)
}

## ? virer => inclure ?a dans la fonction summary(brrprior) uniquement
summary.prior.x.mid.y <- function(y, a, c, d){
	summary.beta.nbinom(a+y,d,c)
}



## prior predictive probability to observe  x and y  simultaneously
dprior.xy <- function(x, y, T, a, b, c, d){
	out <- sapply(x, function(x){
			sapply(y, function(y){
				dprior.y(y, T, a, b)*dprior.x.mid.y(x, y, a, b, c, d)
			})
		})
	out <- matrix(out,ncol=length(x))
	dimnames(out) <- list(paste("y=",y,sep=""),paste("x=",x,sep=""))
return(out)
}

## simulates from the joint prior predictive
rprior.xy <- function(nsims, a, b, c, d, S, T){
	phi <- rprior.phi(nsims, b, c, d, S, T)
	y <- rnbinom(nsims, a, b/(b+T))
	x <- rnbinom(nsims, a+y, (b+T)/(b+T+phi*S))
list(x=x, y=y)
}

#############################################
###### Posterior Predictive #####
################################################
dpost.x <- function(x2, S2, a=0.5, c=0.5, d=0, x1, y1, S1){
        a.post <- a+x1+y1
        c.post <- c+x1
        d.post <- d+a+y1
        out <- sapply(x2, function(x2){
                        prob <- PGIB(x2,a.post,d.post,c.post,S2/S1)
                return(prob)
        })
        names(out) <- paste("x2=",x2,sep="")
return(out)
}

rpost.x <- function(nsims, S2, a, c, d, x1, y1, S1){ # utilise rPGIB
	a.post <- a+x1+y1
	c.post <- c+x1 ; d.post <- a+d+y1
	psi <- rbeta2(nsims,c.post,d.post,1)
	mu0 <- rgamma(nsims, a.post, 1+psi)
	lambda <- mu0*psi/S1
rpois(nsims, lambda*S2)
}




dpost.y <- function(y2, T2, a=0.5, b=0, c=0.5, d=0, x1, y1,  T1){
        a.post <- a+x1+y1
        c.post <- c+x1
        d.post <- d+a+y1
        out <- sapply(y2, function(y2){
                        prob <-  PGIB(y2,a.post,c.post,d.post,T2/(b+T1))
                return(prob)
        })
        names(out) <- paste("y2=",y2,sep="")
return(out)
}


rpost.y <- function(nsims, T2, a=0.5, b=0, c=0.5, d=0, x1, y1, T1){  # utilise rPGIB
	a.post <- a+x1+y1
	c.post <- c+x1 ; d.post <- a+d+y1
	psi <- rbeta2(nsims,c.post,d.post,1)
	mu <- rgamma(nsims, a.post, 1+psi)/(b+T1)
rpois(nsims, mu*T2)
}


dpost.xy <- function(x2, y2, S2, T2, a, b, c, d, x1, y1, S1, T1){
	aa <- a+x1+y1
	cc <- c+x1
	dd <- d+a+y1
	RS <- S2/S1; RT <- T2/T1
	out <- sapply(x2, function(x2){
			sapply(y2, function(y2){
					ccpoch <- exp(lnpoch(aa,x2+y2)+lnpoch(cc,x2)+lnpoch(dd,y2)-lnpoch(cc+dd,x2+y2)-lnfact(x2)-lnfact(y2))
					prob <- ccpoch*RS^x2*RT^y2*
					Gauss2F1(cc+x2,aa+x2+y2,cc+dd+x2+y2,1-((1+RS)/(1+RT)))/(1+RT)^(aa+x2+y2)
				return(prob)
			})
		})
	out <- matrix(out,ncol=length(x2))
	dimnames(out) <- list(paste("y2=",y2,sep=""),paste("x2=",x2,sep=""))
return(out)
}

## posterior predictive of  x mid y
dpost.x.mid.y <- function(x2, y2, S2, T2, a, b, c, d, x1, y1, S1, T1){
	dpost.xy(x2, y2, S2, T2, a, b, c, d, x1, y1, S1, T1)/dpost.y(y2, T2, a, b, c, d, x1, y1, T1)
}


## simulates from the joint posterior predictive
rpost.xy <- function(nsims, S2, T2, a, b, c, d, x1, y1, S1, T1){
	phi.post <- rpost.phi(nsims, a, b, c, d, S1, T1, x1, y1)
	a.post <- a+x1+y1
	b.post <- phi.post*S1+T1+b
	y2 <- rnbinom(nsims, a.post, b.post/(b.post+T2))
	x2 <- rnbinom(nsims, a.post+y2, (b.post+T2)/(b.post+T2+phi.post*S2))
list(x2=x2, y2=y2)
}

## posterior predictive correlation
cor.BPGIB <- function(a,c,d,tau,t){
	sqrt(tau*t*c*d)*(c+d-a)/sqrt((c+d)*(c+d+1)+tau*(d*(a+c+1)+c*(c+1)))/sqrt((c+d)*(c+d+1)+t*(c*(a+d+1)+d*(d+1)))
}

cor.post.xy <- function(S2, T2, a, b, c, d, x1, y1, S1, T1){
 	a.post <- a+x1+y1
	c.post <- c+x1
	d.post <- d+a+y1
	tau <- S2/S1
	t <- T2/(T1+b)
cor.BPGIB(a.post,c.post,d.post,tau,t)
}

############################################################
####### Conditional Prior and Posterior Predictive  ########
############################################################

### conditional prior predictive 
dprior.y.mid.phi <- function(phi, y, T, a, b){
	out <- sapply(y, function(y){
			prob <- dnbinom(y, a, b/(b+T))
		return(prob)
	})
	names(out) <- paste("y=",y,sep="")
return(out)
}

dprior.x.mid.phi<- function(phi, x, S, a, b){
	out <- sapply(x, function(x){
			prob <- dnbinom(x,a,b/(b+phi*S))
		return(prob)
	})
	names(out) <- paste("x=",x,sep="")
return(out)
}

dprior.xy.mid.phi <- function(phi,x,y,S,T,a,b){
	out <- sapply(x, function(x){
			sapply(y, function(y){
				bivariate.nbinom(phi,x,y,S,T,a,b)
			})
		})
	out <- matrix(out,ncol=length(x))
	dimnames(out) <- list(paste("y=",y,sep=""),paste("x=",x,sep=""))
return(out)
}

###### Conditional Posterior Predictive #######

## conditional posterior predictive probability to observe y placebo cases
dpost.y.mid.phi <- function(phi, y2, T2, a, b, x1, y1, S1, T1){
	a.post <- a+x1+y1
	b.post <- phi*S1+T1+b
	out <- sapply(y2, function(y2){
			prob <- dnbinom(y2,a.post,b.post/(b.post+T2))
		return(prob)
	})
	names(out) <- paste("y2=",y2,sep="")
return(out)
}

## conditional posterior predictive probability to observe x vaccine cases
dpost.x.mid.phi <- function(phi, x2, S2, a, b, x1, y1, S1, T1){
	a.post <- a+x1+y1
	b.post <- phi*S1+T1+b
	out <- sapply(x2, function(x2){
			prob <- dnbinom(x2,a.post,b.post/(b.post+phi*S2))
		return(prob)
	})
	names(out) <- paste("x2=",x2,sep="")
return(out)
}

## conditional posterior predictive probability to observe  x and y  simultaneously
dpost.xy.mid.phi <- function(phi, x2, y2, S2, T2, a, b, x1, y1, S1, T1){
	a.post <- a+x1+y1
	b.post <- phi*S1+T1+b
	out <- sapply(x2, function(x2){
			sapply(y2, function(y2){
				prob <- bivariate.nbinom(phi,x2,y2,S2,T2,a.post,b.post)
			return(prob)
			})
		})
	out <- matrix(out,ncol=length(x2))
	dimnames(out) <- list(paste("y2=",y2,sep=""),paste("x2=",x2,sep=""))
return(out)
}
















################################################################
################################################################
######### CREDIBILITY INTERVALS AND ASSOCIATED TEST ############
################################################################
################################################################

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

# brr.posterior <- function(prior, S, T, x, y)

Credibility.Interval <- function(x, y, S, T, a, b, c, d, conf, interval="equi.star"){
	post.icdf <- function(q){
		qpost.phi(q, a, b, c, d, S, T, x, y)
	}
	bounds <- switch(interval, 
			left=c(0, post.icdf(conf)), 
			right=c(post.icdf(1-conf), Inf),
			right.star=c(sign(x)*post.icdf(1-conf), Inf),
			equi=post.icdf(c((1-conf)/2, (1+conf)/2)),
			equi.star=c(sign(x)*post.icdf((1-conf)/2),post.icdf((1+conf)/2)),
			hpd=hpd2(x, y, S, T, a, b, c, d, conf), 
			intrinsic=bounds.intrinsic(x, y, S, T, a, b, c, d, conf)  # mettre l'intrins?que ? part
		)
		if(is.null(bounds)) stop("invalid interval name")
return(bounds)
}

switch.interval <- function(interval){
	switch(interval,
		left="right",
		left.star="right.star",
		right="left",
		equi="equi",
		equi.star="equi.star",
		hpd="hpd",
		intrinsic="intrinsic"
	)
}

Credibility.Interval.VE <- function(x, y, S, T, a, b, c, d, conf, interval="equi.star"){
	interval.phi <- switch.interval(interval)
	if(is.null(interval.phi)) stop("invalid interval name")
	bounds <- Credibility.Interval(x, y, S, T, a, b, c, d, conf, interval.phi)
	bounds <- sort(1-bounds)
return(bounds)
}	 

#Bayesian.test <- function(x,y,q,phi.star,interval,S,T,a,b,c,d){
#	bounds <- Credibility.Interval(x,y,S,T,a,b,c,d,conf=q,interval)
#	out <- if(phi.star<bounds[1] | phi.star>bounds[2]) 1 else 0
#return(out)
#}	

Bayesian.test <- function(x,y,q,phi.star,interval,S,T,a=0.5,b=0,c=0.5,d=0){
	out <- sapply(x, function(x){
		sapply(y, function(y){  
			bounds <- Credibility.Interval(x,y,S,T,a,b,c,d,conf=q,interval)
			test <- if(phi.star<bounds[1] | phi.star>bounds[2]) 1 else 0
			return(test)
			})
		})
	out <- matrix(out,ncol=length(x))
  dimnames(out) <- list(paste("y=",y,sep=""),paste("x=",x,sep=""))
return(t(out))
}	


####################################################################
####################################################################
####################################################################
####################################################################
##################  FREQUENTIST CHARACTERISTICS ####################
####################################################################
####################################################################
####################################################################
####################################################################

####################################################################
####################################################################
##################          POWER CURVE          ###################
####################################################################
#######  power = P(I(X,Y) do not contain phi0 | phi0, mu0  #########
####################################################################
#
power.curve <- function(K, mu0, phi.star, interval, seq.phi0, S, T, a=0.5, b=0, c=0.5, d=0, prec=1e-5, digits=4){
	lambda.max <- max(seq.phi0)*mu0
	upper.x <- qpois(1-prec/2, lambda.max*S)
	upper.y <- qpois(1-prec/2, mu0*T)
	bounds.array <- array(NA, dim=c(upper.x+1, upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x, y, S, T, a, b, c, d, conf=K, interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.phi0)
	Reject <- rep(0,l)
	for(j in 1:l){
		phi0 <- seq.phi0[j]
		lambda0 <- phi0*mu0
		upper.x <- qpois(1-prec/2, lambda0*S)
		for(x in 0:upper.x){
			for(y in 0:upper.y){
				bounds <- bounds.array[x+1, y+1,]
				if(phi.star<bounds[1] | phi.star>bounds[2]){
					Reject[j] <- Reject[j] + dpois(x, S*lambda0)*dpois(y, T*mu0)
				}
			}
		}
	}
list(phi0=seq.phi0, power=round(Reject, digits))
}






############################################################################## # ne pas mettre ?a
############## SIGNIFICANCE LEVEL  ###########################################
significance.curve <- function(K, seq.mu0, phi.star, interval, S, T, a=0.5, b=0, c=0.5, d=0, prec=1e-5, digits=4){
	mu0.max <- max(seq.mu0)
	lambda.max <- phi.star*mu0.max
	upper.x <- qpois(1-prec/2, lambda.max*S)
	upper.y <- qpois(1-prec/2, mu0.max*T)
	bounds.array <- array(NA, dim=c(upper.x+1, upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x, y, S, T, a, b, c, d, conf=K, interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.mu0)
	Reject <- rep(0,l)
	names(Reject) <- paste("mu0=", seq.mu0, sep="")
	phi0 <- phi.star
	for(j in 1:l){
		mu0 <- seq.mu0[j]
		lambda0 <- phi0*mu0
		upper.x <- qpois(1-prec/2, lambda0*S)
		for(x in 0:upper.x){
			for(y in 0:upper.y){
				bounds <- bounds.array[x+1, y+1,]
				if(phi.star<bounds[1] | phi.star>bounds[2]){
					Reject[j] <- Reject[j] + dpois(x, S*lambda0)*dpois(y, T*mu0)
				}
			}
		}
	}
list(mu0=seq.mu0, alpha=round(Reject, digits))
}


####################################################
####################################################
### COVERAGE ARRAY #################################
####################################################
coverage.array <- function(K,seq.psi0, seq.m0, interval){
	l <- length(seq.psi0); ll <- length(seq.m0)
	array <- NULL
	for(k in 1:l){
		psi0 <- seq.psi0[k]
		row <- significance.curve(K, seq.m0, psi0, interval, S=1, T=1)$alpha
		array <- rbind(array, row)
	}
	colnames(array) <- seq.m0
	rownames(array) <- seq.psi0
return(1-array)
}



####################################################
####################################################
### COVERAGE CURVE #################################
####################################################
coverage.curve <- function(K, mu0, interval, seq.phi0, S, T, a, b, c, d, prec=1e-5, digits=4){
	lambda.max <- max(seq.phi0)*mu0
	upper.x <- qpois(1-prec/2, lambda.max*S)
	upper.y <- qpois(1-prec/2, mu0*T)
	bounds.array <- array(NA, dim=c(upper.x+1, upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x, y, S, T, a, b, c, d, conf=K, interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.phi0)
	Coverage <- rep(0,l)
	for(j in 1:l){
		phi0 <- seq.phi0[j]
		lambda0 <- phi0*mu0
		upper.x <- qpois(1-prec/2, lambda0*S)
		for(x in 0:upper.x){
			for(y in 0:upper.y){
				bounds <- bounds.array[x+1, y+1,]
				if(!(phi0<bounds[1] | phi0>bounds[2])){
					Coverage[j] <- Coverage[j] + dpois(x, S*lambda0)*dpois(y, T*mu0)
				}
			}
		}
	}
list(phi0=seq.phi0, coverage=round(Coverage, digits))
}






######################################################################
#################     Marginal Power Curve     #######################
######################################################################

mpower.curve <- function(K,a0,b0,phi.star,interval,seq.phi0,S,T,a,b,c=0.5,d=0,prec=1e-5,digits=4){
	phi0.max <- max(seq.phi0)
	upper.x <- qnbinom(1-prec/2, a0, b0/(b0+phi0.max*S)) 
	upper.y <- qnbinom(1-prec/2, a0, b0/(b0+T)) 
	bounds.array <- array(NA, dim=c(upper.x+1,upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x,y,S,T,a,b,c,d,conf=K,interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.phi0)
	Reject <- rep(0,l)
	for(j in 1:l){
		phi0 <- seq.phi0[j]
		upper.x <-  qnbinom(1-prec/2, a0, b0/(b0+phi0*S))
		for(x in 0:upper.x){
			for(y in 0:upper.y){
				bounds <- bounds.array[x+1, y+1,]
				if(phi.star<bounds[1] | phi.star>bounds[2]){
					Reject[j] <- Reject[j] + bivariate.nbinom(phi0,x,y,S,T,a0,b0)
				}
			}
		}
	}
list(phi0=seq.phi0, power=round(Reject, digits))
}

####################################################
####################################################
### MARGINAL COVERAGE CURVE ########################
####################################################
mcoverage.curve <- function(K, a0, b0, interval, seq.phi0, S, T, a, b, c, d, prec=1e-5, digits=4){
	phi0.max <- max(seq.phi0)
	upper.x <- qnbinom(1-prec/2, a0, b0/(b0+phi0.max*S)) 
	upper.y <- qnbinom(1-prec/2, a0, b0/(b0+T)) 
	bounds.array <- array(NA, dim=c(upper.x+1, upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x, y, S, T, a, b, c, d, conf=K, interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.phi0)
	Coverage <- rep(0,l)
	for(j in 1:l){
		phi0 <- seq.phi0[j]
		upper.x <-  qnbinom(1-prec/2, a0, b0/(b0+phi0*S))
		for(x in 0:upper.x){
			for(y in 0:upper.y){
				bounds <- bounds.array[x+1, y+1,]
				if(!(phi0<bounds[1] | phi0>bounds[2])){
					Coverage[j] <- Coverage[j] + bivariate.nbinom(phi0,x,y,S,T,a0,b0)
				}
			}
		}
	}
list(phi0=seq.phi0, coverage=round(Coverage, digits))
}





################################################
#########  Prior Predictive Power     ##########
################################################
# Computes the prior predictive power using a,b,c,d as hyperparameters 
# and a0,b0,c0,d0 as parameters for sampling distribution
prior.predictive.power <- function(q,phi.star,interval,S,T,a,b,c,d,a0=NULL,b0=NULL,c0=NULL,d0=NULL,nsims=10000){
	if(is.null(a0)){a0<-a; b0<-b; c0<-c; d0<-d}
	sims.xy <- rprior.xy(nsims, a0, b0, c0, d0, S, T)
	sims.decision <- rep(NA,nsims)
	for(i in 1:nsims){
		x <- sims.xy$x[i]; y <- sims.xy$y[i]
		sims.decision[i] <- Bayesian.test(x,y,q,phi.star,interval,S,T,a,b,c,d)
	}
return(mean(sims.decision))
}


######################################################
#########    Posterior Predictive Power     ##########
######################################################
posterior.predictive.power <- function(q,phi.star,interval,x1,y1,S1,T1,S2,T2,a,b,c,d,nsims=10000){
	decision <- function(x,y,S,T){
		Bayesian.test(x,y,q,phi.star,interval,S,T,a,b,c,d)
	}
	sims.x2y2 <- rpost.xy(nsims,S2,T2,a,b,c,d,x1,y1,S1,T1)
	sims.decision <- rep(NA,nsims)
	for(i in 1:nsims){
		x2 <- sims.x2y2$x2[i]; y2 <- sims.x2y2$y2[i]
		sims.decision[i] <- decision(x1+x2,y1+y2,S1+S2,T1+T2)
	}
return(mean(sims.decision))
}

############################################################
#########    Exact Posterior Predictive Power     ##########
#########    when S1/S2=(T1+b)/T2                 ##########
#########    and S1=T1+b  (hence S2=T2)           ##########
############################################################
mean.PGIB <- function(a,c,d,tau){
	tau*a*d/(c+d)
}

# the bivariate PGIB when t=tau
BPGIB11 <- function(x,y,a,c,d,tau){
	out <- array(NA, dim=c(y+1,x+1))
	colnames(out) <- paste("x=",0:x,sep="")
	rownames(out) <- paste("y=",0:y,sep="")
	out[1,1] <- 1/(tau+1)^a
	for(i in 1:y){
		out[(i+1),1] <- out[i,1]*(a+i-1)*(d+i-1)/i/(c+d+i-1)*tau/(tau+1)
	}
	for(j in 1:x){
		out[,(j+1)] <- out[,j]*tau/(tau+1)*(c+j-1)/j*(a+j-2+c(1:(y+1)))/(c+d+j-2+c(1:(y+1)))
	}
out
}

# returns a region [0,x]*[O,y] with y=r*x
# and with BPGIB11 probability > conf
region.BPGIB11 <- function(conf, r, a, c, d, tau){  
	y <- 2*r-1
	out <- NA
	out[1] <- 1/(tau+1)^a
	for(i in 1:y){
		out[(i+1)] <- out[i]*(a+i-1)*(d+i-1)/i/(c+d+i-1)*tau/(tau+1)
	}
	x <- 1
	out. <- out*(c+x-1)/x*(a+x-2+c(1:(y+1)))/(c+d+x-2+c(1:(y+1)))*tau/(tau+1)	
	out <- cbind(out,out.)
	S <- sum(out)
while(S<0.99){
	y.next <- y+r
	for(i in (y+1):y.next){
		temp <- out[dim(out)[1],] # derni?re ligne de out[]
		next.line <- temp*(d+i-1)/i*(a+i-2+c(1:(x+1)))/(c+d+i-2+c(1:(x+1)))*tau/(tau+1)
		out <- rbind(out,next.line)
		S <- S+sum(next.line)
	}
	y <- y.next
	x <- x+1
	temp <- out[,dim(out)[2]] 
	next.col <-temp*(c+x-1)/x*(a+x-2+c(1:(y+1)))/(c+d+x-2+c(1:(y+1)))*tau/(tau+1)	
	out <- cbind(out,next.col)
	S <- S+sum(next.col)		
}
	colnames(out) <- paste("x=",0:x,sep="")
	rownames(out) <- paste("y=",0:y,sep="")
return(list(probs=out, xy=c(x,y), P=S)) 	
}


# tau=S2/S1=T2/(T1+b)
# and assuming S1=T1+b  (hence S2=T2)
posterior.predictive.power2 <- function(q,phi.star,interval,x1,y1,tau,a,b,c,d,prec=1e-5,digits=4){
	### determines a region with post. pred. probability > q for the higher tau ###
		c.post<-c+x1 ; d.post<-a+d+y1 ; a.post <- a+x1+y1
	# determines r to get a region [0,x]*[O,y] with y=r*x
	mean.x2 <- mean.PGIB(a.post,d.post,c.post,tau)
	mean.y2 <-  mean.PGIB(a.post,c.post,d.post,tau)
	r <- max(1,floor(mean.y2/mean.x2))  
	region <- region.BPGIB11(1-prec, r, a.post, c.post, d.post,max(tau))
	upper.xy <- region$xy
	upper.x <- upper.xy[1]; upper.y <- upper.xy[2]
	bounds.array <- array(NA, dim=c(upper.x+1,upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			# here we use the assumption S1=T1+b
			bounds <- Credibility.Interval(x1+x,y1+y,S=1,T=1,a,b=0,c,d,conf=q,interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(tau)
	Reject <- rep(0,l)
	for(j in 1:l){
		Probs <- BPGIB11(upper.x,upper.y,a.post,c.post,d.post,tau[j])
		for(x2 in 0:upper.x){
			for(y2 in 0:upper.y){
				bounds <- bounds.array[x2+1, y2+1,]
				if(phi.star<bounds[1] | phi.star>bounds[2]){
					prob <- Probs[y2+1,x2+1]
					Reject[j] <- Reject[j] + prob
				}
			}
		}
	}
list(tau=tau, power=round(Reject, digits)) 
}





################################################
######## Conditional Prior Predictive Power #### 
################################################ 

prior.predictive.power.mid.phi <- function(q,phi.star,interval,seq.phi0,S,T,a,b,c,d){
	mpower.curve(q,a,b,phi.star,interval,seq.phi0,S,T,a,b,c,d)
}


####################################################
######## Conditional Posterior Predictive Power ####
####################################################

posterior.predictive.power.mid.phi <- function(K,phi.star,interval,seq.phi0,x1,y1,S1,T1,S2,T2,a,b,c,d,prec=1e-5,digits=4){
	a.post <- a+x1+y1
	b.post <- function(phi) phi*S1+T1+b
	phi0.max <- max(seq.phi0)
	upper.y <- qnbinom(1-prec/2, a.post, min(b.post(seq.phi0)/(b.post(seq.phi0)+T2)))
	upper.x <- qnbinom(1-prec/2, a.post, min(b.post(seq.phi0)/(b.post(seq.phi0)+seq.phi0*S2))) 
 	bounds.array <- array(NA, dim=c(upper.x+1,upper.y+1,2))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			bounds <- Credibility.Interval(x1+x,y1+y,S1+S2,T1+T2,a,b,c,d,conf=K,interval)
			bounds.array[x+1, y+1,1] <- bounds[1]
			bounds.array[x+1, y+1,2] <- bounds[2]
		}
	}
	l <- length(seq.phi0)
	Reject <- rep(0,l)
	for(j in 1:l){
		phi0 <- seq.phi0[j]
		upper.y <- qnbinom(1-prec/2, a.post, b.post(phi0)/(b.post(phi0)+T2))
		upper.x <- qnbinom(1-prec/2, a.post, b.post(phi0)/(b.post(phi0)+phi0*S2))
		for(x2 in 0:upper.x){
			for(y2 in 0:upper.y){
				bounds <- bounds.array[x2+1, y2+1,]
				if(phi.star<bounds[1] | phi.star>bounds[2]){
					prob <- dpost.xy.mid.phi(phi0,x2,y2,S2,T2,a,b,x1,y1,S1,T1)
					Reject[j] <- Reject[j] + prob
				}
			}
		}
	}
list(phi0=seq.phi0, power=round(Reject, digits))
}























########################################################
################ intrinsic discrepancy #################
########################################################
rho <- function(phi, phi0, S, T){
      coef <- S/T
      gamma <- phi*coef
      gamma0 <- phi0*coef
      G <- gamma * log(gamma/gamma0) + log((gamma0+1)/(gamma+1))*(gamma+1)
      H <- gamma+1 - exp(gamma0/(gamma0+1)*log(gamma/gamma0))*(gamma0+1)
      result <- pmin(G,H)
      result
}

intrinsic.discrepancy <- function(phi, phi0, mu, S, T){
      result <- mu * T * rho(phi, phi0, S, T)
      result
}

#############################################################
##################### posterior cost ######################## # faire aussi avec simulations 
#############################################################
posterior.cost <- function(phi0, x, y,  S, T, a=0.5, b=0, c=0.5, d=0, subd=10000, tol=1e-6){
      post.c <- x+c
      post.d <- y+a+d
      post.a <- x+y+a
      lambda <- (T+b)/S
      K <- post.a*post.d/(post.c+post.d)*T/(T+b)
      integrande <- function(u){
            phi <- lambda * u/(1-u)
            rho(phi, phi0, S, T)*dbeta(u, post.c, post.d+1)
      }
	i <- -3
	old.value <-0
	value <- Inf
	while(abs(value-old.value)>tol){
		old.value <- value
		i <- i-1
		M <- qbeta(1-10^i,post.c, post.d+1)
	      value <- integrate(integrande, 0, M, subdivisions=subd)$value
	}
K*value
}

#############################################################
################### set posterior cost ######################
#######  alternative="less" for  H1: phi0 < phi.star  #######
######  alternative="gretaer" for  H1: phi0 > phi.star  #####
#############################################################
set.posterior.cost <- function(phi.star, alternative, x, y, S, T, a=0.5, b=0, c=0.5, d=0, subd=10000, tol=1e-6){
      post.c <- x+c
      post.d <- y+a+d
      post.a <- x+y+a
      lambda <- (T+b)/S
      K <- post.a*post.d/(post.c+post.d)*T/(T+b)
      integrande <- function(u){
            phi <- lambda * u/(1-u)
            rho(phi, phi.star, S, T)*dbeta(u, post.c, post.d+1)
      }
	psi.star <- phi.star/lambda
	u.star <- psi.star/(1+psi.star)
	bounds <- switch(alternative, less=c(0,u.star), greater=c(u.star, 1))
	value <- integrate(integrande, bounds[1], bounds[2], subdivisions=subd, rel.tol=tol)$value
	K*value
}

###########################################################
### significance curve #################
##########################################
intrinsic.scurve <- function(K, seq.mu0, phi.star, alternative, S, T, a=0.5, b=0, c=0.5, d=0, prec=1e-5, digits=4){
        mu0.max <- max(seq.mu0)
        lambda.max <- phi.star*mu0.max
        upper.x <- qpois(1-prec/2, lambda.max*S)
        upper.y <- qpois(1-prec/2, mu0.max*T)
 	spcost.array <- array(NA, dim=c(upper.x+1, upper.y+1))
	for(x in 0:upper.x){
		for(y in 0:upper.y){
			spcost.array[x+1,y+1] <- set.posterior.cost(phi.star, alternative, x, y, S, T, a, b, c, d)
		}
	}
        l <- length(seq.mu0)
        Reject <- rep(0,l)
        names(Reject) <- paste("mu0=", seq.mu0, sep="")
        phi0 <- phi.star
        for(j in 1:l){
                mu0 <- seq.mu0[j]
                lambda0 <- phi0*mu0
                upper.x <- qpois(1-prec/2, lambda0*S)
                for(x in 0:upper.x){
                        for(y in 0:upper.y){
					test <- spcost.array[x+1, y+1]
					if(test > K){
                                        Reject[j] <- Reject[j] + dpois(x, S*lambda0)*dpois(y, T*mu0)
                                }
                        }
                }
        }
list(mu0=seq.mu0, alpha=round(Reject, digits))
}

##################################################################
##################### intrinsic estimator ########################
##################################################################
intrinsic.estimate<-function(x, y, S, T, a, b, c, d, subd=10000, tol = 1e-08){
      post.cost <- function(u0){
            phi0 <- u0/(1-u0)
            posterior.cost(phi0, x, y, S, T, a, b, c, d, subd)
      }
      optimize <- optimize(post.cost, c(0, 1), tol=tol)
      u0.min <- optimize$minimum
      estimator <- u0.min/(1-u0.min)
      loss <- optimize$objective
      out <- list(estimator, loss)
      names(out) <- c("estimator", "loss")
out
}

##################################################################
##################### intrinsic credible interval ################
##################################################################
bounds.intrinsic<-function(x, y, S, T, a, b, c, d, conf, parameter="phi", subd=10000, tol = 1e-08){
      post.cost <- function(phi0){
            posterior.cost(phi0, x, y, S, T, a, b, c, d, subd)
      }
      post.icdf <- function(p){
            qpost.phi(p, a, b, c, d, S, T, x, y)
      }
      conf <- min(conf, 1 - conf)
      f <- function(p, post.icdf, conf){
            u.phi <- post.icdf(1 - conf + p) 
            l.phi <- post.icdf(p)
            (post.cost(u.phi)-post.cost(l.phi))^2
      }
      minimize <- optimize(f, c(0, conf), post.icdf = post.icdf, 
          conf = conf, tol=tol)$minimum
	out <- switch(parameter, 
			phi=c(post.icdf(minimize), post.icdf(1 - conf + minimize)),
			VE = sort(1-c(post.icdf(minimize), post.icdf(1 - conf + minimize))))
out
}



########################### examples ##########################################
a<-0.5; b<-0; c<-1/2; d<-0; S<-100; T<-S; x<-0; y<-20
intrinsic.estimate(x, y, S, T, a, b, c, d)
bounds<-bounds.intrinsic(x, y, S, T, a, b, c, d, conf=0.95); bounds
ppost.phi(bounds[2], a, b, c, d, S, T,  x, y)- ppost.phi(bounds[1], a, b, c, d, S, T, x, y)
################################################################################






