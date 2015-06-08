library(brr)


# posterior mu  -------------------------------------------------

a <- 0.5; b <- 0; c <- 0.5; d <- 0
x <- 10; y <- 10
S <- 1; T <- 1

# centered around 10
seq(0,20,len=100) %>% {plot(., dpost_mu(.,a,b,c,d,T,x,y), type="l")}

# decreases as T increases 
seq(0,20,len=100) %>% {plot(., dpost_mu(.,a,b,c,d,T=2,x,y), type="l")}

# highly concentrated prior around 2 (crash for a=150)
seq(0,5,len=100) %>% {plot(., dpost_mu(.,a=140,b=70,c,d,T,x,y), type="l")}
abline(v=2)
seq(0,5,len=100) %>% {lines(., dpost_mu(.,a=140,b=70,c,d,T,x,y=2), col="red")}


# posterior lambda --------------------------------------------------------

# centered around 10
seq(0,20,len=100) %>% {plot(., dpost_lambda(.,a,c,d,S,x,y), type="l")}

# decreases as S increases
seq(0,20,len=100) %>% {plot(., dpost_lambda(.,a,c,d,S=2,x,y), type="l")}

# highly concentrated prior around 1 x 2
seq(0,5,len=100) %>% {plot(., dprior_phi(., b=50, c=7, d=80,S=4,T), type="l")}
seq(0,5,len=100) %>% {plot(., dprior_lambda(.,a=100,b=50,c=7,d=80,S=4,T), type="l")}
abline(v=2)

# crash de dprior_lambda :
dprior_lambda(3.3,a=100,b=50,c=4,d=120,S,T)
dprior_lambda(3.4,a=100,b=50,c=4,d=120,S,T)
a=100;b=50;c=4;d=120
lambda=3.4
library(gsl)
log(hyperg_U(a+d,a-c+1,b*S/(b+T)*lambda))
                                                                
seq(0,5,len=100) %>% {plot(., dpost_lambda(.,a=100,c=7,d=80,S=4,x=8,y=2), type="l")}
abline(v=2)



# posterior predictive x --------------------------------------------------


