## brr package for R
Bayesian inference on the ratio of two Poisson rates.


### Install ###

The `brr` package is still under development. You can install it from `github` using the `devtools` package. 

```r
devtools::install_github('brr', 'stla')
```

### Basic usage ###

Create a `brr` object with the `Brr` function to set the prior parameters `a`, `b`, `c`, `d`, the two Poisson counts `x` and `y` and the samples sizes (times at risk) `S` and `T`. Simply do not set the prior parameters to use the non-informative prior:

```r
model <- Brr(x=2, S=17877, y=9, T=16674)
``` 

Plot the posterior distribution of the rate ratio `phi`:

```r
plot(model, dpost(phi))
```

Get a credibility interval about `phi`:

```r
confint(model)
```

Get the probability that `phi>1`:

```r
ppost(model, "phi", 1, lower.tail=FALSE)
```

Update the `brr` object to include new sample sizes and get a summary of the posterior predictive distribution of `x`:

```r
model <- model(Snew=10000, Tnew=1000)
spost(model, "x", type="pandoc")
```
