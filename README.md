## brr package for R

**Bayesian inference on the ratio of two Poisson rates.**

### What does it do ? ###

Suppose you have two counts of events and, assuming each count follows a Poisson distribution with an unknown incidence rate, you are interested in the ratio of the two rates (or *relative risk*).  The `brr` package allows to perform the Bayesian analysis of the relative risk using the natural semi-conjugate family of prior distributions, with a default non-informative prior (see references).


### Install ###

You can install:

- the latest released version from CRAN with 

```r
install.packages("brr")
```

- the latest development version from `github` using the `devtools` package:

```r
devtools::install_github('stla/brr', build_vignettes=TRUE)
```

### Basic usage ###

Create a `brr` object with the `Brr` function to set the prior parameters `a`, `b`, `c`, `d`, the two Poisson counts `x` and `y` and the samples sizes (times at risk) `S` and `T` in the two groups. Simply do not set the prior parameters to use the non-informative prior:

```r
model <- Brr(x=2, S=17877, y=9, T=16674)
``` 

Plot the posterior distribution of the rate ratio `phi`:

```r
plot(model, dpost(phi))
```

Get credibility intervals about `phi`:

```r
confint(model)
```

Get the posterior probability that `phi>1`:

```r
ppost(model, "phi", 1, lower.tail=FALSE)
```

Update the `brr` object to include new sample sizes and get a summary of the posterior predictive distribution of `x`:

```r
model <- model(Snew=10000, Tnew=10000)
spost(model, "x", output="pandoc")
```

### To learn more ###

Look at the vignettes:

```r
browseVignettes(package = "brr")
```

### Find a bug ? Suggestion for improvment ?

Please report at https://github.com/stla/brr/issues

### References ###

S. Laurent, C. Legrand: *A Bayesian framework for the ratio of two Poisson rates in the context of vaccine efficacy trials.* ESAIM, Probability \& Statistics 16 (2012), 375--398.

S. Laurent: *Some Poisson mixtures  distributions with a hyperscale parameter.* Brazilian Journal of Probability and Statistics 26 (2012), 265--278.

S. Laurent: *Intrinsic Bayesian inference on a Poisson rate and on the ratio of two Poisson rates.* Journal of Statistical Planning and Inference 142 (2012), 2656--2671.
