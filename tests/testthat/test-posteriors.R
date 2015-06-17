context("Posteriors")

test_that("Check stochastic monotonicity of mu", {
  a <- 3; b <- 7; c <- 4; d <- 5; S <- 10; T <- 13; x <- 1; y <- 0
  p1 <- integrate(function(mu) dpost_mu(mu,a,b,c,d,T,x,y), lower=0, upper=.1)
  p2 <- integrate(function(mu) dpost_mu(mu,a,b,c,d,T+1,x,y), lower=0, upper=.1)
  expect_true(p1$message=="OK" && p2$message=="OK")
  expect_less_than(p1$value, p2$value)
})