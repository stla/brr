context("Test miscellanous")

test_that("dd_moment", {
  set.seed(666)
  tol <- 1e-10
  # Poisson 
  lambda <- rnorm(3, 5)
  m <- sapply(lambda, function(lambda) dd_moment(dpois, lambda=lambda, accuracy=tol))
  expect_equal(m, lambda, tolerance=tol)
  m2 <- sapply(lambda, function(lambda) dd_moment(dpois, k=2, lambda=lambda, accuracy=tol))
  expect_equal(m2-m^2, lambda, tolerance=tol)
  # beta_nbinom
  accuracy <- 1e-12
  tol <- 1e-9
  a <- rnorm(3, 5)
  c <- 6; d <- 3
  m <- sapply(a, function(a) dd_moment(dbeta_nbinom, a=a, c=c, d=d, accuracy=accuracy))
  m2 <- sapply(a, function(a) dd_moment(dbeta_nbinom, k=2, a=a, c=c, d=d, accuracy=accuracy))
  m1 <- sapply(a, function(a) summary_beta_nbinom(a=a, c=c, d=d)$mean) 
  SD <- sapply(a, function(a) summary_beta_nbinom(a=a, c=c, d=d)$sd) 
  expect_equal(m, m1, tolerance=tol)
  expect_equal(m2-m^2, SD^2, tolerance=tol)
})