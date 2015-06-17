context("Test miscellanous")

test_that("dd_moment", {
  set.seed(666)
  lambda <- rnorm(3, 5)
  tol=1e-10
  m <- sapply(lambda, function(lambda) dd_moment(dpois, lambda=lambda, accuracy=1e-10))
  expect_equal(m, lambda, tolerance=tol)
  m2 <- sapply(lambda, function(lambda) dd_moment(dpois, k=2, lambda=lambda, accuracy=1e-10))
  expect_equal(m2-m^2, lambda, tolerance=tol)    
})