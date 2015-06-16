context("Distributions")

test_that("dPGIB returns positive values", {
  expect_that(all(dPGIB(0:2,3,4,2,2.5)>0), is_true())
})

test_that("Check moment GB2", {
  a <- 3; c <- 4; d <- 5; tau <- 10
  expect_that((moment_GB2(1, a, c, d, tau) - summary_gamma(a, 1)$mean*summary_beta2(d, c, 1/tau)$mean) < .Machine$double.eps, 
              is_true())
})

