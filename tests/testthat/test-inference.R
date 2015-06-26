context("Inference")

test_that("Credibility intervals", {
  a <- 0.5; b <- 0; c <- 0.5; d <- 0; S <- 10; T <- 5; x <- 10; y <- 10
  intervals <- brr_intervals(x, y, S, T, a, b, c, d, intervals=c("left","equi"))
  expect_is(intervals, "list")
  expect_identical(unname(intervals$equi[1]), qpost_phi(.025, a, b, c, d, S, T, x, y))
})