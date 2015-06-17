context("Intrinsic loss")

test_that("Figure 10", {
  loss <- intrinsic_phi0(phi0=0.5, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)
  expect_less_than(loss, log(10)/2)
  phi0 <- c(0.12, 0.1202, .994, .995)
  loss <- sapply(phi0, function(phi0) intrinsic_phi0(phi0=phi0, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)) 
  expect_true(loss[3] < log(10) & loss[4] > log(10))
  expect_true(loss[1] > log(10) & loss[2] < log(10))
  expect_equal(intrinsic_estimate(x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)$estimate, .4090697, tol=1e-7)
})