context("Intrinsic loss")

test_that("Figure 10 JSPI", {
  loss <- intrinsic_phi0(phi0=0.5, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)
  expect_less_than(loss, log(1.713))
  expect_less_than(log(1.7129), loss)
  phi0 <- c(0.12, 0.1202, .994, .995)
  loss <- sapply(phi0, function(phi0) intrinsic_phi0(phi0=phi0, x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)) 
  expect_true(loss[3] < log(10) & loss[4] > log(10))
  expect_true(loss[1] > log(10) & loss[2] < log(10))
  expect_equal(intrinsic_estimate(x=6, y=15, S=1, T=1, a = 0.5, b = 0, c = 0.5, d = 0)$estimate, .4090697, tol=1e-7)
})

test_that("Figure 15 Master", {
  x=10; y=20; S=1; T=1; a=0.5; b=0; c=0.5; d=0
  bounds <- intrinsic_bounds(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d, conf=90/100)
  loss1 <- intrinsic_phi0(bounds[1], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
  loss2 <- intrinsic_phi0(bounds[2], x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
  expect_equal(loss1, loss2)
  expect_equal(loss1, 1.680312, tolerance=1e-6)
  estimate <- intrinsic_estimate(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)
  expect_equal(estimate$estimate, 0.5055056, tolerance=1e-6)
  expect_equal(estimate$loss, 0.4693617, tolerance=1e-6)
  mode <- spost_phi(x=x, y=y, S=S, T=T, a=a, b=b, c=c, d=d)$mode
  expect_less_than(mode, estimate$estimate)
})

