context("Simulations")

test_that("Check PGB2 using simulations",{
  set.seed(666)
  a <- 2 ; c <- 5 ; d <- 30; tau <- 6
  sims <- rPGB2(1e6, a, c, d, tau)
  p1 <- ecdf(sims)(1:3)
  p2 <- pPGB2(1:3, a, c, d, tau)
  expect_equal(p1, p2, tolerance=1e-3)
})

test_that("Check PGIB using simulations",{
  set.seed(666)
  a <- 2 ; alpha <- 5 ; beta <- 30; rho <- 0.6
  sims <- rPGIB(1e6, a, alpha, beta, rho)
  p1 <- ecdf(sims)(1:3)
  p2 <- pPGIB(1:3, a,  alpha, beta, rho)
  expect_equal(p1, p2, tolerance=1e-3, scale=1)
})

