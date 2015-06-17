context("Quantiles")

test_that("Quantiles PGIB", {
  p <- 0.5
  q <- qPGIB(p, 10, 2, 12, 2.1) 
  expect_true(pPGIB(q, 10, 2, 12, 2.1) > p && pPGIB(q-1, 10, 2, 12, 2.1) < p)
  
})