context("Quantiles")

test_that("Quantiles PGIB and PGB2", {
  p <- 0.5
  q <- qPGIB(p, 10, 2, 12, 2.1) 
  expect_true(pPGIB(q, 10, 2, 12, 2.1) > p && pPGIB(q-1, 10, 2, 12, 2.1) < p)
  q <- qPGB2(p, 10, 2, 12, 2.1) 
  expect_true(pPGB2(q, 10, 2, 12, 2.1) > p && pPGB2(q-1, 10, 2, 12, 2.1) < p)
})