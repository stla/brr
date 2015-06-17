context("Brr")

test_that("Test Brr function", {
  model <- Brr(a=2, b=4)
  expect_that(
    identical(model(), list(a=2,b=4)), 
    is_true()
  )
  model <- model(a=3)
  expect_that(
    identical(model(), list(a=3,b=4)), 
    is_true()
  )
  model <- model(c=4)
  expect_that(
    identical(model()[c("a","b","c")], list(a=3,b=4,c=4)), 
    is_true()
  )
})

test_that("Test dprior function", {
  model <- Brr(a=2, b=4)
  expect_that(
    all(dprior(model, "mu", 1:3) == dprior_mu(mu=1:3, a=2, b=4)), 
    is_true()
  )
  model <- Brr(a=2, b=4, c=3, d=10, T=10)
  expect_that(
    all(dprior(model, "x", 1:3) == dprior_x(x=1:3, a=2, b=4, c=3, d=10, T=10)), 
    is_true()
  )
})

test_that("Test dpost function", {
  # mu
  model <- Brr(a=2, b=4)
  expect_error(
    dpost(model, "mu", 1)
  )
  model <- model(c=3, d=4, x=5, y=6, T=8)
  expect_that(
    all(dpost(model, "mu", 1:3) == dpost_mu(mu=1:3, a=2, b=4, c=3, d=4, x=5, y=6, T=8)), 
    is_true()
  )  
  ## the same :
  expect_equal(dpost(model, "mu", 1:3),  dpost_mu(mu=1:3, a=2, b=4, c=3, d=4, x=5, y=6, T=8))
  # phi
  expect_error(
    dpost(model, "phi", 1)
  )
  model <- model(S=10)
  expect_that(
    all(dpost(model, "phi", 1:3) == dpost_phi(phi=1:3, a=2, b=4, c=3, d=4, x=5, y=6, S=10, T=8)), 
    is_true()
  )  
#   expect_that(
#     all(dpost(model, "x", 1:3) == dprior_x(x2=1:3, a=2, b=4, c=3, d=10, T=10)), 
#     is_true()
#   )
})