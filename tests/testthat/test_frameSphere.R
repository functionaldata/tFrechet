library(testthat)

test_that("Combining with x yields an orthogonal matrix if the input x is a unit vector", {
  d <- 5
  set.seed(1)
  x <- rnorm(d)
  x <- x / sqrt(sum(x^2))
  A <- cbind(x,frameSphere(x))
  expect_lt(sum(abs( A %*% t(A) - diag(d) )), 1e-10)
})
