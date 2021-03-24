library(testthat)

test_that("A few examples", {
  expect_true(sum(abs(pol2car(c(1, 0, pi/4)) - c(1,0,1)/sqrt(2))) < 1e-10)
  expect_true(sum(abs(pol2car(c(1, pi, 0)) - c(-1,0,0))) < 1e-10)
})
