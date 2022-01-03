library(testthat)

test_that("A few examples", {
  expect_true(sum(abs( c(1, 0, pi/4) - car2pol(c(1,0,1)/sqrt(2)) )) < 1e-10)
  expect_true(sum(abs( c(1, pi, 0) - car2pol(c(-1,0,0)) )) < 1e-10)
})
