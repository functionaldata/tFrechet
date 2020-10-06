library(testthat)

test_that("Works with for spherical data", {
  y1 <- c(0,0,1)
  y2 <- c(0,1,0)
  dist <- SpheGeoDist(y1,y2)
  expect_true(abs(dist-pi/2) < 1e-10)
})

test_that("Stops when y1 and y2 have different lengths", {
  y1 <- c(0,0,1,0)
  y2 <- c(0,1,0)
  expect_error(SpheGeoDist(y1,y2), "y1 and y2 should be of the same length.")
})

test_that("Stops when y1 or y2 are not unit vectors", {
  y1 <- c(0,0,1.1)
  y2 <- c(0,1,0)
  expect_error(SpheGeoDist(y1,y2), "y1 is not a unit vector.")
})

test_that("Stops when y1 or y2 are not numeric vectors", {
  y1 <- c("A","b",0)
  y2 <- c(0,1,0)
  expect_error(SpheGeoDist(y1,y2))
})
