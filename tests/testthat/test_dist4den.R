require(testthat)

test_that("Works with for density input", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$y <- dnorm(d1$x)
  d2 <- list(x = d1$x + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  dist <- dist4den(d1 = d1,d2 = d2)
  expect_true(abs(dist-1) < 1e-10)
})

test_that("Works with for quantile input", {
  qSup <- seq(0,1,length.out = 201)
  q1 <- qunif(qSup, 0,1)
  q2 <- qunif(qSup, 1,2)
  dist <- dist4den(list(x = qSup, y = q1), list(x = qSup, y = q2), fctn_type = "quantile")
  expect_true(abs(dist-1) < 1e-10)
})

test_that("Stops when missing one of d1, d2", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$y <- dnorm(d1$x)
  expect_error(dist4den(d1 = d1), "Requires the input of both d1 and d2")
})

test_that("Stops when fctn_type is mis-specified", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$y <- dnorm(d1$x)
  d2 <- list(x = d1$x + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2, fctn_type = "pdf"), "Unrecognized value of fctn_type")
})

test_that("Stops when d1 is not a list", {
  qSup <- seq(0,1,length.out = 201)
  d1 <- qunif(qSup, 0,1)
  d2 <- list(x = seq(-6,6,0.01) + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2),"d1 and d2 should be lists.")
})

test_that("Stops when d1 does not have elements x or y", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$z <- dnorm(d1$x)
  d2 <- list(x = d1$x + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2), "d1 should consist of two elements x and y.")
})

test_that("Stops when d1$x and d1$y have different lengths", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$y <- c(dnorm(d1$x),0)
  d2 <- list(x = d1$x + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2),"d1\\$x and d1\\$y should have the same length.")
})

test_that("Stops when d1 does not integrate to 1 with tolerance of 1e-05", {
  d1 <- list(x = seq(-3,3,0.01))
  d1$y <- dnorm(d1$x)
  d2 <- list(x = seq(-6,6,0.01) + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2),"d1 should be a density function")
})

test_that("Stops when d1$y is not all non-negative", {
  d1 <- list(x = seq(-6,6,0.01))
  d1$y <- dnorm(d1$x)
  d1$y[1] <- -1
  d2 <- list(x = seq(-6,6,0.01) + 1)
  d2$y <- dnorm(d2$x, mean = 1)
  expect_error(dist4den(d1 = d1,d2 = d2),"d1 should be a density function")
})

test_that("Stops when quantile functions given in d1 is not monotonic", {
  qSup <- seq(0,1,length.out = 201)
  q1 <- qunif(qSup, 0,1)
  q2 <- qunif(qSup, 1,2)
  q2[1:3] <- q2[3:1]
  expect_error(
    dist4den(list(x = qSup, y = q1), list(x = qSup, y = q2), fctn_type = "quantile"),
    "Quantile functions given in d1 and d2 are not monotonic"
  )
})
