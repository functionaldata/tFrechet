library(testthat)

test_that("works for location differences", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  delta <- 0.75
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = delta, sd = 0.5)
  mu2 <- rnorm(n2, mean = 0, sd = 0.5)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- qSup
  res <-DenCPD(qin = Ly, supin = Lx, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})

test_that("works for scale differences", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  r <- 0.5
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 1)
  mu2 <- rnorm(n2, mean = 0, sd = r)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- qSup
  res <- DenCPD(qin = Ly, supin = Lx, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})

test_that("works if the two populations are the same", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 0.5)
  mu2 <- rnorm(n2, mean = 0, sd = 0.5)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- qSup
  res <- DenCPD(qin = Ly, supin = Lx, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy > .05 & res$pvalBoot > .05, TRUE)
})

test_that("works for density samples", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  alpha1 <- runif(n1, min = 1, max = 2)
  beta1 <- runif(n1, min = 1, max = 2)
  alpha2 <- runif(n2, min = 2, max = 3)
  beta2 <- runif(n2, min = 2, max = 3)
  Lx <- seq(0, 1, 0.01)
  Y1 <- lapply(1:n1, function(i) dbeta(Lx, shape1 = alpha1[i], 
                                       shape2 = beta1[i]))
  Y2 <- lapply(1:n2, function(i) dbeta(Lx, shape1 = alpha2[i], 
                                       shape2 = beta2[i]))
  Ly <- c(Y1, Y2)
  # normalize each Ly[[i]] to ensure that it integrates to 1
  Ly <- t(fdadensity::normaliseDensities(matrix(unlist(Ly), nrow = n1 + n2, 
                                                byrow = TRUE), Lx))
  Ly <- split(Ly, col(Ly))
  res <- DenCPD(din = Ly, supin = Lx, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})

test_that("the number of support grids in supin is not equal to the number of observed distributions in qin", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 0.5)
  mu2 <- rnorm(n2, mean = 0, sd = 0.5)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- rep(list(qSup), n1 + n2 + 1)
  expect_error(DenCPD(qin = Ly, supin = Lx), 
               "the number of support grids in supin is not equal to the number of observed distributions in qin")
})

test_that("the number of support points must be equal to the number of observations for each quantile function in qin", {
  set.seed(1)
  n1 <- 100
  n2 <- 200
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 0.5)
  mu2 <- rnorm(n2, mean = 0, sd = 0.5)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- rep(list(qSup), n1 + n2)
  Lx[[1]] <- c(Lx[[1]], 1)
  expect_error(DenCPD(qin = Ly, supin = Lx), 
               "the number of support points must be equal to the number of observations for each quantile function in qin")
})
