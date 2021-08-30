library(testthat)

test_that("works for location differences between the populations", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  delta <- 1
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 0.5)
  mu2 <- rnorm(n2, mean = delta, sd = 0.5)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- qSup
  group <- c(rep(1, n1), rep(2, n2))
  res <- DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile", boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})

test_that("works for scale differences between the populations", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  r <- 0.5
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  mu1 <- rnorm(n1, mean = 0, sd = 0.2)
  mu2 <- rnorm(n2, mean = 0, sd = 0.2*r)
  Y1 <- lapply(1:n1, function(i) {
    qnorm(qSup, mu1[i], sd = 1)
  })
  Y2 <- lapply(1:n2, function(i) {
    qnorm(qSup, mu2[i], sd = 1)
  })
  Ly <- c(Y1, Y2)
  Lx <- qSup
  group <- c(rep(1, n1), rep(2, n2))
  res <- DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile", boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})

test_that("works if the two populations are the same", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
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
  group <- c(rep(1, n1), rep(2, n2))
  res <- DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile", boot = TRUE))
  expect_equal(res$pvalAsy > .05 & res$pvalBoot > .05, TRUE)
})

test_that("Ly and Lx should have the same length", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
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
  group <- c(rep(1, n1), rep(2, n2))
  expect_error(DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile")), 
               "Ly and Lx should have the same length")
})

test_that("Ly and group should have the same length", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
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
  group <- c(rep(1, n1), rep(2, n2), 1)
  expect_error(DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile")), 
               "Ly and group should have the same length")
})

test_that("each vector in Ly and its corresponding vector in Lx should have the same length", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
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
  group <- c(rep(1, n1), rep(2, n2))
  expect_error(DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile")), 
               "each vector in Ly and its corresponding vector in Lx should have the same length")
})
