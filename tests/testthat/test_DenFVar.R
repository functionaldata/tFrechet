library(testthat)

test_that("works for density samples", {
  set.seed(1)
  n <- 100
  alpha <- runif(n, min = 1, max = 2)
  beta <- runif(n, min = 1, max = 2)
  Lx <- seq(0, 1, 0.01)
  Ly <- lapply(1:n, function(i) dbeta(Lx, shape1 = alpha[i], shape2 = beta[i]))
  # normalize each Ly[[i]] to ensure that it integrates to 1
  Ly <- t(fdadensity::normaliseDensities(matrix(unlist(Ly), nrow = n, byrow = TRUE), Lx))
  Ly <- split(Ly, col(Ly))
  res <- DenFVar(din = Ly, supin = Lx)
  expect_equal((res$DenFVar - 0.0044)^2 < 1e-6, TRUE)
})

test_that("works for quantile samples", {
  set.seed(1)
  n <- 100
  mu <- rnorm(n, mean = 0, sd = 0.5)
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  Ly <- lapply(1:n, function(i) qnorm(qSup, mu[i], sd = 1))
  Lx <- qSup
  res <- DenFVar(qin = Ly, supin = Lx)
  DenFMean <- rowMeans(matrix(unlist(Ly), nrow = length(Ly[[1]]), ncol = n))
  DenFVar <- mean(sapply(Ly, function(Lyi) {
    pracma::trapz(qSup, (Lyi - DenFMean)^2)
  }))
  expect_equal((res$DenFVar - DenFVar)^2 < 1e-6, TRUE)
})

test_that("Works with noisy yin input", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  sd <- 0.01
  yin <- lapply(xin, function(x) {
    rnorm(1000, x, sd)
  })
  qSup <- qbeta((1:99)/100,1/2,1/2)
  res <- DenFVar(yin = yin, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  DenFMean <- colMeans(qtrue)
  DenFVar <- mean(apply(qtrue, 1, function(qtruei) pracma::trapz(qSup, (qtruei - DenFMean)^2)))
  expect_true(res$DenFVar - DenFVar < 1e-4)
})
