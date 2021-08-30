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
  res <- DenFVar(Ly, Lx)
  expect_equal((res$DenFVar - 0.0044)^2 < 1e-6, TRUE)
})

test_that("works for quantile samples", {
  set.seed(1)
  n <- 100
  mu <- rnorm(n, mean = 0, sd = 0.5)
  qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
  Ly <- lapply(1:n, function(i) qnorm(qSup, mu[i], sd = 1))
  Lx <- qSup
  res <- DenFVar(Ly, Lx, optns = list(fctn_type = "quantile"))
  DenFMean <- rowMeans(matrix(unlist(Ly), nrow = length(Ly[[1]]), ncol = n))
  DenFVar <- mean(sapply(Ly, function(Lyi) {
    pracma::trapz(qSup, (Lyi - DenFMean)^2)
  }))
  expect_equal((res$DenFVar - DenFVar)^2 < 1e-6, TRUE)
})
