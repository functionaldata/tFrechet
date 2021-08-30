library(testthat)

test_that("works for graph Laplacian matrices", {
  set.seed(1)
  n <- 100
  gamma <- 2
  Ly <- lapply(1:n, function(i) {
    igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma, 
                                               directed = FALSE), 
                             sparse = FALSE)
  })
  mup <- rowMeans(matrix(unlist(Ly), ncol = n))
  Vp <- mean(sapply(Ly, function(Lyi) sum((Lyi - mup)^2)))
  res <- NetFVar(Ly)
  expect_equal((res$NetFVar-Vp)^2 < 1e-3, TRUE)
})

test_that("works for covariance matrices", {
  set.seed(1)
  n <- 100
  U <- pracma::randortho(10)
  Ly <- lapply(1:n, function(i) {
    U %*% diag(rexp(10, (1:10)/2)) %*% t(U)
  })
  mup <- rowMeans(matrix(unlist(Ly), ncol = n))
  Vp <- mean(sapply(Ly, function(Lyi) sum((Lyi - mup)^2)))
  res <- NetFVar(Ly)
  expect_equal((res$NetFVar-Vp)^2 < 1e-3, TRUE)
})

test_that("works for correlation matrices", {
  set.seed(1)
  n <- 100
  m <- 10
  d <- m * (m - 1) / 2
  alpha <- 2
  beta <- 2
  Ly <- lapply(1:n, function(i) {
    yVec <- rbeta(d, shape1 = alpha, shape2 = beta)
    y <- matrix(0, nrow = m, ncol = m)
    y[lower.tri(y)] <- yVec
    y <- y + t(y)
    diag(y) <- 1
    y
  })
  mup <- rowMeans(matrix(unlist(Ly), ncol = n))
  Vp <- mean(sapply(Ly, function(Lyi) sum((Lyi - mup)^2)))
  res <- NetFVar(Ly)
  expect_equal((res$NetFVar-Vp)^2 < 1e-3, TRUE)
})

test_that("each matrix in Ly should be of the same dimension", {
  set.seed(1)
  n <- 100
  U <- pracma::randortho(10)
  Ly <- lapply(1:n, function(i) {
    U %*% diag(rexp(10, (1:10)/2)) %*% t(U)
  })
  Ly[[1]] <- matrix(0, nrow = 20, ncol = 20)
  expect_error(NetFVar(Ly), "each matrix in Ly should be of the same dimension")
})

test_that("each matrix in Ly should be a square matrix", {
  set.seed(1)
  n <- 100
  U <- pracma::randortho(10)
  Ly <- lapply(1:n, function(i) {
    y <- U %*% diag(rexp(10, (1:10)/2)) %*% t(U)
    matrix(y, nrow = 5)
  })
  expect_error(NetFVar(Ly), "each matrix in Ly should be a square matrix")
})
