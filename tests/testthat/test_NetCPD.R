library(testthat)
# note: the following tests is very time-consuming 

test_that("works for networks from the stochastic block model with different preferential matrices", {
  set.seed(1)
  n1 <- 30
  n2 <- 30
  C1 <- cbind(c(0.8, 0.3, 0.2), c(0.3, 0.6, 0.1), c(0.2, 0.1, 0.75)) # 3 groups in sbm
  C2 <- matrix(c(0.9, 0.2, 0.2, 0.6), 2, 2) # 2 groups in sbm
  Y1 <- lapply(1:n1, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C1, 
                                                block.sizes = c(5, 5, 10)), 
                             sparse = FALSE)
  })
  Y2 <- lapply(1:n2, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C2, 
                                                block.sizes = c(10, 10)), 
                             sparse = FALSE)
  })
  Ly <- c(Y1, Y2)
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy < 1e-5, TRUE)
})

test_that("works for scale-free networks from the Barabasi-Albert model with different powers of the preferential attachment", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  gamma1 <- 2
  gamma2 <- 3
  Y1 <- lapply(1:n1, function(i) {
    igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma1, 
                                               directed = FALSE), 
                             sparse = FALSE)
  })
  Y2 <- lapply(1:n2, function(i) {
    igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma2, 
                                               directed = FALSE), 
                             sparse = FALSE)
  })
  Ly <- c(Y1, Y2)
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy < 1e-3, TRUE)
})

test_that("works if the two populations are the same", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  C1 <- cbind(c(0.8, 0.3, 0.2), c(0.3, 0.6, 0.1), c(0.2, 0.1, 0.75)) # 3 groups in sbm
  C2 <- C1
  Y1 <- lapply(1:n1, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C1, 
                                                block.sizes = c(5, 5, 10)), 
                             sparse = FALSE)
  })
  Y2 <- lapply(1:n2, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C2, 
                                                block.sizes = c(5, 5, 10)), 
                             sparse = FALSE)
  })
  Ly <- c(Y1, Y2)
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy > .05, TRUE)
})

test_that("works for covariance matrices", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  U <- pracma::randortho(10)
  Y1 <- lapply(1:n1, function(i) {
    U %*% diag(rexp(10, (1:10)/2)) %*% t(U)
  })
  Y2 <- lapply(1:n2, function(i) {
    U %*% diag(rexp(10, 1:10)) %*% t(U)
  })
  Ly <- c(Y1, Y2)
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy < 1e-1, TRUE)
})

test_that("works for correlation matrices", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  m <- 10
  d <- m * (m - 1) / 2
  alpha1 <- 1
  beta1 <- 1
  alpha2 <- 2
  beta2 <- 2
  Y1 <- lapply(1:n1, function(i) {
    yVec <- rbeta(d, shape1 = alpha1, shape2 = beta1)
    y <- matrix(0, nrow = m, ncol = m)
    y[lower.tri(y)] <- yVec
    y <- y + t(y)
    diag(y) <- 1
    y
  })
  Y2 <- lapply(1:n2, function(i) {
    yVec <- rbeta(d, shape1 = alpha2, shape2 = beta2)
    y <- matrix(0, nrow = m, ncol = m)
    y[lower.tri(y)] <- yVec
    y <- y + t(y)
    diag(y) <- 1
    y
  })
  Ly <- c(Y1, Y2)
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy < 1e-5, TRUE)
})

test_that("works for array input", {
  set.seed(1)
  n1 <- 30
  n2 <- 30
  C1 <- cbind(c(0.8, 0.3, 0.2), c(0.3, 0.6, 0.1), c(0.2, 0.1, 0.75)) # 3 groups in sbm
  C2 <- matrix(c(0.9, 0.2, 0.2, 0.6), 2, 2) # 2 groups in sbm
  Y1 <- lapply(1:n1, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C1, 
                                                block.sizes = c(5, 5, 10)), 
                             sparse = FALSE)
  })
  Y2 <- lapply(1:n2, function(i) {
    igraph::laplacian_matrix(igraph::sample_sbm(n = 20, pref.matrix = C2, 
                                                block.sizes = c(10, 10)), 
                             sparse = FALSE)
  })
  Ly <- array(unlist(c(Y1, Y2)), c(20, 20, n1 + n2))
  res <- NetCPD(Ly)
  expect_equal(res$pvalAsy < 1e-5, TRUE)
})
