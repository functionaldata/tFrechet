library(testthat)

test_that("works for networks from the stochatic block model with different preferential matrices", {
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
  group <- c(rep(1, n1), rep(2, n2))
  res <- NetANOVA(Ly, group, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
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
  group <- c(rep(1, n1), rep(2, n2))
  res <- NetANOVA(Ly, group, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-3 & res$pvalBoot < 1e-3, TRUE)
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
  group <- c(rep(1, n1), rep(2, n2))
  res <- NetANOVA(Ly, group, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy > .05 & res$pvalBoot > .05, TRUE)
})

test_that("Ly and group should have the same length", {
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
  group <- c(rep(1, n1), rep(2, n2), 1)
  expect_error(NetANOVA(Ly, group), 
               "Ly and group should have the same length")
})

test_that("works for more than two groups", {
  set.seed(1)
  n1 <- 100
  n2 <- 100
  n3 <- 50
  gamma1 <- 2
  gamma2 <- 2.5
  gamma3 <- 3
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
  Y3 <- lapply(1:n3, function(i) {
    igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma3, 
                                               directed = FALSE), 
                             sparse = FALSE)
  })
  Ly <- c(Y1, Y2, Y3)
  group <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  res <- NetANOVA(Ly, group, optns = list(boot = TRUE))
  expect_equal(res$pvalAsy < 1e-5 & res$pvalBoot < 1e-5, TRUE)
})