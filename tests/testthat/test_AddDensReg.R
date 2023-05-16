library(testthat)


test_that("Gives outputs invariant to shifting dSup", {
  
  library(MASS)
  
  # additive component functions
  g1 <- function (u, x1) sin(2*pi*u) * (2*x1 - 1)
  g2 <- function (u, x2) sin(2*pi*u) * sin(2*pi*x2)
  
  g <- function (u, x) g1(u, x[1]) + g2(u, x[2])
  
  # generating random samples from conditional quantile functions
  GenLqdNoise <- function (u, e) e[1]*sin(pi*u) + e[2]*sin(2*pi*u) 
  GenQdensResp <- function (u, x, e) exp(g(u, x) + GenLqdNoise(u, e))
  
  GenQuantileResp <- function (u, x, e) {
    
    tmp1 <- integrate(GenQdensResp, lower = 0, upper = u, x = x, e = e)$value
    tmp2 <- integrate(GenQdensResp, lower = 0, upper = 1, x = x, e = e)$value
    
    return (tmp1 / tmp2)
  }
  
  set.seed(999)
  n <- 150
  N <- 250
  
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
  X <- pnorm(mvrnorm(n, rep(0, 2), Sigma))
  
  Ly <- Ly_shift <- list()
  for (i in 1:n) {
    U_i <- runif(N)
    E_i <- c(rnorm(1, 0, 0.1), rnorm(1, 0, 0.05))
    Ly[[i]] <- sapply(1:N, function (l) GenQuantileResp(U_i[l], X[i,], E_i))
    Ly_shift[[i]] <- 2 * Ly[[i]] - 1 
  }
  
  M <- 51
  x1 <- x2 <- seq(0, 1, length.out = M)
  x <- cbind(x1, x2)
  
  hu <- 0.05
  hx <- c(0.075, 0.075)
  dSup <- c(0, 1)
  dSup_shift <- c(-1, 1)
  
  # estimating the functional additive model
  estAddDensReg <- frechet::AddDensReg(Ly = Ly, X = X, x = x, hu = hu, hx = hx, dSup = dSup)
  estAddDensReg_shift <- frechet::AddDensReg(Ly = Ly_shift, X = X, x = x, hu = hu, hx = hx, dSup = dSup_shift)
  
  # LQD component function estimates
  g0Sbf <- estAddDensReg$lqdSbfMean
  gjSbf <- estAddDensReg$LlqdSbfComp
  
  g0Sbf_shift <- estAddDensReg_shift$lqdSbfMean
  gjSbf_shift <- estAddDensReg_shift$LlqdSbfComp
  
  err <- max(abs(g0Sbf - g0Sbf_shift), abs(gjSbf[[1]] - gjSbf_shift[[1]]), abs(gjSbf[[2]] - gjSbf_shift[[2]]))
  
  expect_true(err < 1e-10)
})



test_that("Satisfies the asymptotic rate of convergence", {
  
  library(MASS)
  
  # additive component functions
  g1 <- function (u, x1) sin(2*pi*u) * (2*x1 - 1)
  g2 <- function (u, x2) sin(2*pi*u) * sin(2*pi*x2)
  
  g <- function (u, x) g1(u, x[1]) + g2(u, x[2])
  
  # generating random samples from conditional quantile functions
  GenLqdNoise <- function (u, e) e[1]*sin(pi*u) + e[2]*sin(2*pi*u) 
  GenQdensResp <- function (u, x, e) exp(g(u, x) + GenLqdNoise(u, e))
  
  GenQuantileResp <- function (u, x, e) {
    
    tmp1 <- integrate(GenQdensResp, lower = 0, upper = u, x = x, e = e)$value
    tmp2 <- integrate(GenQdensResp, lower = 0, upper = 1, x = x, e = e)$value
    
    return (tmp1 / tmp2)
  }
  
  set.seed(999)
  n <- 150
  N <- 250
  
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
  X <- pnorm(mvrnorm(n, rep(0, 2), Sigma))
  
  Ly <- list()
  for (i in 1:n) {
    U_i <- runif(N)
    E_i <- c(rnorm(1, 0, 0.1), rnorm(1, 0, 0.05))
    Ly[[i]] <- sapply(1:N, function (l) GenQuantileResp(U_i[l], X[i,], E_i))
  }
  
  M <- 51
  x1 <- x2 <- seq(0, 1, length.out = M)
  x <- cbind(x1, x2)
  
  hu <- 0.05
  hx <- c(0.075, 0.075)
  dSup <- c(0, 1)
  
  # estimating the functional additive model
  estAddDensReg <- frechet::AddDensReg(Ly = Ly, X = X, x = x, hu = hu, hx = hx, dSup = dSup)
  
  # true LQD component functions
  g1Eval <- g2Eval <- matrix(nrow = length(estAddDensReg$lqdGrid), ncol = M)
  for (l in seq(estAddDensReg$lqdGrid)) {
    for (m in seq(M)) {
      g1Eval[l,m] <- g1(estAddDensReg$lqdGrid[l], x1[m])
      g2Eval[l,m] <- g2(estAddDensReg$lqdGrid[l], x2[m])
    }
  }
  
  # LQD component function estimates
  gjSbf <- estAddDensReg$LlqdSbfComp
  
  err1 <- max(abs(g1Eval - gjSbf[[1]]))
  err2 <- max(abs(g2Eval - gjSbf[[2]]))
  
  densGrid <- estAddDensReg$densGrid
  bw_dens <- lapply(1:n,
                  function (i) {
                    bw_i <- frechet::CreateDensity(y = Ly[[i]], 
                                                   optns = list(outputGrid = densGrid))$bw
                    return (bw_i)
                  }
  )
  bw_err <- unlist(bw_dens)
    
  crit1 <- max(bw_err^2 + sqrt(log(N) / (N * bw_err)), 
               hu^2 + sqrt(log(n) / (n * hu)), 
               hx[1]^2 + sqrt(log(n) / (n * hx[1])))
  
  crit2 <- max(bw_err^2 + sqrt(log(N) / (N * bw_err)), 
               hu^2 + sqrt(log(n) / (n * hu)), 
               hx[2]^2 + sqrt(log(n) / (n * hx[2])))
  
  expect_true(err1 < crit1)
  expect_true(err2 < crit2)
  
})
