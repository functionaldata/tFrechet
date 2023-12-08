library(testthat)

test_that("truth: great circle (no noise)", {
  n <- 101
  xin <- seq(-1,1,length.out = n)
  theta_true <- rep(pi/2,n)
  phi_true <- (xin + 1) * pi / 4
  err_sd <- 0
  ytrue <- apply( cbind( 1, phi_true, theta_true ), 1, pol2car )
  if (err_sd > 0) {
    basis <- apply( ytrue, 2, frameSphere ) # [1:3,] - 1st basis vector; [4:6,] - 2nd basis vector
    set.seed(1)
    yin_tg <- basis[1:3,] * rnorm(n, mean = 0, sd = err_sd) + 
      basis[4:6,] * rnorm(n, mean = 0, sd = err_sd)
    yin <- t( sapply( seq_len(n), function(i) expSphere( base = ytrue[,i], tg = yin_tg[,i] ) ) )
  } else {
    yin <- t( ytrue )
  }
  
  xout <- xin
  res <- GloSpheReg(xin=xin, yin=yin, xout=xout)
  expect_lt( mean(sapply(seq_len(n), function(i) SpheGeoDist(res$yout[i,], ytrue[,i]))), 5e-3 )
})

# How to generate data according to a global model?
test_that("truth: great circle + noise", {
  n <- 101
  xin <- seq(-1,1,length.out = n)
  theta_true <- rep(pi/2,n)
  phi_true <- (xin + 1) * pi / 4
  err_sd <- 0.02
  ytrue <- apply( cbind( 1, phi_true, theta_true ), 1, pol2car )
  if (err_sd > 0) {
    basis <- apply( ytrue, 2, frameSphere ) # [1:3,] - 1st basis vector; [4:6,] - 2nd basis vector
    set.seed(1)
    yin_tg <- basis[1:3,] * rnorm(n, mean = 0, sd = err_sd) + 
      basis[4:6,] * rnorm(n, mean = 0, sd = err_sd)
    yin <- t( sapply( seq_len(n), function(i) expSphere( base = ytrue[,i], tg = yin_tg[,i] ) ) )
  } else {
    yin <- t( ytrue )
  }
  
  xout <- xin
  res <- GloSpheReg(xin=xin, yin=yin, xout=xout)
  expect_lt( mean(sapply(seq_len(n), function(i) SpheGeoDist(res$yout[i,], ytrue[,i]))), 5e-3 )
})
