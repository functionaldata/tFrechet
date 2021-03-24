library(testthat)

test_that("Works with specified bandwith, random design", {
  set.seed(1)
  n <- 200
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- runif(n)
  err_sd <- 0.2
  xout <- seq(0.1,0.9,0.02)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- t(apply(cbind(rep(1,length(xin)), phi_true, theta_true), 1, pol2car))
  basis <- t(apply(ytrue, 1, frameSphere))
  yin_tg <- basis[,1:3] * rnorm(n, mean = 0, sd = err_sd) + 
    basis[,4:6] * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expSphere(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.15, kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- t(apply(cbind(rep(1,length(xout)), phi_true, theta_true), 1, pol2car))
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})

test_that("Works with specified bandwith, fixed design", {
  set.seed(1)
  n <- 201
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- seq(0,1,length.out = n)
  err_sd <- 0.2
  xout <- seq(0.1,0.9,0.02)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- t(apply(cbind(rep(1,length(xin)), phi_true, theta_true), 1, pol2car))
  basis <- t(apply(ytrue, 1, frameSphere))
  yin_tg <- basis[,1:3] * rnorm(n, mean = 0, sd = err_sd) + 
    basis[,4:6] * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expSphere(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.15, kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- t(apply(cbind(rep(1,length(xout)), phi_true, theta_true), 1, pol2car))
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})

test_that("Works with CV-chosen bandwith, random design", {
  set.seed(1)
  n <- 200
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- runif(n)
  err_sd <- 0.2
  xout <- seq(0.1,0.9,0.02)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- t(apply(cbind(rep(1,length(xin)), phi_true, theta_true), 1, pol2car))
  basis <- t(apply(ytrue, 1, frameSphere))
  yin_tg <- basis[,1:3] * rnorm(n, mean = 0, sd = err_sd) + 
    basis[,4:6] * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expSphere(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- t(apply(cbind(rep(1,length(xout)), phi_true, theta_true), 1, pol2car))
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})
