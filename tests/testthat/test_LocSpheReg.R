library(testthat)

test_that("Works with specified bandwith, random design", {
  set.seed(1)
  n <- 200
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- runif(n)
  err_sd <- 0.2
  xout <- seq(0,1,length.out = 51)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  basis <- bss_tgsp(phi = phi_true, theta = theta_true)
  yin_tg <- basis$b1 * rnorm(n, mean = 0, sd = err_sd) + 
    basis$b2 * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expS2(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.15, kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})

test_that("Works with specified bandwith, fixed design", {
  set.seed(1)
  n <- 201
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- seq(0,1,length.out = n)
  err_sd <- 0.2
  xout <- seq(0,1,length.out = 51)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  basis <- bss_tgsp(phi = phi_true, theta = theta_true)
  yin_tg <- basis$b1 * rnorm(n, mean = 0, sd = err_sd) + 
    basis$b2 * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expS2(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.15, kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})

test_that("Works with CV-chosen bandwith, random design", {
  set.seed(1)
  n <- 200
  # simulate the data according to the simulation in Petersen & Müller (2019)
  xin <- runif(n)
  err_sd <- 0.2
  xout <- seq(0,1,length.out = 51)
  
  phi_true <- acos(xin)
  theta_true <- pi * xin
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  basis <- bss_tgsp(phi = phi_true, theta = theta_true)
  yin_tg <- basis$b1 * rnorm(n, mean = 0, sd = err_sd) + 
    basis$b2 * rnorm(n, mean = 0, sd = err_sd)
  yin <- t(sapply(seq_len(n), function(i) expS2(base = ytrue[i,], tg = yin_tg[i,])))
  
  res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(kernel = "epan"))
  phi_true <- acos(xout)
  theta_true <- pi * xout
  ytrue <- pol2car(phi = phi_true, theta = theta_true)
  expect_true(mean(sapply(seq_along(xout), function(i) SpheGeoDist(res$yout[i,], ytrue[i,])^2)) < 3e-3)
})
