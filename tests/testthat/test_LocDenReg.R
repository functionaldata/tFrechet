library(testthat)

test_that("Works with fully observed distributions; p=1", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  qSup <- qbeta((1:99)/100,1/2,1/2)
  sd <- 0.1
  qin <- t(sapply(xin, function(x) {
    qnorm(c(1e-6,qSup,1-1e-6), x, sd)
  }))
  xout <- xin
  res <- LocDenReg(xin=xin, qin=qin, xout=xout, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(qin) < 1e-4)
})

test_that("Works with discrete noisy measurements; p=1", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  sd <- 0.1
  yin <- lapply(xin, function(x) {
    rnorm(1000, x, sd)
  })
  xout <- xin
  qSup <- qbeta((1:99)/100,1/2,1/2)
  res <- LocDenReg(xin=xin, yin=yin, xout=xout, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(unlist(yin)) < 2e-2)
})

test_that("Works with specifying outputGrid; p=1", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  yin <- lapply(xin, function(x) {
    rnorm(100, x, 0.01)
  })
  xout <- xin
  dSup <- seq(-0.5,1.5,0.01)
  res <- LocDenReg(xin=xin, yin=yin, xout=xout, optns=list(outputGrid = dSup))
  expect_true("dSup" %in% names(res))
  expect_equal(sum(abs(res$dSup - dSup)), 0)
})

test_that("Generates warnings when more than one of the three, yin, hin, and qin, are specified and priority order is yin, hin, qin.; p=1", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  yin <- lapply(xin, function(x) {
    rnorm(100, x, 0.01)
  })
  xout <- xin
  hin <- lapply(yin, function(y) hist(y, breaks = 50))
  qSup <- seq(0,1,0.01)
  qin <- t(sapply(xin, function(x) qbeta(qSup, x*10 + 1, 1)))
  expect_warning(res_yh <- LocDenReg(xin=xin, yin=yin, hin=hin, xout=xout))
  expect_warning(res_yq <- LocDenReg(xin=xin, yin=yin, qin=qin, xout=xout))
  expect_warning(res_hq <- LocDenReg(xin=xin, hin=hin, qin=qin, xout=xout))
  res_y <- LocDenReg(xin=xin, yin=yin, xout=xout)
  res_h <- LocDenReg(xin=xin, hin=hin, xout=xout)
  expect_equal(sum(abs(res_y$qin - res_yh$qin)),0)
  expect_equal(sum(abs(res_y$qin - res_yq$qin)),0)
  expect_equal(sum(abs(res_h$qin - res_hq$qin)),0)
})

test_that("Works with fully observed distributions and providing bandwidth; p=2", {
  set.seed(1)
  xin <- cbind(runif(200),runif(200))
  qSup <- qbeta((1:99)/100,1/2,1/2)
  sd <- 0.1
  qin=matrix(0,nrow=nrow(xin),ncol=(length(qSup)+2))
  qtrue=matrix(0,nrow=nrow(xin),ncol=(length(qSup)))
  for(i in 1:nrow(xin)){
    qin[i,]=qnorm(c(1e-6,qSup,1-1e-6), mean=xin[i,1]^2+xin[i,2], sd=sd)
    qtrue[i,]=qnorm(qSup, mean=xin[i,1]^2+xin[i,2], sd=sd)
  }
  xout <- xin
  res <- LocDenReg(xin=xin, qin=qin, xout=xout, optns = list(qSup = c(0,qSup,1)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(qin) < 1e-4)
})

test_that("Works with fully observed distributions; p=2", {
  set.seed(1)
  xin <- cbind(runif(200),runif(200))
  qSup <- qbeta((1:99)/100,1/2,1/2)
  sd <- 0.1
  qin=matrix(0,nrow=nrow(xin),ncol=(length(qSup)+2))
  qtrue=matrix(0,nrow=nrow(xin),ncol=(length(qSup)))
  for(i in 1:nrow(xin)){
    qin[i,]=qnorm(c(1e-6,qSup,1-1e-6), mean=xin[i,1]^2+xin[i,2], sd=sd)
    qtrue[i,]=qnorm(qSup, mean=xin[i,1]^2+xin[i,2], sd=sd)
  }
  xout <- xin
  res <- LocDenReg(xin=xin, qin=qin, xout=xout, optns = list(qSup = c(0,qSup,1),bwReg=c(0.05,0.05)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(qin) < 1e-2)
})
