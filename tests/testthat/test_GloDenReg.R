library(testthat)

test_that("Works with fully observed distributions", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  qSup <- qbeta((1:99)/100,1/2,1/2)
  sd <- 0.1
  qin <- t(sapply(xin, function(x) {
    qnorm(c(1e-6,qSup,1-1e-6), x, sd)
  }))
  xout <- xin
  res <- GloDenReg(xin=xin, qin=qin, xout=xout, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(qin) < 1e-4)
})

test_that("Works with discrete noisy measurements", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  sd <- 0.1
  yin <- lapply(xin, function(x) {
    rnorm(1000, x, sd)
  })
  xout <- xin
  qSup <- qbeta((1:99)/100,1/2,1/2)
  res <- GloDenReg(xin=xin, yin=yin, xout=xout, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(mean(sqrt(apply((qtrue - res$qout[,-c(1,length(res$qSup))])^2,
                              1, pracma::trapz, x = qSup))) / mean(unlist(yin)) < 2e-2)
})

test_that("Works with specifying outputGrid", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  yin <- lapply(xin, function(x) {
    rnorm(100, x, 0.01)
  })
  xout <- xin
  dSup <- seq(-0.5,1.5,0.01)
  res <- GloDenReg(xin=xin, yin=yin, xout=xout, optns=list(outputGrid = dSup))
  expect_true("dSup" %in% names(res))
  expect_equal(sum(abs(res$dSup - dSup)), 0)
})

test_that("Generates warnings when more than one of the three, yin, hin, and qin, are specified and priority order is yin, hin, qin.", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  yin <- lapply(xin, function(x) {
    rnorm(100, x, 0.01)
  })
  xout <- xin
  hin <- lapply(yin, function(y) hist(y, breaks = 50))
  qSup <- seq(0,1,0.01)
  qin <- t(sapply(xin, function(x) qbeta(qSup, x*10 + 1, 1)))
  expect_warning(res_yh <- GloDenReg(xin=xin, yin=yin, hin=hin, xout=xout))
  expect_warning(res_yq <- GloDenReg(xin=xin, yin=yin, qin=qin, xout=xout))
  expect_warning(res_hq <- GloDenReg(xin=xin, hin=hin, qin=qin, xout=xout))
  res_y <- GloDenReg(xin=xin, yin=yin, xout=xout)
  res_h <- GloDenReg(xin=xin, hin=hin, xout=xout)
  expect_equal(res_y, res_yh)
  expect_equal(res_y, res_yq)
  expect_equal(res_h, res_hq)
})

test_that("Generates warnings when R squared is required with qin but no qSup.", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  xout <- xin
  qSup <- seq(0,1,0.01)
  qin <- t(sapply(xin, function(x) qbeta(qSup, x*10 + 1, 1)))
  expect_warning(GloDenReg(xin=xin, qin=qin, xout=xout, optns=list(Rsq = TRUE)))
})
