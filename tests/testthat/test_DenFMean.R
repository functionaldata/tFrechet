library(testthat)

test_that("Works with noisy yin input", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  sd <- 0.01
  yin <- lapply(xin, function(x) {
    rnorm(1000, x, sd)
  })
  qSup <- qbeta((1:99)/100,1/2,1/2)
  res <- DenFMean(yin=yin, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(sum((res$qout[2:100] - colMeans(qtrue))^2) < 1e-4)
})

test_that("Works with qin input", {
  set.seed(1)
  xin <- seq(0,1,0.01)
  qSup <- qbeta((1:99)/100,1/2,1/2)
  sd <- 0.1
  qin <- t(sapply(xin, function(x) {
    qnorm(c(1e-6,qSup,1-1e-6), x, sd)
  }))
  res <- DenFMean(qin=qin, optns = list(qSup = c(0,qSup,1)))
  qtrue <- t(sapply(xin, function(x) qnorm(qSup, mean = x, sd = sd)))
  expect_true(sum((res$qout[2:100] - colMeans(qtrue))^2) < 1e-8)
})


test_that("Weighted barycenter unequal weights works", {
  set.seed(1)
  n=50
  xin <- seq(0,1,length.out=n)
  qSup=seq(1e-6,1-1e-6,length.out=100)
  qin=t(sapply(xin,function(x){qnorm(qSup,mean=0,sd=1)+x}))
  weightsF=n:1
  weightsF=weightsF/sum(weightsF)
  res=DenFMean(qin=qin, optns = list(qSup = c(0,qSup[2:(length(qSup)-1)],1),weights=weightsF))
  barycenter_weighted=colSums(diag(weightsF)%*%qin[,2:(length(qSup)-1)])
  expect_true(sum((res$qout[2:(length(qSup)-1)] - barycenter_weighted)^2) < 1e-8)
})

test_that("Weighted barycenter using equal weights works", {
  set.seed(1)
  n=50
  xin <- seq(0,1,length.out=n)
  qSup=seq(1e-6,1-1e-6,length.out=100)
  qin=t(sapply(xin,function(x){qnorm(qSup,mean=0,sd=1)+x}))
  weightsF=rep(1/n,n)
  res=DenFMean(qin=qin, optns = list(qSup = c(0,qSup[2:(length(qSup)-1)],1),weights=weightsF))
  barycenter_weighted=colSums(diag(weightsF)%*%qin[,2:(length(qSup)-1)])
  expect_true(sum((res$qout[2:(length(qSup)-1)] - barycenter_weighted)^2) < 1e-8)
})



