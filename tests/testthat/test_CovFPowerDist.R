
require(testthat)

test_that('Check M as input works', {
  #Example M input as array
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  for (i in 1:2){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Covdist=CovFPowerDist(M=M,optns=list(metric="frobenius"))$dist
  expect_equal(length(Covdist),1)
})

test_that('Check case M as input works: power metric case', {
  #Example M input as array
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  for (i in 1:2){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Covdist=CovFPowerDist(M=M,optns=list(metric="power",alpha=2))$dist
  expect_equal(length(Covdist),1)
})

test_that('Check case M as input works: power metric case', {
  #Example M input as array
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  for (i in 1:2){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Covdist=CovFPowerDist(M=M,optns=list(metric="power",alpha=0))$dist
  expect_equal(length(Covdist),1)
})

test_that('Check power metric computation: frobenius case with power call and providing alpha 1', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  Covdist=CovFPowerDist(M=M,optns=list(metric="power",alpha=1))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})

test_that('Check power metric computation: frobenius case with power call without providing alpha', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  Covdist=CovFPowerDist(M=M,optns=list(metric="power"))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})

test_that('Check power metric computation: frobenius case with power call and providing alpha 1', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  Covdist=CovFPowerDist(M=M,optns=list(metric="frobenius"))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})

test_that('Check power metric computation: alpha 2.5', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  Covdist=CovFPowerDist(M=M,optns=list(metric="power",alpha=2.5))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)^2.5-c(2.5,1.2,4.8)^2.5)^2))/2.5,0)
})

test_that('Check power metric computation: alpha 0', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  Covdist=CovFPowerDist(M=M,optns=list(metric="power",alpha=0))$dist
  expect_equal(Covdist-sqrt(sum((log(c(1,2,3))-log(c(2.5,1.2,4.8)))^2)),0)
})




