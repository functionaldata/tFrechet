
require(testthat)

test_that('Check M as input works', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Fmean=CovFMean(M=M,optns=list(metric="frobenius"))
  expect_equal(length(Fmean$Mout),1)
})

test_that('Check case M as input works: power metric case', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Fmean=CovFMean(M=M,optns=list(metric="power",alpha=2))
  expect_equal(length(Fmean$Mout),1)
})

test_that('Check case M as input works: cholesky metric case', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Fmean=CovFMean(M=M,optns=list(metric="cholesky",alpha=2))
  expect_equal(length(Fmean$Mout),1)
})

test_that('Check case M as input works: log-cholesky metric case', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  Fmean=CovFMean(M=M,optns=list(metric="log_cholesky",alpha=2))
  expect_equal(length(Fmean$Mout),1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target)', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  Fmean=CovFMean(M=M,optns=list(metric="frobenius"))
  aux1=sum(abs(Fmean$Mout[[1]]-diag(c(2,2))))
  if(aux1<=0.05){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) power case', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,])^(1/3))
  }
  Fmean=CovFMean(M=M,optns=list(metric="power",alpha=3))

  aux1=sum(abs(Fmean$Mout[[1]]-diag(c(2,2)^(1/3))))
  if(aux1<=0.01){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) log_cholesky case', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(rnorm(n),rnorm(n))
  for (i in 1:n){
    M[,,i]<- diag(exp(x[i,]))
  }

  xout=cbind(0,0)
  M0 = diag(exp(as.vector(xout)))
  Fmean=CovFMean(M=M,optns=list(metric="log_cholesky"))
  aux1=sum(abs(Fmean$Mout[[1]]-M0))

  if(aux1<=0.05){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) cholesky case', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag((1+x[i,])^(2))
  }
  xout=cbind(0,0)
  M0 <- diag((1+as.vector(xout))^(2))
  Fmean=CovFMean(M=M,optns=list(metric="cholesky"))
  aux1=sum(abs(Fmean$Mout[[1]]-M0))

  if(aux1<=0.05){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})




test_that('Check unweighted Frobenius case works', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  weightsF=rep(1/n,n)
  Fmean=CovFMean(M=M,optns=list(metric="frobenius",weights=weightsF))
  cont=matrix(0,nrow=m,ncol=m)
  for(i in 1:n){
    cont=cont+M[,,i]
  }
  cont=cont/n
  expect_true(sum(abs(cont-Fmean$Mout[[1]]))<1e-8)
})

test_that('Check weighted Frobenius case works', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  weightsF=n:1
  weightsF=weightsF/sum(weightsF)
  Fmean=CovFMean(M=M,optns=list(metric="frobenius",weights=weightsF))
  cont=matrix(0,nrow=m,ncol=m)
  for(i in 1:n){
    cont=cont+M[,,i]*weightsF[i]
  }
  expect_true(sum(abs(cont-Fmean$Mout[[1]]))<1e-5)
})

test_that('Check weighted Frobenius using power input case works', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  weightsF=n:1
  weightsF=weightsF/sum(weightsF)
  Fmean=CovFMean(M=M,optns=list(metric="power",alpha=1,weights=weightsF))
  cont=matrix(0,nrow=m,ncol=m)
  for(i in 1:n){
    cont=cont+M[,,i]*weightsF[i]
  }
  expect_true(sum(abs(cont-Fmean$Mout[[1]]))<1e-5)
})

test_that('Check weighted general power input case works', {
  set.seed(1234321)
  n=500 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  idx=1:n%%2
  weightsF=idx/sum(idx)
  Fmean=CovFMean(M=M,optns=list(metric="power",alpha=2,weights=weightsF))
  cont=matrix(0,nrow=m,ncol=m)
  for(i in 1:n){
    cont=cont+(M[,,i]%*%M[,,i])*weightsF[i]
  }
  expect_true(sum(abs(sqrt(cont)-Fmean$Mout[[1]]))<1e-4)
})

test_that('Check unweighted Frobenius using power input case works', {
  set.seed(1234321)
  n=5000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  weightsF=rep(1/n,n)
  Fmean=CovFMean(M=M,optns=list(metric="power",alpha=1,weights=weightsF))
  cont=matrix(0,nrow=m,ncol=m)
  for(i in 1:n){
    cont=cont+M[,,i]*weightsF[i]
  }
  expect_true(sum(abs(cont-Fmean$Mout[[1]]))<1e-5)
})
