require(testthat)
test_that('error: x and xout must have same number of columns', {
  n=30 #sample size
  m=20 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(rnorm(n),rnorm(n))
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=optns),"x and xout must have same number of columns")
})

test_that('error: bandwidth must be positive', {
  n=30 #sample size
  m=20 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  bwCov = -1.0
  x=cbind(rnorm(n),rnorm(n))
  xout = matrix(c(0.25,0.5,0.9,0.5),2) #output predictor levels
  optns = list(bwCov = bwCov, metric='log_cholesky')
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=optns),"bandwidth must be positive")
})

test_that('error: M must be an array or a list', {
  #alpha=2.5
  n=30 #sample size
  m=20 # dimension of covariance matrices
  M <- matrix(1,n,3)
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=optns),"M must be an array or a list")
})


test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  n=30 #sample size
  m=20 # dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=optns),"the number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('Check correlation matrix output in the case p=2 with cross validation', {
  n=20 #sample size
  m=6 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(rnorm(n),rnorm(n))
  xout =cbind(runif(3),runif(3))#output predictor levels
  optns = list(kernel ='gauss',corrOut=TRUE,metric='log_cholesky')
  aux=LFRCovCholesky(x=x,M=M,xout=xout,optns)
  Mout=aux[[2]]
  expect_equal(sum(diag(Mout[[2]])),m)
})

test_that('Check Local Regression Simulated Setting Works (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=100000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))#cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag(exp(x[i,]))
  }
  xout=cbind(0.5,1)
  M0 = diag(exp(as.vector(xout)))
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="log_cholesky",bwCov=c(0.5,0.5)))
  aux1 = sum(abs(Cov_est$Mout[[1]]- M0))
  
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Local Regression Simulated Setting Works (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=100 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))#cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag(exp(x[i,]))
  }
  xout=cbind(0.5,1)
  M0 = diag(exp(as.vector(xout)))
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="log_cholesky"))#using CV
  aux1 = sum(abs(Cov_est$Mout[[1]]- M0))
  
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Local Regression Simulated Setting Works (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=100000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag((1+x[i,])^(2))#diag((2+x[i,])^(1/3))
  }
  
  xout=cbind(0.5,0.5)
  M0 <- diag((1+as.vector(xout))^(2))
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="cholesky",bwCov=c(0.5,0.5)))
  aux1 = sum(abs(Cov_est$Mout[[1]]-M0))
  
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})


test_that('Check Local Regression Simulated Setting Works (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=100 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag((1+x[i,])^(2))#diag((2+x[i,])^(1/3))
  }
  
  xout=cbind(0.5,0.5)
  M0 <- diag((1+as.vector(xout))^(2))
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="cholesky")) #using CV
  aux1 = sum(abs(Cov_est$Mout[[1]]-M0))
  
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})


