
require(testthat)
test_that('error: x and xout must have same number of columns', {
  n=200             # sample size
  t=seq(0,1,length.out=100)       # length of data
  x = cbind(rnorm(n),rnorm(n))
  theta1 = theta2 = array(0,n)
  for(i in 1:n){
    theta1[i] = rnorm(1,x[i],x[i]^2)
    theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
  }
  y = matrix(0,n,length(t))
  phi1 = sqrt(3)*t
  phi2 = sqrt(6/5)*(1-t/2)
  y = theta1%*%t(phi1) + theta2 %*% t(phi2)
  xout = matrix(c(0.25,0.5,0.75),3)
  expect_error(GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3)),"x and xout must have the same number of columns")
})

test_that('error: x and xout must have same number of columns', {
  n=200             # sample size
  t=seq(0,1,length.out=100)       # length of data
  x = cbind(rnorm(n),rnorm(n))
  theta1 = theta2 = array(0,n)
  for(i in 1:n){
    theta1[i] = rnorm(1,x[i],x[i]^2)
    theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
  }
  y = matrix(0,n,length(t))
  phi1 = sqrt(3)*t
  phi2 = sqrt(6/5)*(1-t/2)
  y = theta1%*%t(phi1) + theta2 %*% t(phi2)
  xout = matrix(c(0.25,0.5,0.75),3)
  expect_error(GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE)),"x and xout must have the same number of columns")
})



test_that('error: x and xout must have same number of columns', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=0)),"x and xout must have the same number of columns")
})

test_that('error: x and xout must have same number of columns', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE)),"x and xout must have the same number of columns")
})

test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3)),"The number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE)),"The number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('Check correlation matrix output in the case p=2 with cross validation', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  aux=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=TRUE,metric="power",alpha=3))$Mout
  expect_equal(sum(diag(aux[[1]])),m)
})

test_that('Check correlation matrix output in the case p=2 with cross validation', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=TRUE))
  aux=Cov_est$Mout
  expect_equal(sum(diag(aux[[1]])),m)
})


test_that('Check case y as input works', {
  n=200             # sample size
  t=seq(0,1,length.out=100)       # length of data
  x = matrix(runif(n),n)
  theta1 = theta2 = array(0,n)
  for(i in 1:n){
    theta1[i] = rnorm(1,x[i],x[i]^2)
    theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
  }
  y = matrix(0,n,length(t))
  phi1 = sqrt(3)*t
  phi2 = sqrt(6/5)*(1-t/2)
  y = theta1%*%t(phi1) + theta2 %*% t(phi2)
  xout = matrix(c(0.25,0.5,0.75),3)
  Cov_est=GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case y as input works', {
  n=200             # sample size
  t=seq(0,1,length.out=100)       # length of data
  x = matrix(runif(n),n)
  theta1 = theta2 = array(0,n)
  for(i in 1:n){
    theta1[i] = rnorm(1,x[i],x[i]^2)
    theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
  }
  y = matrix(0,n,length(t))
  phi1 = sqrt(3)*t
  phi2 = sqrt(6/5)*(1-t/2)
  y = theta1%*%t(phi1) + theta2 %*% t(phi2)
  xout = matrix(c(0.25,0.5,0.75),3)
  Cov_est=GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case M as input works', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case M as input works', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  xout=cbind(runif(3),runif(3)) #output predictor levels
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target)', {
  set.seed(1234321)
  n=100000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  
  xout=cbind(0,0.5)
  Cov_est=GFRCov(x=x,M=M,xout=xout,optns=list(corrOut=FALSE))
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  Cov_est=GFRCovPower(x=x,M=M,xout=xout,optns=list(alpha=1,corrOut=FALSE))
  aux2=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  if(aux1+aux2<=0.001){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) on main Global function', {
  set.seed(1234321)
  n=100000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,])^(1/3))
  }
  
  xout=cbind(0,0.5)
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric='power',alpha=3))
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5)^(1/3))))
  if(aux1<=0.001){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})
