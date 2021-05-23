require(testthat)

test_that('error: x and xout must have same number of columns', {
  set.seed(1234321)
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
  expect_error(LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3)),"x and xout must have the same number of columns")
})

test_that('error: x and xout must have same number of columns', {
  set.seed(1234321)
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
  expect_error(LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="frobenius")),"x and xout must have the same number of columns")
})

test_that('error: x and xout must have same number of columns', {
  set.seed(1234321)
  #Example M input
  n=30 #sample size
  m=30 #dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
   y0=rnorm(m)
   aux<-15*diag(m)+y0%*%t(y0)
   M[,,i]<-aux
  }
  x=cbind(rnorm(n),rnorm(n))
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=0)),"x and xout must have the same number of columns")
})


test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  set.seed(1234321)
  n=30 #sample size
  m=30 #dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
   y0=rnorm(m)
   aux<-15*diag(m)+y0%*%t(y0)
   M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=0)),"The number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  set.seed(1234321)
  n=30 #sample size
  m=30 #dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  expect_error(LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power")),"The number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('Check correlation matrix output in the case p=1 with cross validation', {
  set.seed(1234321)
  n=30 #sample size
  m=30 #dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
   y0=rnorm(m)
   aux<-15*diag(m)+y0%*%t(y0)
   M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  aux=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=TRUE,metric="power",alpha=0))$Mout
  expect_equal(sum(diag(aux[[1]])),m)
})

test_that('Check correlation matrix output in the case p=1 with cross validation', {
  set.seed(1234321)
  n=30 #sample size
  m=15 #dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  aux=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=TRUE))$Mout
  expect_equal(sum(diag(aux[[1]])),m)
})

test_that('Check case y as input works power case', {
  set.seed(1234321)
  n=100             # sample size
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
  Cov_est=LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3,bwCov=0.1))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case y as input works frobenius', {
  set.seed(1234321)
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
  Cov_est=LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="frobenius",bwCov=0.1))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case y as input works frobenius with CV', {
  set.seed(1234321)
  n=200             # sample size
  t=seq(0,1,length.out=30)       # length of data
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
  Cov_est=LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="frobenius"))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case M as input works', {
  set.seed(1234321)
  n=30 #sample size
  m=15 #dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
   y0=rnorm(m)
   aux<-15*diag(m)+y0%*%t(y0)
   M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check Local 1D Regression Frobenius metric (accurate estimate to the true target)', {
  set.seed(1234321)
  n=100 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- runif(n,min=-1,max=1)
  for (i in 1:n){
    M[,,i]<-diag(c(2+x[i],4-x[i]))
  }
  xout=0.5
  
  Cov_est=LFRCovPower(x=x,M=M,xout=xout,optns=list(alpha=1,corrOut=FALSE)) #using CV
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2.5,3.5))))
  
  Cov_est=LFRCovPower(x=x,M=M,xout=xout,optns=list(alpha=1,corrOut=FALSE,bwCov=c(0.5,0.5)))
  aux1=aux1+sum(abs(Cov_est$Mout[[1]]-diag(c(2.5,3.5))))
  
  Cov_est=LFRCov(x=x,M=M,xout=xout,optns=list(corrOut=FALSE)) #using CV
  aux2=sum(abs(Cov_est$Mout[[1]]-diag(c(2.5,3.5))))
  
  Cov_est=LFRCov(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,bwCov=c(0.5,0.5)))
  aux2=aux2+sum(abs(Cov_est$Mout[[1]]-diag(c(2.5,3.5))))
  
  if(aux1+aux2<=1e-10){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Local 2D Regression Frobenius metric (accurate estimate to the true target)', {
  set.seed(1234321)
  n=200 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  
  xout=cbind(0,0.5)
  Cov_est=LFRCovPower(x=x,M=M,xout=xout,optns=list(alpha=1,corrOut=FALSE)) #using CV
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  Cov_est=LFRCov(x=x,M=M,xout=xout,optns=list(corrOut=FALSE)) #using CV
  aux2=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  
  if(aux1+aux2<=1e-10){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Local Regression frobenius metric (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=200 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,]))
  }
  
  xout=cbind(0,0.5)
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="frobenius",bwCov=c(0.5,0.5)))
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="frobenius")) #using CV
  aux1=aux1+sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5))))
  
  if(aux1<=1e-10){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Local Regression power metric (accurate estimate to the true target) on main Local function', {
  set.seed(1234321)
  n=200 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<-diag((2+x[i,])^(1/3))
  }
  xout=cbind(0,0.5)
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3,bwCov=c(0.5,0.5)))
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5)^(1/3))))
  
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3)) #using CV
  aux1=aux1+sum(abs(Cov_est$Mout[[1]]-diag(c(2,2.5)^(1/3))))
  
  if(aux1<=1e-10){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check case y as input works fractional power', {
  set.seed(1234321)
  n=100             # sample size
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
  Cov_est=LocCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=1/2,bwCov=0.1))
  expect_equal(length(Cov_est$Mout),3)
})

test_that('Check case M as input works with fractional power', {
  set.seed(1234321)
  n=30 #sample size
  m=15 #dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-15*diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=matrix(rnorm(n),n)
  xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
  Cov_est=LocCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=1/2))
  expect_equal(length(Cov_est$Mout),3)
})
