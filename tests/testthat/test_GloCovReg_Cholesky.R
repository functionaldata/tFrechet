require(testthat)
test_that('error: x must be a matrix', {
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  for (i in 1:n){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=rnorm(n) #vector of predictor values
  xout=matrix(runif(3)) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=optns),"x must be a matrix")
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
  x=cbind(rnorm(n),rnorm(n)) #vector of predictor values
  xout=matrix(runif(3)) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=optns),"x and xout must have same number of columns")
})


test_that('error: M must be an array or a list', {
  #Example array as input for M
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- matrix(1,10,5)
  x=matrix(rnorm(n)) #vector of predictor values
  xout=matrix(runif(3)) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=optns),"M must be an array or a list")
})


test_that('error: the number of rows of x must be the same as the number of covariance matrices in M', {
  #Example array as input for M
  n=10 #sample size
  m=5 # dimension of covariance matrices
  M <- array(0,c(m,m,n+1))
  for (i in 1:(n+1)){
    y0=rnorm(m)
    aux<-diag(m)+y0%*%t(y0)
    M[,,i]<-aux
  }
  x=matrix(rnorm(n)) #vector of predictor values
  xout=matrix(runif(3)) #output predictor levels
  optns = list(metric='log_cholesky')
  expect_error(GloCovReg(x=x,M=M,xout=xout,optns=optns),"the number of rows of x must be the same as the number of covariance matrices in M")
})

test_that('Check returned covariance matrices list has the right dimension', {
  n=30 #sample size
  m=5  #dimension of covariance matrices
  x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
  M <- array(0,c(m,m,n))
  a = rnorm(m); b = rnorm(m)
  A = diag(m)+a%*%t(a);
  B = diag(m)+3*b%*%t(b);
  for (i in 1:n){
    aux <- x[i,1]*A + x[i,2]**2*B
    M[,,i] <- aux %*% t(aux)
  }
  xout=cbind(runif(5),runif(5)) #output predictor levels
  aux <- GFRCovCholesky(x=x, M=M, xout=xout)$Mout
  expect_equal(length(aux),dim(xout)[1])
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) on main Global function', {
  set.seed(1234321)
  n=10000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(rnorm(n),rnorm(n))
  for (i in 1:n){
    M[,,i]<- diag(exp(x[i,]))
  }
  
  xout=cbind(0.5,1)
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(metric='log_cholesky',corrOut=FALSE))
  aux1=sum(abs(Cov_est$Mout[[1]]-diag(exp(as.vector(xout))) ))
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

test_that('Check Global Regression Simulated Setting Works (accurate estimate to the true target) on main Global function', {
  set.seed(1234321)
  n=10000 #sample size
  m=2 # dimension of covariance matrices
  M <- array(0,c(m,m,n))
  x<- cbind(runif(n,min=-1,max=1),runif(n,min=-1,max=1))
  for (i in 1:n){
    M[,,i]<- diag((1+x[i,])^(2))
  }
  
  xout=cbind(0.5,0.5)
  M0 <- diag((1+as.vector(xout))^(2))
  Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric='cholesky',alpha=3))
  aux1 = sum(abs(Cov_est$Mout[[1]]- M0 ))
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})








