
require(testthat)

test_that('Check power metric computation: frobenius case with power call and providing alpha 1', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  A=M[,,1]
  B=M[,,2]
  Covdist=dist4cov(A=A,B=B,optns=list(metric="power",alpha=1))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})


test_that('Check power metric computation: frobenius case with power call without providing alpha', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  A=M[,,1]
  B=M[,,2]
  Covdist=dist4cov(A=A,B=B,optns=list(metric="power"))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})

test_that('Check power metric computation: frobenius case with power call and providing alpha 1', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  A=M[,,1]
  B=M[,,2]
  Covdist=dist4cov(A=A,B=B,optns=list(metric="frobenius"))$dist
  expect_equal(Covdist-sqrt(sum((c(1,2,3)-c(2.5,1.2,4.8))^2)),0)
})

test_that('Check power metric computation: alpha 0', {
  #Example M input as array
  m=3 # dimension of covariance matrices
  M <- array(0,c(m,m,2))
  M[,,1]=diag(c(1,2,3))
  M[,,2]=diag(c(2.5,1.2,4.8))
  A=M[,,1]
  B=M[,,2]
  Covdist=dist4cov(A=A,B=B,optns=list(metric="power",alpha=0))$dist
  expect_equal(Covdist-sqrt(sum((log(c(1,2,3))-log(c(2.5,1.2,4.8)))^2)),0)
})

test_that('Check The Distance Works (accurate estimate to the true target) for Cholesky metric', {
  set.seed(1234321)
  p <- 3
  a <- matrix(rnorm(p*p),p,p)
  a <- chol(t(a) %*% (a))
  #a[lower.tri(matrix(rnorm(p*p),p,p))] <- 0
  
  b <- matrix(rnorm(p*p),p,p)
  b <- chol(t(b) %*% (b))
  #b[lower.tri(matrix(rnorm(p*p),p,p))] <- 0
  
  true.dist <- sqrt(sum((a-b)^2))
  A <- t(a) %*% a
  B <- t(b) %*% b
  optns = list(metric='cholesky')
  res <- dist4cov(A=A,B=B,optns=optns)
  aux1 <- abs(res$dist-true.dist)
  aux1
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})


test_that('Check The Distance Works (accurate estimate to the true target) for Log Cholesky metric', {
  set.seed(1234321)
  p <- 3
  a <- abs(diag(rnorm(p)))
  b <- abs(diag(rnorm(p)))
  
  true.dist <- sqrt(sum((log(diag(a))-log(diag(b)))^2))
  A <- t(a) %*% a
  B <- t(b) %*% b
  optns = list(metric='log_cholesky')
  res <- dist4cov(A=A,B=B,optns=optns)
  aux1 <- abs(res$dist-true.dist)
  aux1
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})



