require(testthat)
test_that('A and B must be of matrix class', {
  A <- c(1,10)
  B <- matrix(1,10,10)
  optns = list(metric='log_cholesky')
  expect_error(DistCholesky(A,B,optns=optns),"A and B must be of matrix class")
})

test_that('Both A and B must have the same dimension', {
  A <- matrix(1,9,9)
  B <- matrix(1,10,10)
  optns = list(metric='log_cholesky')
  expect_error(DistCholesky(A,B,optns=optns),"Both A and B must have the same dimension")
})

test_that('Both A and B must be square matrices', {
  A <- matrix(1,10,9)
  B <- matrix(1,10,9)
  optns = list(metric='log_cholesky')
  expect_error(DistCholesky(A,B,optns=optns),"Both A and B must be square matrices")
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
  res <- DistCholesky(A,B,optns=optns)
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
  res <- DistCholesky(A,B,optns=optns)
  aux1 <- abs(res$dist-true.dist)
  aux1
  if(aux1<=1e-3){
    flag=1
  }else{
    flag=0
  }
  expect_equal(flag,1)
})

