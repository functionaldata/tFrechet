require(testthat)

test_that('error: metric choice not supported', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(GloCorReg(x, y, xOut, optns = list(metric = 'procrustes')), 
               'metric choice not supported')
})

test_that('error: alpha must be non-negative', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(GloCorReg(x, y, xOut, optns = list(metric = 'frobenius', alpha = -1)), 
               'alpha must be non-negative')
})

test_that('error: x and xOut must have the same number of columns', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- cbind(seq(0.1, 0.9, length.out = 9), seq(0.1, 0.9, length.out = 9))
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(GloCorReg(x, y, xOut, optns = list(metric = 'frobenius')), 
               'x and xOut must have the same number of columns')
})

test_that('error: the number of rows in x must be the same as the number of correlation matrices in y', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  x <- runif(n+1, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(GloCorReg(x, y, xOut, optns = list(metric = 'frobenius')), 
               'the number of rows in x must be the same as the number of correlation matrices in y')
})

test_that('fitted matrices are correlation matrices in the case p=2 with Frobenius metric', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- cbind(seq(0.1, 0.9, length.out = 9), seq(0.1, 0.9, length.out = 9))
  x <- cbind(runif(n, min = 0, max = 1), runif(n, min = 0, max = 1))
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i, 1], shape2 = 1-x[i, 1])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  res <- GloCorReg(x, y, xOut)
  expect_equal(unique(sapply(res$fit, function(fiti) sum(diag(fiti)))), q)
})

test_that('fitted matrices are correlation matrices in the case p=2 with Euclidean power metric', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- cbind(seq(0.1, 0.9, length.out = 9), seq(0.1, 0.9, length.out = 9))
  x <- cbind(runif(n, min = 0, max = 1), runif(n, min = 0, max = 1))
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i, 2], shape2 = 1-x[i, 2])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  res <- GloCorReg(x, y, xOut, optns = list(metric = 'power'))
  expect_equal(unique(sapply(res$fit, function(fiti) sum(diag(fiti)))), q)
})

test_that('global regression simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- cbind(seq(0.1, 0.9, length.out = 9), seq(0.1, 0.9, length.out = 9))
  x <- cbind(runif(n, min = 0, max = 1), runif(n, min = 0, max = 1))
  y <- list()
  for(i in 1:n){
    yVec <- c(rep(x[i, 1], 22), rep(x[i, 2], 23))
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  res <- GloCorReg(x, y, xOut)
  eps <- sqrt(sum((res$fit[[1]]-y[[1]])^2))
  expect_equal(eps<=1e-3, TRUE)
})

test_that('global regression simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- cbind(seq(0.1, 0.9, length.out = 9), seq(0.1, 0.9, length.out = 9))
  x <- cbind(runif(n, min = 0, max = 1), runif(n, min = 0, max = 1))
  y <- list()
  for(i in 1:n){
    yVec <- c(rep(x[i, 1], 22), rep(x[i, 2], 23))
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  y <- array(unlist(y), c(q, q, n))
  res <- GloCorReg(x, y, xOut)
  eps <- sqrt(sum((res$fit[[1]]-y[, , 1])^2))
  expect_equal(eps<=1e-3, TRUE)
})

test_that('global regression simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- c(0.1, 0.2)
  x <- cbind(runif(n, min = 0, max = 1), runif(n, min = 0, max = 1))
  y <- list()
  for(i in 1:n){
    yVec <- c(rep(x[i, 1], 22), rep(x[i, 2], 23))
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  y <- array(unlist(y), c(q, q, n))
  res <- GloCorReg(x, y, xOut)
  eps <- sqrt(sum((res$fit[[1]]-y[, , 1])^2))
  expect_equal(eps<=1e-3, TRUE)
})
