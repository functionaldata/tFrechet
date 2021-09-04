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
  expect_error(LocCorReg(x, y, xOut, optns = list(metric = 'procrustes')), 
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
  expect_error(LocCorReg(x, y, xOut, optns = list(metric = 'frobenius', alpha = -1)), 
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
  expect_error(LocCorReg(x, y, xOut, optns = list(metric = 'frobenius')), 
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
  expect_error(LocCorReg(x, y, xOut, optns = list(metric = 'frobenius')), 
               'the number of rows in x must be the same as the number of correlation matrices in y')
})

test_that('fitted matrices are correlation matrices in the case p=2 with Frobenius metric', {
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
  res <- LocCorReg(x, y, xOut, optns = list(bw = 0.08))
  expect_equal(unique(sapply(res$fit, function(fiti) sum(diag(fiti)))), q)
})

test_that('fitted matrices are correlation matrices in the case p=2 with Euclidean power metric', {
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
  res <- LocCorReg(x, y, xOut, optns = list(metric = 'power', alpha = 0.5, bw = 0.08))
  expect_equal(unique(sapply(res$fit, function(fiti) sum(diag(fiti)))), q)
})

test_that('local regression simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rep(x[i], d)
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  res <- LocCorReg(x, y, xOut)
  eps <- mean(sapply(1:n, function(i) sqrt(sum((res$fit[[i]]-y[[i]])^2))))
  expect_equal(eps<=1e-3, TRUE)
})

test_that('error: local method is designed to work in low dimensional case (p is either 1 or 2)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  xOut <- cbind(xOut, xOut, xOut)
  x <- runif(n, min = 0, max = 1)
  x <- cbind(x, x, x)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i, 1], shape2 = 1-x[i, 1])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(LocCorReg(x, y, xOut), 
               'local method is designed to work in low dimensional case')
})

test_that('error: dimension of bandwidth does not agree with x', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  xOut <- cbind(xOut, xOut)
  x <- runif(n, min = 0, max = 1)
  x <- cbind(x, x)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i, 1], shape2 = 1-x[i, 1])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(LocCorReg(x, y, xOut, optns = list(bw = 0.08)), 
               'dimension of bandwidth does not agree with x')
})

test_that('error: bandwidth must be positive', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  expect_error(LocCorReg(x, y, optns = list(metric = 'power', bw = -0.08)), 
               'bandwidth must be positive')
})

test_that('local regression simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  n <- 100
  q <- 10
  d <- q*(q-1)/2
  xOut <- seq(0.1, 0.9, length.out = 9)
  x <- runif(n, min = 0, max = 1)
  y <- list()
  for(i in 1:n){
    yVec <- rep(x[i], d)
    y[[i]] <- matrix(0, nrow = q, ncol = q)
    y[[i]][lower.tri(y[[i]])] <- yVec
    y[[i]] <- y[[i]] + t(y[[i]])
    diag(y[[i]]) <- 1
  }
  y <- array(unlist(y), c(q, q, n))
  res <- LocCorReg(x, y, xOut)
  eps <- sqrt(sum((res$fit[[1]]-y[, , 1])^2))
  expect_equal(eps<=1e-3, TRUE)
})
