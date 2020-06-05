library(testthat)


testthat::test_that('user provided bw for mean works',{
  set.seed(1)
  x0<-seq(0,1,length.out = 100)
  Y<-rbeta(100,5,3)
  bwd<-0.1
  f1 <- CreateDensity(y=Y,optns = list(userBwMu=bwd,outputGrid=x0))
  f2 <- CreateDensity(histogram = hist(Y),optns = list(userBwMu=bwd,outputGrid=x0))
  binY <- c(0,sort(runif(19)),1)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(userBwMu=bwd,outputGrid=x0))
  expect_equal(f1$bw,bwd)
  expect_equal(f2$bw,bwd)
  expect_equal(f3$bw,bwd)
  #expect_equal(true_dens,f1$y,tolerance = 0.25)
})

testthat::test_that('user provided output grid works',{
  set.seed(1)
  x0<-c(seq(0,0.5,length.out = 50),seq(0.51,1,length.out = 25))
  Y<-rbeta(100,4,5)
  binY <- c(0,sort(runif(9)),1)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  
  f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
  f2 <- CreateDensity(histogram = hist(Y),optns = list(outputGrid=x0))
  f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
  expect_equal(f1$x,x0)
  expect_equal(f2$x,x0)
  expect_equal(f3$x,x0)
})

testthat::test_that('incorrect input data format gives errors',{
  set.seed(1)
  x0<-seq(0,1,length.out = 100)
  Y<-rbeta(100,4,5)
  binY <- c(0,sort(runif(9)),1)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  
  expect_error(CreateDensity(y=hist(Y),optns = list(outputGrid=x0)))
  expect_error(CreateDensity(histogram = Y,optns = list(outputGrid=x0)))
  expect_error(CreateDensity(freq=freqY,optns = list(outputGrid=x0)))
})

testthat::test_that('bw selection works/get reasonable density estimates (beta distribution)',{
  set.seed(1)
  x0=seq(0,1,length.out = 100)
  Y=rbeta(100,3,2)
  binY <- c(0,sort(runif(19,0,1)),1)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  
  f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
  f2 <- CreateDensity(histogram = hist(Y),optns = list(outputGrid=x0))
  f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
  
  true_dens=dbeta(x0,3,2)
  expect_equal(f1$y,true_dens,tolerance=0.2)
  expect_equal(f2$y,true_dens,tolerance=0.2)
  expect_equal(f3$y,true_dens,tolerance=0.2)
})

testthat::test_that('bw selection works/get reasonable density estimates (normal distribution)',{
  set.seed(1)
  x0=seq(-3,3,length.out = 100)
  Y=rnorm(200)
  binY <- c(-3,sort(runif(19,-3,3)),3)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  
  f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
  f2 <- CreateDensity(histogram = hist(Y),optns = list(outputGrid=x0))
  f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
  
  true_dens=dnorm(x0,0,1)
  expect_equal(f1$y,true_dens,tolerance=0.2)
  expect_equal(f2$y,true_dens,tolerance=0.2)
  expect_equal(f3$y,true_dens,tolerance=0.2)
})

testthat::test_that('bw selection works/get reasonable density estimates (exponential distribution)',{
  set.seed(1)
  x0=seq(0,5,length.out = 100)
  Y=rexp(200)
  binY <- c(0,sort(runif(19,0,5)),5)
  freqY <- c()
  for (i in 1:(length(binY)-1)) {
    freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
  }
  
  f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
  f2 <- CreateDensity(histogram = hist(Y),optns = list(outputGrid=x0))
  f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
  
  true_dens=dexp(x0)
  expect_equal(f1$y,true_dens,tolerance=0.2)
  expect_equal(f2$y,true_dens,tolerance=0.2)
  expect_equal(f3$y,true_dens,tolerance=0.2)
})
