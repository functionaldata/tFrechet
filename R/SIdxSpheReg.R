library(manifold)
library(trust) 

#### Main functions ####

#### Main implementation function that returns 
#### the estimated direction parameter (unit vector) 
SIdxSpheReg <- function(xin, yin, bw = NULL, M = NULL, ker = ker_gauss, iter = 1000,
                        verbose = T) {
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## yin: n by m matrix of output (m: dimension of sphere; S^{m-1})
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  ## iter: number of direction generation
  ## verbose: print the iteration counts?
  
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  if (!is.matrix(xin)){
    stop("xin should be matrix.")
  }
  if(!is.matrix(yin)){
    stop("yin should be matrix.")
  }
  
  p <- ncol(xin)
  
  
  ## Parameter (bandwidth, bin size) choice using cross-validation
  needParam <- (is.null(M) | is.null(bw))
  
  if (needParam) {
    param <- SpheTuning(xin, yin, normalize(rep(1,p)))
  }
  
  # bw2 depends on bw
  if (is.null(bw)) {
    bw2 <- param[1]
  } else {
    bw2 <- bw
  }
  
  # M2 depends on M
  if (is.null(M)) {
    M2 <- param[2]
  } else {
    M2 <- M
  }
  
  
  ## Generate Equal Grid over the angles
  coords_mat <- matrix(rnorm(iter * p), nrow = iter, ncol = p)
  
  # Normalize each row to have unit norm
  coords_mat <- coords_mat / sqrt(rowSums(coords_mat^2))
  coords_mat[,1] = abs(coords_mat[,1])
  
  fdi_curr = Inf
  
  for(i in 1:iter){
    
    direc_new = coords_mat[i,]
    

    binned_dat <- SpheBinned_data(xin, yin, direc_new, M2)
    proj_binned <- binned_dat$binned_xmean %*% direc_new
    
    err <- 0
    for (l in 1:M2) {
      
      res <- SpheDirLocLin(
        xin, yin,
        direc_new, 
        proj_binned[l],
        bw2, ker = ker_gauss)
      err <- err + SpheGeoDist(res, binned_dat$binned_ymean[l, ])^2
      
    }
    fdi_new <- err / M2
    
    if (fdi_new < fdi_curr) {
      fdi_curr = fdi_new
      direc_curr = direc_new
      bw_curr = bw2
      M_curr = M2      
    }
    
    if(verbose){
      if(i %% 100 == 0){
        print(paste("Iteration number:", i,"/",iter))
      }
    }
  }
  return(list(est = normalize(direc_curr), bw = bw_curr, M = M_curr))
}


#### Directional local F-regression for spherical data, 
#### Direction along which to compute projection, and bandw choice
SpheDirLocLin <- function(xin, yin, direc, xout, bw, ker = ker_gauss) {
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## yin: n by m matrix of output (m: dimension of sphere; S^{m-1})
  ## direc: directional vector, length p 
  ## xout: vector of t
  ## bw: bandwidth b
  ## ker: ker_gauss, ker_unif, ker_epan
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  
  n <- nrow(xin)
  p <- ncol(xin)
  if (n <3) {
    stop("The number of observations is too small")
  }
  projec <- xin %*% direc 
  xin_eff <- projec
  
  mu0 <- mean(ker((xin_eff - xout) / bw))
  mu1 <- mean(ker((xin_eff - xout) / bw) * (xin_eff - xout))
  mu2 <- mean(ker((xin_eff - xout) / bw) * ((xin_eff - xout)^2))
  s <- ker((xin_eff - xout) / bw) * (mu2 - mu1 * (xin_eff - xout)) / (mu0 * mu2 - mu1^2)
  s <- as.vector(s) # vector of length n
  m <- ncol(yin)
  
  # Initial guess
  y0 <- colMeans(yin * s) 
  y0 <- y0 / l2norm(y0)
  
  # Check and adjust initial guess
  if (sum(sapply(1:n, function(i) {sum(yin[i,] * y0)})[which(ker((xout - xin_eff) / bw) > 0)] > 1 - 1e-8) != 0) {
    y0 <- y0 + rnorm(m) * 1e-3
    y0 <- y0 / l2norm(y0)
  }
  
  # Compute value
  objfun = function(y){
    y <- y / l2norm(y)
    if(all(sapply(y, l2norm) == 1)){
      f <- Inf
    } else {
      f <- mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
    }
    g <- 2 * colMeans(t(
      sapply(1:n, function(i) {
        SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y)})) * s)
    res = sapply(1:n, function(i){
      grad_i = SpheGeoGrad(yin[i,], y)
      return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
    }, simplify = "array")
    h = 2 * apply(res, 1:2, mean)
    return(list(value=f, gradient=g, hessian=h))
  }
  
  res <- trust(objfun, parinit = y0, rinit = 0.1, rmax = 1e5, minimize = T)
  return(normalize(res$argument))
}

#### Implements the selcetion of bandw for the local fr reg and, 
#### for the  optimal choice of bandw, selects the optimal bin size
SpheTuning <- function(xin, yin, direc, ker = ker_gauss){
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## yin: n by m matrix of output (m: dimension of sphere; S^{m-1})
  ## direc: directional vector, length p 
  ## ker: ker_gauss, ker_unif, ker_epan
  
  
  ## CV function to select the bandwidth
  bwCV <- function(bw, xin, yin, direc, ker = ker_gauss) {
    if (is.vector(xin)){
      stop("The number of observations is too small")
    } 
    
    n <- nrow(xin)
    #m <- ncol(yin)
    
    projec <- xin %*% direc
    ind_cv <- split(1:n, rep(1:5, length.out = n))
    cv_err <- 0
    for (i in 1:5) {
      xin_eff <- xin[-ind_cv[[i]], ]
      yin_eff <- yin[-ind_cv[[i]], ]
      
      for (k in 1:length(ind_cv[[i]])) {
        res <- SpheDirLocLin(xin_eff, yin_eff, direc, projec[ind_cv[[i]][k]], bw, ker)
        cv_err <- cv_err + SpheGeoDist(res, yin[ind_cv[[i]][k], ])^2
      }
      cv_err <- cv_err / length(ind_cv[[i]])
    }
    return(cv_err / 5)
  }
  
  ## CV function to select M (number of bins)
  bwCV_M <- function(xin, yin, direc, M, bw, ker = ker_gauss) {
    binned_dat <- SpheBinned_data(xin, yin, direc, M)
    xin_binned <- binned_dat$binned_xmean
    yin_binned <- binned_dat$binned_ymean
    proj_binned <- xin_binned %*% direc
    
    cv_err <- 0
    
    for (i in 1:M) {
      xin_eff <- xin_binned[-i, ]
      yin_eff <- yin_binned[-i, ]
      res <- SpheDirLocLin(xin_eff, yin_eff, direc, proj_binned[i], bw, ker = ker_gauss)
      cv_err <- cv_err + SpheGeoDist(res, yin_binned[i, ])^2
    }
    return(cv_err / M)
  }
  
  
  n <- nrow(xin)
  projec = xin %*% direc
  
  xinSt = unique(sort(projec))
  bw_min = max(c(diff(xinSt)))*1.1
  bw_max = (max(projec) - min(projec))/3
  if (bw_max < bw_min){
    if (bw_min > bw_max * 3/2){
      warning("Data is too sparse.")
      bw_max = bw_min * 1.01
    } else{
      bw_max = bw_max * 3/2
    }
  }
  
  ## bandwidth choice
  bw = optim(par = runif(1, min = bw_min, max = bw_max), #runif initial value
             fn = bwCV, xin = xin, yin = yin, direc = direc,
             method = "Brent", 
             lower = bw_min, upper = bw_max)$par
  
  ## M choice
  M_range = ceiling(n/c(2:30))
  M_range = unique(M_range[60 > M_range & M_range > 15])
  
  if (length(M_range) >0){
    
    cv_err_curr = Inf
    for(M in M_range){
      
      cv_err_new = bwCV_M(xin, yin, direc, M, bw)
      if (cv_err_new < cv_err_curr){
        
        cv_err_curr <- cv_err_new
        M_curr <- M
        
      }
      
    }
    
  } else{
    M = 15
  }
  
  return(c(bw, M))
}

### Binning step: given data and direction bins the support of the projection 
### and returns a representative point for the data (xin and yin)
#### Depends on the number of bins: M
SpheBinned_data <- function(xin, yin, direc, M) {
  if (M < 4){
    stop("The number of binned data should be greater than 3.")
  }
  
  n <- nrow(xin)
  p <- ncol(xin)
  m <- ncol(yin)
  
  if(n < M){
    stop("The number of binned data cannot exceed the number of observations.")
  }
  
  projec <- xin %*% direc
  range_of_projec <- seq(min(projec), max(projec), length.out = M)
  
  binned_xmean <- matrix(NA, M, p)
  binned_xmean[1, ] <- xin[which.min(projec), ]
  
  binned_ymean <- matrix(NA, M, m)
  binned_ymean[1, ] <- yin[which.min(projec), ]
  
  for (l in 2:(M - 1)) {
    idx = (n*l) %/% M
    idx_set = which(projec == sort(projec)[idx])
    binned_xmean[l, ] <- xin[idx_set[1], ]
    binned_ymean[l, ] <- yin[idx_set[1], ]
  }
  binned_xmean[M, ] <- xin[which.max(projec), ]
  binned_ymean[M, ] <- yin[which.max(projec), ]
  
  return(list(projec = projec, binned_xmean = binned_xmean, binned_ymean = binned_ymean))
}


#### Additional functions ####
l2norm <- function(x){
  sqrt(sum(x^2))
}

normalize <- function(x){
  x / sqrt(sum(x^2))
}

SpheGeoDist <- function(y1, y2){
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if(sum(y1 * y2) > 1){
    return(0)
  }else if(sum(y1 * y2) < -1){
    return(pi)
  }else{
    return(acos(sum(y1 * y2)))
  }
}

SpheGeoGrad <- function(y1, y2){
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  tmp = 1 - sum(y1 * y2)^2
  if (tmp <= 0) {
    return(- Inf * y1)
  } else{
    return(-(tmp)^(-1/2) + y1)
  }
}

SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  x = x / l2norm(x)
  y = y / l2norm(y)
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}

#### Kernels for local Fr reg
ker_gauss <- function(x) {
  return(exp(-x^2 / 2) / sqrt(2 * pi))
}

ker_unif <- function(x) {
  return(as.integer((x <= 1) & (x >= -1)))
}

ker_epan <- function(x, n = 1) {
  return((2 * n + 1) / (4 * n) * (1 - x^(2 * n)) * as.integer(abs(x) <= 1))
}


#### estimation performance calculator for estimates of the parameter, 
#### which are unit vectors
bias_calc2 <- function(est, true_beta, reps) {
  # true_beta: d vector
  # esti: reps*d matrix
  d <- length(true_beta)
  esti <- matrix(unlist(est), nrow = reps, byrow = T)
  extrinsic_mean <- colMeans(esti)
  
  if (sum(extrinsic_mean == 0) == d) {
    cat("Too few points\n")
    return(c(200, 200, 200))
  } else {
    intrinsic_mean <- as.vector(frechetMean(createM("Sphere"), t(esti)))
    bias_angle <- acos(sum(intrinsic_mean * true_beta))
    #dev <- mean(apply(esti, 1, function(x) sum(intrinsic_mean * x)))
    angle <- var(apply(esti, 1, function(x) acos(sum(intrinsic_mean * x))))
    #angle is dev?
    return(c(bias_angle, angle))
  }
}

#### Test ####
#### Sample data generation (Satarupa's version)
SpheGenerate_data <- function(n, sd, true_beta, link){
  d <- length(true_beta)
  xin <- matrix(runif(n * d), nrow = n)
  xin <- t(apply(xin, 1, normalize))
  link_proj = link(sapply(1:n, function(i){sum(true_beta * xin[i,])}))
  link_proj = (link_proj - min(link_proj))/(max(link_proj) - min(link_proj))
  err_sd = sd
  phi_true = acos(link_proj)
  theta_true = pi * link_proj
  m = 3
  ytrue = cbind(sin(phi_true) * cos(theta_true),
                sin(phi_true) * sin(theta_true),
                cos(phi_true))
  basis = list(
    b1 = cbind(cos(phi_true) * cos(theta_true),
               cos(phi_true) * sin(theta_true),
               -sin(phi_true)),
    b2 = cbind(sin(theta_true), -cos(theta_true), 0)
  )
  yin_tg = basis$b1 * rnorm(n, mean = 0, sd = err_sd) +
    basis$b2 * rnorm(n, mean = 0, sd = err_sd)
  yin = matrix(0, nrow = n, ncol = m)
  for(i in 1:n){
    tgNorm = sqrt(sum(yin_tg[i]^2))
    if(tgNorm < 1e-10){
      yin[i,] <- ytrue[i,]
    } else {
      yin[i,] <- sin(tgNorm) * yin_tg[i] / tgNorm + cos(tgNorm) * ytrue[i,]
    }
  }
  return(list(xin = xin, yin = yin))
}


#set.seed(100)
#b <- c(3, -1.3, -3, 1.7)
#b0 <- normalize(b)
#b0 #0.6313342 -0.2735781 -0.6313342  0.3577560


#dat <- SpheGenerate_data(100, 0, b0, function(x) x)
#res_sphe <- SIdxSpheReg(xin = dat$xin, yin = dat$yin)
#res_sphe

