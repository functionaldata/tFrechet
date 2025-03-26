library(Matrix)
library(osqp)


# Main Function
SIdxDenReg = function(xin, qin, bw=NULL, M=NULL, ker = ker_gauss, lower = -Inf, upper = Inf,
                      iter = 1000, verbose = T) {
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## qin: n by m matrix of quantile function  (m: length of quantile)
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  ## iter: used for generation of directions.
  ## verbose: print the iteration counts?
  
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  if (!is.matrix(xin)){
    stop("xin should be matrix.")
  }
  if(!is.matrix(qin)){
    stop("qin should be matrix.")
  }
  
  p <- ncol(xin)
  
  ## Parameter (bandwidth, bin size) choice using cross-validation
  needParam <- (is.null(M) | is.null(bw))
  
  if (needParam) {
    param <- DenTuning(xin, qin, normalize(rep(1, p)))
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
  
  ## Find the single index
  fdi_curr = Inf
  for(i in 1:nrow(coords_mat)){
    
    direc_new = coords_mat[i,]
    
    binned_dat <- DenBinned_data(xin, qin, direc_new, M2)
    proj_binned <- binned_dat$binned_xmean %*% direc_new
    
    err <- 0
    for (l in 1:M2) {
      
      res <- DenDirLocLin(
        xin, qin,
        direc_new, 
        proj_binned[l],
        bw2, ker = ker_gauss, lower = -Inf, upper = Inf)
      err <- err + mean((sort(res) - sort(binned_dat$binned_qmean[l, ]))^2)
    
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
        print(paste("Iteration number:", i,"/",nrow(coords_mat)))
      }
    }
    
  }
  
  return(list(est = normalize(direc_curr), bw = bw_curr, M = M_curr))
  
}

#### Directional local F-regression given probability density, 
#### Direction along which to compute projection, and bandw choice
DenDirLocLin <- function(xin, qin, direc, xout, bw, ker = ker_gauss, 
                      lower = -Inf, upper = Inf) {
  
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## qin: n by m matrix of quantile function  (m: length of quantile)
  ## direc: directional vector, length p
  ## xout: A k by p matrix with output measurements of the predictors. Default is \code{xin}.
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  
  n <- nrow(xin)

  if (n < 3) {
    stop("The number of observations is too small")
  }
  projec <- xin %*% direc
  xin_eff <- projec

  mu0 <- mean(ker((xin_eff - xout) / bw))
  mu1 <- mean(ker((xin_eff - xout) / bw) * (xin_eff - xout))
  mu2 <- mean(ker((xin_eff - xout) / bw) * ((xin_eff - xout)^2))
  s <- ker((xin_eff - xout) / bw) * (mu2 - mu1 * (xin_eff - xout)) / (mu0 * mu2 - mu1^2)
  s <- as.vector(s)
  m <- ncol(qin)
  b0 <- c(lower, rep(0, m - 1), -upper)
  
  Pmat <- Diagonal(n = m)
  Amat <- bandSparse(m+1, m, k = c(0, -1), diag = list(rep(1, m), rep(-1, m)))
  
  gx <- colMeans(qin * s)
  
  prob <- osqp(Pmat, -gx, Amat, b0, pars = osqpSettings(verbose = FALSE))
  results <- prob$Solve()$x
  return(results)
}


#### Implements the selcetion of bandw for the local F-reg and, 
#### for the optimal choice of bandw, selects the optimal bin size
DenTuning <- function(xin, qin, direc, ker = ker_gauss){
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## qin: n by m matrix of quantile function  (m: length of quantile)
  ## direc: directional vector, length p
  ## ker: ker_gauss, ker_unif, ker_epan
  
  ## CV function to select the bandwidth
  bwCV <- function(bw, xin, qin, direc, ker = ker_gauss, lower = -Inf, upper = Inf) {
    
    if (is.vector(xin)){
      stop("The number of observations is too small")
    } 
    
    n <- nrow(xin)
    p <- ncol(xin)
    m <- ncol(qin)
    
    projec <- xin %*% direc
    ind_cv <- split(1:n, rep(1:5, length.out = n))
    cv_err <- 0
    
    for (i in 1:5) {
      xin_eff <- xin[-ind_cv[[i]], ]
      qin_eff <- qin[-ind_cv[[i]], ]
      
      for (k in 1:length(ind_cv[[i]])) {
        res <- DenDirLocLin(xin_eff, qin_eff , direc, projec[ind_cv[[i]][k]], bw, ker, lower = -Inf, upper = Inf)
        cv_err <- cv_err + mean((sort(res) - sort(qin[ind_cv[[i]][k], ]))^2)
      }
      cv_err <- cv_err / length(ind_cv[[i]])
    }
    return(cv_err / 5)
  }
  
  ## CV function to select M
  bwCV_M <- function(xin, qin, direc, M, bw, ker = ker_gauss, lower = -Inf, upper = Inf) {
    binned_dat <- DenBinned_data(xin, qin, direc, M)
    xin_binned <- binned_dat$binned_xmean
    qin_binned <- binned_dat$binned_qmean
    proj_binned <- xin_binned %*% direc
    
    m <- ncol(qin_binned)
    
    cv_err <- 0
    
    for (i in 1:M) {
      xin_eff <- xin_binned[-i, ]
      qin_eff <- qin_binned[-i, ]
      dat_cv <- list(xin = xin_eff, qin = qin_eff)
      res <- DenDirLocLin(xin_eff, qin_eff, direc, proj_binned[i], bw, ker = ker_gauss, lower = -Inf, upper = Inf)
      cv_err <- cv_err + mean((sort(res) - sort(qin_binned[i, ]))^2)
    }
    
    if(!is.nan(cv_err)){
      return(cv_err / M)
    } else{
      return(Inf)
    }
    
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
  
  ## bandwidth choice using bwCV
  bw = optim(par = runif(1, min = bw_min, max = bw_max), 
             fn = bwCV, xin = xin, qin = qin, direc = direc,
             method = "Brent", 
             lower = bw_min, upper = bw_max)$par
  
  M_range = ceiling(n/c(2:30))
  M_range = unique(M_range[60 > M_range & M_range > 15])
  
  if (length(M_range) >0){
    
    cv_err_curr = Inf
    for(M in M_range){
      
      cv_err_new = bwCV_M(xin,qin, direc, M, bw)
      if (cv_err_new < cv_err_curr){
      
        cv_err_curr <- cv_err_new
        M_curr <- M
      
      }
      
    }
    
  } else{
    M = 15
  }
  
  #end
  return(c(bw, M))
  
}

#### Binning step: given data and direction bins the support of the projection
#### Returns a representative point for the data (xin and qin)
DenBinned_data <- function(xin, qin, direc, M) {
  
  if (M < 4){
    stop("The number of binned data should be greater than 3.")
  }
  
  n <- nrow(xin)
  p <- ncol(xin)
  m <- ncol(qin)
  
  if(n < M){
    stop("The number of binned data cannot exceed the number of observations.")
  }
  
  projec <- xin %*% direc
  range_of_projec <- seq(min(projec), max(projec), length.out = M)
  
  binned_xmean <- matrix(NA, M, p)
  binned_xmean[1, ] <- xin[which.min(projec), ]
  
  binned_qmean <- matrix(NA, M, m)
  binned_qmean[1, ] <- qin[which.min(projec), ]
  
  for (l in 2:(M - 1)) {
    idx = (n*l)%/%M
    idx_set = which(projec == sort(projec)[idx])
    binned_xmean[l, ] <- xin[idx_set[1], ]
    
    binned_qmean[l, ] <- qin[idx_set[1], ]
  }
  
  binned_xmean[M, ] <- xin[which.max(projec), ]
  binned_qmean[M, ] <- qin[which.max(projec), ]
  
  return(list(projec = projec, binned_xmean = binned_xmean, binned_qmean = binned_qmean))
}


#### Additional functions ####
ker_gauss <- function(x) {
  return(exp(-x^2 / 2) / sqrt(2 * pi))
}

ker_unif <- function(x) {
  return(as.integer((x <= 1) & (x >= -1)))
}

ker_epan <- function(x, n = 1) {
  return((2 * n + 1) / (4 * n) * (1 - x^(2 * n)) * as.integer(abs(x) <= 1))
}

normalize <- function(x){
  x / sqrt(sum(x^2))
}

DenGen_data_setting = function(n, true_beta, link){
  
  d = length(true_beta)
  rho = 1/4
  xin = 2*pnorm(MASS::mvrnorm(n, mu = rep(0,d), Sigma = (1-rho)*diag(d) + matrix(rho,d,d) ))-1
  qSup = qbeta(seq(0.01, 0.99, by = 0.01), 1/2, 1/2)
  m = length(qSup) + 2
  
  qin = matrix(0, n, m)
  for(i in 1:n){
    
    proj = sum(true_beta * xin[i, ])
    qin[i,] = qnorm(c(1e-6, qSup, 1 - 1e-6), 
                    mean = link(proj) ,
                    sd = 1)
                    #sd=exp(proj)/(1+exp(proj)))
    
  }
  
  return(list(xin = xin, qin = qin))
  
}
