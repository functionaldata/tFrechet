library(frechet)
library(Matrix)

# Main Function : Single Index F-regression with covariance response with Frobenius metric
SIdxCovReg = function(xin, Min, bw=NULL, M=NULL, ker = ker_gauss, lower = -Inf, upper = Inf, iter =  500,
                      verbose = T){
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## Min: q by q by n array where \code{M[,,i]} contains the i-th covariance matrix of dimension q by q
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  ## iter: generation of directions.  
  ## verbose: print the iteration counts?
  
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  if (!is.matrix(xin)){
    stop("xin should be matrix.")
  }
  if(!is.array(Min)){
    stop("Min should be array.")
    
    if(length(dim(Min)) != 3){
      stop("Min should be 3-dimensional array.")
    }
    
  }
  
  p <- ncol(xin)
  
  fdi_curr = Inf
  
  for(i in 1:iter){
    
    #set.seed(i)
    direc_new = normalize(rnorm(n = p))
    
    if(direc_new[1] < 0){
      direc_new = -1 * direc_new
    }
    
    ## Parameter (bandwidth, bin size) choice using cross-validation
    if(is.null(bw) | is.null(M)){
      
      if(is.null(bw)){
        
        param = CovTuning(xin, Min, direc_new)
        bw2 = param[1]
        
        if(!is.null(M)){
          
          M2 = M
          if (M < 4){stop("The number of binned data should be greater than 3")}
          
        } else{
          
          M2 = param[2]
          
        }
        
      }
      
    }else {
      bw2 = bw
      M2 = M
    }
    
    binned_dat <- CovBinned_data(xin, Min, direc_new, M2)
    proj_binned <- binned_dat$binned_xmean %*% direc_new
    
    err <- 0
    for (l in 1:M2) {
      
      res <- CovDirLocLin(
        xin, Min,
        direc_new, 
        proj_binned[l],
        bw2, ker = ker_gauss, lower = -Inf, upper = Inf)
      err <- err + sum((res - binned_dat$binned_Mmean[,,l])^2)
      
    }
    fdi_new <- err / M2
    
    if (fdi_new < fdi_curr) {
      
      fdi_curr = fdi_new
      direc_curr = direc_new
      bw_curr = bw2
      M_curr = M2
      
    }
    
    if(verbose){
      if(i %% 10 == 0){
        print(paste("Iteration number:", i,"/",iter))
      }
    }
  }
  
  return(list(est = normalize(direc_curr), bw = bw_curr, M = M_curr))
  
}


#### Directional local F-regression given covariance response with Frobenius metric, 
#### Direction along which to compute projection, and bandw choice
CovDirLocLin <- function(xin, Min, direc, xout, bw, ker = ker_gauss, 
                         lower = -Inf, upper = Inf) {
  
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## Min: q by q by n array where \code{M[,,i]} contains the i-th covariance matrix of dimension q by q
  ## direc: directional vector, length p
  ## xout: A k by p matrix with output measurements of the predictors. Default is \code{xin}.
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  
  if (is.vector(xin)){
    stop("The number of observations is too small")
  }
  
  if(!is.array(Min)){
    stop("Min should be array.")
    
    if(length(dim(Min)) != 3){
      stop("Min should be 3-dimensional array.")
    }
    
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
  q = dim(Min)[1]
  
  sL = sum(s)
  
  M_res = matrix(0, nrow = q, ncol = q)
  for(i in 1:n){
    M_res = M_res + s[i]*Min[,,i]/sL
  }
  
  M_res = as.matrix(Matrix::nearPD(M_res)$mat)
  M_res = Matrix::forceSymmetric(M_res)

  return(M_res)
}

#### Implements the selcetion of bandw for the local F-reg and, 
#### for the optimal choice of bandw, selects the optimal bin size
CovTuning <- function(xin, Min, direc, ker = ker_gauss){
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## Min: q by q by n array where \code{M[,,i]} contains the i-th covariance matrix of dimension q by q
  ## direc: directional vector, length p
  ## ker: ker_gauss, ker_unif, ker_epan
  
  ## CV function to select the bandwidth
  bwCV <- function(bw, xin, Min, direc, ker = ker_gauss, lower = -Inf, upper = Inf) {
    
    if (is.vector(xin)){
      stop("The number of observations is too small")
    } 
    
    n <- nrow(xin)
    p <- ncol(xin)
    q <- dim(Min)[1]
    
    projec <- xin %*% direc
    ind_cv <- split(1:n, rep(1:5, length.out = n))
    cv_err <- 0
    
    for (i in 1:5) {
      xin_eff <- xin[-ind_cv[[i]], ]
      Min_eff <- Min[,,-ind_cv[[i]]]
      
      for (k in 1:length(ind_cv[[i]])) {
        res <- CovDirLocLin(xin_eff, Min_eff , direc, projec[ind_cv[[i]][k]], bw, ker, lower = -Inf, upper = Inf)
        cv_err <- cv_err + sum((res - Min[,,ind_cv[[i]][k]])^2) 
        
      }
      cv_err <- cv_err / length(ind_cv[[i]])
    }
    return(cv_err / 5)
  }
  
  ## CV function to select M
  bwCV_M <- function(xin, Min, direc, M, bw, ker = ker_gauss, lower = -Inf, upper = Inf) {
    binned_dat <- CovBinned_data(xin, Min, direc, M)
    xin_binned <- binned_dat$binned_xmean
    Min_binned <- binned_dat$binned_Mmean
    proj_binned <- xin_binned %*% direc
    
    cv_err <- 0
    
    for (i in 1:M) {
      xin_eff <- xin_binned[-i, ]
      Min_eff <- Min_binned[,,-i]
      res <- CovDirLocLin(xin_eff, Min_eff, direc, proj_binned[i], bw, ker = ker_gauss, lower = -Inf, upper = Inf)
      cv_err <- cv_err + sum((res - Min_binned[,,i])^2)
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
             fn = bwCV, xin = xin, Min = Min, direc = direc,
             method = "Brent", 
             lower = bw_min, upper = bw_max)$par
  
  
  M_range = ceiling(n^(1/c(2:7)))
  M_range = unique(M_range[M_range > 3])
  if (length(M_range) >0){
    
    cv_err_curr = Inf
    for(M in M_range){
      
      cv_err_new = bwCV_M(xin,Min, direc, M, bw)

      if (cv_err_new < cv_err_curr){
        
        cv_err_curr <- cv_err_new
        M_curr <- M
        
      }
      
    }
    
  } else{
    M = n
  }
  
  #end
  return(c(bw, M))
  
}


#### Binning step: given data and direction bins the support of the projection
#### Returns a representative point for the data (xin and Min)
CovBinned_data <- function(xin, Min, direc, M) {
  
  if (M < 4){
    stop("The number of binned data should be greater than 3.")
  }
  
  n <- nrow(xin)
  p <- ncol(xin)
  q <- dim(Min)[1]
  
  if(n < M){
    stop("The number of binned data cannot exceed the number of observations.")
  }
  
  projec <- xin %*% direc
  range_of_projec <- seq(min(projec), max(projec), length.out = M)
  
  binned_xmean <- matrix(NA, M, p)
  binned_xmean[1, ] <- xin[which.min(projec), ]
  
  binned_Mmean <- array(NA, dim=c(q,q,M))
  binned_Mmean[,, 1] <- Min[,, which.min(projec)]
  
  for (l in 2:(M - 1)) {
    idx = (n*l)%/%M
    idx_set = which(projec == sort(projec)[idx])
    binned_xmean[l, ] <- xin[idx_set[1], ]
    
    binned_Mmean[,, l] <- Min[,,idx_set[1]]
  }
  
  binned_xmean[M, ] <- xin[which.max(projec), ]
  binned_Mmean[,,M] <- Min[,,which.max(projec)]
  
  return(list(projec = projec, binned_xmean = binned_xmean, binned_Mmean = binned_Mmean))
}

#### Additional functions ####
ker_gauss <- function(x) {
  return(exp(-x^2 / 2) / sqrt(2 * pi))
}

normalize <- function(x){
  x / sqrt(sum(x^2))
}

CovGen_data_setting = function(n, true_beta, link){
  
  d = length(true_beta)
  rho = 1/4
  xin = 2*pnorm(MASS::mvrnorm(n, mu = rep(0,d), Sigma = (1-rho)*diag(d) + matrix(rho,d,d) ))-1
  
  q = 3
  Min = array(0, c(q,q,n))
  for(i in 1:n){
    
    proj = sum(true_beta * xin[i, ])
    eig = c(exp(link(proj)),exp(link(proj)/2), exp(-link(proj)))
    Min[,,i] = diag(eig) 
    
  }
  
  return(list(xin = xin, Min = Min))
  
}

#### Test

set.seed(100)
b <- c(3, -1.3, -3, 1.7)
b0 <- normalize(b)
b0 #0.6313342 -0.2735781 -0.6313342  0.3577560

for(rep in 1:10){
  dat <- CovGen_data_setting(100, b0, function(x) x)
  res_cov <- SIdxCovReg(dat$xin, dat$Min, iter = 100)
  print(res_cov$est)
}

save(res_cov, file = "res_cov.RData")

