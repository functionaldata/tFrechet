library(frechet)
library(Matrix)

# Main Function: Single Index F-regression with network response using Frobenius metric
SIdxNetReg <- function(xin, Min, bw = NULL, M = NULL, ker = ker_gauss, iter = 1000, verbose = TRUE) {
  ## xin: n by p matrix of input (n: number of inputs, p: dimension of predictors)
  ## Min: m by m by n 3-d array of graph Laplacian  (m: number of nodes)
  ## bw: bandwidth b
  ## M: size of binning
  ## ker: ker_gauss, ker_unif, ker_epan
  ## iter: used for generation of directions.
  ## verbose: print the iteration counts?
  
  # Input validation
  if (!is.matrix(xin)) stop("xin should be a matrix.")
  if (!is.array(Min) || length(dim(Min)) != 3) stop("Min should be a 3-dimensional array.")
  
  p <- ncol(xin)
  
  ## Parameter (bandwidth, bin size) choice using cross-validation
  needParam <- (is.null(M) | is.null(bw))
  
  if (needParam) {
    param <- NetTuning(xin, Min, normalize(rep(1,p)))
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
  
  fdi_curr <- Inf
  
  for (i in 1:iter) {
    
    direc_new = coords_mat[i,]
    
    binned_dat <- NetBinned_data(xin, Min, direc_new, M2)
    proj_binned <- binned_dat$binned_xmean %*% direc_new
    
    res <- NetDirLocLin(xin, Min, direc_new, proj_binned, bw2)
    
    err <- 0
    for (l in 1:M2) {
      err <- err + sum((res$predict[[l]] - binned_dat$binned_Mmean[,,l])^2)
    }
    
    fdi_new <- err / M2
    if (fdi_new < fdi_curr) {
      fdi_curr <- fdi_new
      direc_curr <- direc_new
      bw_curr <- bw2
      M_curr <- M2
    }
    
    if (verbose && i %% 100 == 0) {
      print(paste("Iteration number:", i, "/", iter))
    }
  }
  
  return(list(est = normalize(direc_curr), bw = bw_curr, M = M_curr))
}

# Directional local F-regression given network response (Frobenius metric)
NetDirLocLin <- function(xin, Min, direc, xout, bw) {
  
  # Input validation
  if (!is.matrix(xin)) stop("xin should be a matrix.")
  if (!is.array(Min) || length(dim(Min)) != 3) stop("Min should be a 3-dimensional array.")
  
  projec <- xin %*% direc
  
  Min_list = list()
  for(i in 1:dim(Min)[3]){
    
    Min_list[[i]] = Min[,,i]
    
  }
  
  lnr_fit = lnr(gl = Min_list, x = projec, xOut = xout, optns = list(bw = bw))
  
  return(lnr_fit)
}


# Bandwidth tuning function using cross-validation
NetTuning <- function(xin, Min, direc, ker = ker_gauss) {
  
  bwCV <- function(bw, xin, Min, direc, ker = ker_gauss) {
    projec <- xin %*% direc
    ind_cv <- split(1:nrow(xin), rep(1:5, length.out = nrow(xin)))
    cv_err <- 0
    
    for (i in 1:5) {
      xin_eff <- xin[-ind_cv[[i]], ]
      Min_eff <- Min[,,-ind_cv[[i]]]
      res <- NetDirLocLin(xin_eff, Min_eff, direc, projec, bw)
      
      for (k in ind_cv[[i]]) {
        cv_err <- cv_err + sum((res$predict[[k]] - Min[,,k])^2)
      }
      cv_err <- cv_err / length(ind_cv[[i]])
    }
    return(cv_err / 5)
  }
  
  bwCV_M <- function(xin, Min, direc, M, bw, ker = ker_gauss) {
    
    # Bin the data along the direction vector
    binned_dat <- NetBinned_data(xin, Min, direc, M)
    xin_binned <- binned_dat$binned_xmean
    Min_binned <- binned_dat$binned_Mmean
    proj_binned <- xin_binned %*% direc
    
    cv_err <- 0
    
    for (i in 1:M) {
      xin_eff <- xin_binned[-i, , drop = FALSE]
      Min_eff <- Min_binned[,,-i, drop = FALSE]
      
      res <- NetDirLocLin(xin_eff, Min_eff, direc, proj_binned[i], bw)
      cv_err <- cv_err + sum((res$predict[[1]] - Min_binned[,,i])^2)
    }
    
    # Return error; avoid returning NaN or Inf in case of numerical instability
    if (!is.nan(cv_err)) {
      return(cv_err / M)
    } else {
      return(Inf)
    }
  }
  
  n <- nrow(xin)

  bw_min <- max(diff(sort(xin %*% direc))) * 1.1
  bw_max <- (max(xin %*% direc) - min(xin %*% direc)) / 3
  
  if (bw_max < bw_min) {
    bw_max <- max(bw_min * 1.01, bw_min * 3/2)
  }
  
  bw <- optim(par = runif(1, bw_min, bw_max), fn = bwCV, xin = xin, Min = Min, direc = direc,
              method = "Brent", lower = bw_min, upper = bw_max)$par
  
  M_range = ceiling(n/c(2:30))
  M_range = unique(M_range[60 > M_range & M_range > 15])
  
  if (length(M_range) >0){
    
    cv_err_curr = Inf
    for(M in M_range){
      
      cv_err_new = bwCV_M(xin, Min, direc, M, bw, ker)
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

# Binning function
NetBinned_data <- function(xin, Min, direc, M) {
  
  if (M < 4) stop("The number of binned data should be greater than 3.")
  if (nrow(xin) < M) stop("The number of binned data cannot exceed the number of observations.")
  
  p = ncol(xin)
  
  projec <- xin %*% direc
  range_of_projec <- seq(min(projec), max(projec), length.out = M)
  binned_xmean <- matrix(NA, M, p)
  binned_xmean[1, ] <- xin[which.min(projec), ]
  
  binned_Mmean <- array(NA, dim = c(dim(Min)[1], dim(Min)[2], M))
  binned_Mmean[,, 1] <- Min[,, which.min(projec)]
  
  for (l in 2:(M-1)) {
    idx <- (nrow(xin) * l) %/% M
    idx_set <- which(projec == sort(projec)[idx])
    binned_xmean[l, ] <- xin[idx_set[1], ]
    binned_Mmean[,, l] <- Min[,, idx_set[1]]
  }
  
  binned_xmean[M, ] <- xin[which.max(projec), ]
  binned_Mmean[,,M] <- Min[,,which.max(projec)]
  
  return(list(binned_xmean = binned_xmean, binned_Mmean = binned_Mmean))
}

# Gaussian kernel
ker_gauss <- function(x) {
  return(exp(-x^2 / 2) / sqrt(2 * pi))
}

# Vector normalization function
normalize <- function(x) {
  return(x / sqrt(sum(x^2)))
}

lnr <- function(gl = NULL, x = NULL, xOut = NULL, optns = list()) {
  if (is.null(gl) | is.null(x)) {
    stop("requires the input of both gl and x")
  }
  if (is.null(optns$metric)) {
    optns$metric <- "frobenius"
  }
  if (!(optns$metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    optns$alpha <- 1
  }
  if (optns$alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$kernel)) {
    optns$kernel <- "gauss"
  }
  if (is.null(optns$bw)) {
    optns$bw <- NA
  }
  if (is.null(optns$digits)) {
    optns$digits <- NA
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (p > 2) {
    stop("local method is designed to work in low dimensional case (p is either 1 or 2)")
  }
  if (!is.na(sum(optns$bw))) {
    if (sum(optns$bw <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw) != p) {
      stop("dimension of bandwidth does not agree with x")
    }
  }
  if (!is.list(gl)) {
    if (is.array(gl)) {
      gl <- lapply(seq(dim(gl)[3]), function(i) gl[, , i])
    } else {
      stop("gl must be a list or an array")
    }
  }
  if (length(gl) != n) {
    stop("the number of rows in x must be the same as the number of graph Laplacians in gl")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }
  nodes <- colnames(gl[[1]])
  m <- ncol(gl[[1]]) # number of nodes
  if (is.null(nodes)) nodes <- 1:m
  glVec <- matrix(unlist(gl), ncol = m^2, byrow = TRUE) # n by m^2
  if (substr(optns$metric, 1, 1) == "p") {
    glAlpha <- lapply(gl, function(gli) {
      eigenDecom <- eigen(gli)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^optns$alpha) %*% t(U)
    })
    glAlphaVec <- matrix(unlist(glAlpha), ncol = m^2, byrow = TRUE) # n by m^2
  }
  
  # initialization of OSQP solver
  W <- 2^32 # bound on the weights of edges in the graph
  nConsts <- m^2 # number of constraints
  l <- c(rep.int(0, m * (m + 1) / 2), rep.int(-W, m * (m - 1) / 2))
  u <- rep.int(0, nConsts)
  q <- rep.int(0, m^2)
  P <- diag(m^2)
  consts <- matrix(0, nrow = nConsts, ncol = m^2)
  k <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
      consts[k, (i - 1) * m + j] <- -1
    }
  }
  for (i in 1:m) {
    consts[k + i, ((i - 1) * m + 1):(i * m)] <- rep(1, m)
  }
  k <- k + m
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
    }
  }
  model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))
  
  # select kernel
  Kern <- kerFctn(optns$kernel)
  K <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * Kern(x[i] / h[i])
    }
    return(as.numeric(k))
  }
  
  # choose bandwidth by cross-validation
  if (is.na(sum(optns$bw))) {
    hs <- matrix(0, p, 20)
    for (l in 1:p) {
      hs[l, ] <- exp(seq(
        from = log(n^(-1 / (1 + p)) * (max(x[, l]) - min(x[, l])) / 10),
        to = log(5 * n^(-1 / (1 + p)) * (max(x[, l]) - min(x[, l]))),
        length.out = 20
      ))
    }
    cv <- array(0, 20^p)
    for (k in 0:(20^p - 1)) {
      h <- array(0, p)
      for (l in 1:p) {
        kl <- floor((k %% (20^l)) / (20^(l - 1))) + 1
        h[l] <- hs[l, kl]
      }
      for (j in 1:n) {
        a <- x[j, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * (xi - a)))
          mu2 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * ((xi - a) %*% t(xi - a))))
        }
        skip <- FALSE
        tryCatch(solve(mu2), error = function(e) skip <<- TRUE)
        if (skip) {
          cv[k + 1] <- Inf
          break
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(as.matrix(x[-j, ]), 1, function(xi) {
          K(xi - a, h) * (1 - wc %*% (xi - a))
        }) # weight
        if (substr(optns$metric, 1, 1) == "f") {
          qNew <- apply(glVec[-j, ], 2, weighted.mean, w) # m^2
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits) # round
          fitj[fitj > 0] <- 0 # off diagonal should be negative
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
        } else if (substr(optns$metric, 1, 1) == "p") {
          bAlpha <- matrix(apply(glAlphaVec[-j, ], 2, weighted.mean, w), ncol = m) # m by m
          eigenDecom <- eigen(bAlpha)
          Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
          U <- eigenDecom$vectors
          qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits) # round
          fitj[fitj > 0] <- 0 # off diagonal should be negative
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
          # eigenDecom <- eigen(fitj)
          # Lambda <- pmax(Re(eigenDecom$values), 0)
          # U <- eigenDecom$vectors
          # fitjAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
          # cv[k+1] <- cv[k+1] + sum((glAlpha[[j]]-fitjAlpha)^2)/n# using Euclidean power metric
        }
      }
    }
    bwi <- which.min(cv)
    optns$bw <- array(0, p)
    for (l in 1:p) {
      kl <- floor((bwi %% (20^l)) / (20^(l - 1))) + 1
      optns$bw[l] <- hs[l, kl]
    }
  }
  
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  if (substr(optns$metric, 1, 1) == "f") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
      }) # weight
      qNew <- apply(glVec, 2, weighted.mean, w) # m^2
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
    }
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # weight
        qNew <- apply(glVec, 2, weighted.mean, w) # m^2
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  } else if (substr(optns$metric, 1, 1) == "p") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
      }) # weight
      bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m by m
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
      # eigenDecom <- eigen(fit[[i]])
      # Lambda <- pmax(Re(eigenDecom$values), 0)
      # U <- eigenDecom$vectors
      # fitiAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
      # residuals[i] <- sqrt(sum((glAlpha[[i]]-fitiAlpha)^2))# using Euclidean power metric
    }
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # weight
        bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  }
  class(res) <- "nr"
  res
}

kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}



#set.seed(1)
#n = 100
#m=3

# Set Parameters
#X = matrix(c(runif(n, -1, 1), 
#             runif(n, -1, 1),
#             runif(n, 1, 2),
             
#             rgamma(n, 3, 1),
#             rgamma(n, 4, 1),
#             rgamma(n, 5, 1),

#             rbinom(n, 1, 0.2),
#             rbinom(n, 1, 0.3),
#             rbinom(n, 1, 0.5)), n) 

#y = lapply(1:n, function(i){
#  a = 2 * sin(pi*X[i,1])^2*X[i,7] + cos(pi*X[i,2])^2*(1-X[i,7])
#  b = X[i,4]*X[i,8]+X[i,5]*(1-X[i,8])
  
#  Vec = -rbeta(m*(m-1)/2, shape1 = a, shape2 = b)
#  temp = matrix(0, nrow = m, ncol = m)
#  temp[lower.tri(temp)] = Vec
#  temp <- temp + t(temp)
#  diag(temp) = -colSums(temp)
#  return(temp)
#})

#y_mat = array(0, c(m, m, n))
#for(j in 1:n){
#  y_mat[,,j] = y[[j]]
#}

#res_net = SIdxNetReg(xin = X, Min = y_mat)
