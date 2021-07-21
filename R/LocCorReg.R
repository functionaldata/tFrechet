#'@title Local Fr$\'{e}$chet regression for correlation matrices
#'@description Local Fr$\'{e}$chet regression for correlation matrices with Euclidean predictors.
#'@param x an n by p matrix of predictors.
#'@param y a list (length n) of $m$ by $m$ correlation matrix.
#'@param xOut an nOut by p matrix of output predictor levels.
#'@param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'@details Available control options are
#'\describe{
#'\item{metric}{choice of metric. 'frobenius' and 'power' are supported, which corresponds to Frobenius metric and Euclidean power metric, repectively. Default is Frobenius metric}
#'\item{alpha}{the power for Euclidean power metric. Default is $1$ which corresponds to Frobenius metric.}
#'\item{kernel}{choice of kernel. 'gauss', 'uniform', and 'epan' are supported, corresponding to Gaussian kernel, uniform kernel, and Epanechnikov kernel, respectively. Default is 'gauss'.}
#'\item{bw}{bandwidth for local Fr\'{e}chet regression, if not entered it would be chosen from cross validation.}
#'\item{digits}{the integer indicating the number of decimal places (round) to be kept in the output. Default is NULL, which means no round operation.}
#'}
#'@return A \code{corReg} object --- a list containing the follwing fields:
#'\item{fit}{a list of estimated correlation matrices at \code{x}.}
#'\item{predict}{a list of estimated correlation matrices at \code{xOut}. Included if \code{xOut} is not \code{NULL}.}
#'\item{residuals}{Frobenius distance between the true and fitted correlation matrices.}
#'\item{y}{the response used.}
#'\item{x}{the predictor used.}
#'\item{xOut}{the output predictor level used.}
#'\item{optns}{the control options used.}
#'@examples
#'# Generate simulation data
#'n <- 100
#'m <- 10
#'d <- m*(m-1)/2
#'xOut <- seq(0.1, 0.9, length.out = 9)
#'x <- runif(n, min = 0, max = 1)
#'y <- list()
#'for(i in 1:n){
#'  yVec <- rbeta(d, shape1 = sin(pi*x[i]), shape2 = 1-sin(pi*x[i]))
#'  y[[i]] <- matrix(0, nrow = m, ncol = m)
#'  y[[i]][lower.tri(y[[i]])] <- yVec
#'  y[[i]] <- y[[i]] + t(y[[i]])
#'  diag(y[[i]]) <- 1
#'}
#'# Frobenius metric
#'fit1 <- LocCorReg(x, y, xOut, 
#'                  optns = list(metric = 'frobenius', digits = 2))
#'# Euclidean power metric
#'fit2 <- LocCorReg(x, y, xOut, 
#'                  optns = list(metric = 'power', alpha = .5, 
#'                               kernel = 'epan', bw = 0.08))
#'@references
#'\itemize{
#'\item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#'}
#'@export

LocCorReg <- function(x, y, xOut = NULL, optns = list()){
  if(is.null(optns$metric)){
    metric <- "frobenius"
  } else {
    metric <- optns$metric
  }
  if(!(metric %in% c("frobenius", "power"))){
    stop("metric choice not supported.")
  }
  if(is.null(optns$kernel)){
    kernel <- 'gauss'
  } else{
    kernel <- optns$kernel
  }
  if(is.null(optns$bw)){
    bw <- NA
  } else{
    bw <- optns$bw
  }
  if(is.null(optns$alpha)){
    alpha <- 1
  } else{
    alpha <- optns$alpha
  }
  if(alpha<0){
    stop('alpha must be non-negative')
  }
  if(is.null(optns$digits)){
    digits <- NA
  } else {
    digits <- optns$digits
  }
  if(!is.matrix(x)){
    if(is.data.frame(x) | is.vector(x)) x <- as.matrix(x)
    else stop('x must be a matrix or a data frame')
  }
  if(!is.null(xOut)){
    if(!is.matrix(xOut)){
      if(is.data.frame(xOut) | is.vector(xOut)) xOut <- as.matrix(xOut)
      else stop('xOut must be a matrix or a data frame')
    }
    nOut <- nrow(xOut)# number of predictions
  }
  else{
    nOut <- 0
  }
  if(ncol(x) != ncol(xOut)){
    stop('x and xout must have the same number of columns')
  }
  if(nrow(x) != length(y)){
    stop('x and y must have the same number of observations')
  }
  n <- nrow(x)# number of observations
  p <- ncol(x)# number of covariates
  if(!is.na(sum(bw))){
    if(sum(bw<=0)>0){
      stop("bandwidth must be positive")
    }
    if(length(bw) != p){
      stop('Dimension of bandwidth does not agree with x')
    }
  }
  m <- ncol(y[[1]])# dimension of the correlation matrix
  yVec <- matrix(unlist(y), ncol = m^2, byrow = TRUE)# n by m^2
  if(substr(metric, 1, 1)=="p"){
    yAlpha <- lapply(y, function(yi) {
      eigenDecom <- eigen(yi)
      Lambda <- pmax(Re(eigenDecom$values), 0)# exclude 0i
      U <- eigenDecom$vectors
      U%*%diag(Lambda^alpha)%*%t(U)
    })
    yAlphaVec <- matrix(unlist(yAlpha), ncol = m^2, byrow = TRUE)# n by m^2
  }
  if(p > 3){
    warning('Local method is designed to work in low dimensional case, the result might be unstable.')
  }
  # define different kernels
  K_Gaussian <- function(x, h){
    x <- as.vector(x)
    p <- length(x)
    k <- 1
    for(i in 1:p){
      k <- k*1/sqrt(2*pi)*exp(-(x[i]/h[i])^2/2)/h[i]
    }
    return(as.numeric(k))
  }
  K_Uniform <- function(x, h){
    x <- as.vector(x)
    p <- length(x)
    k <- 1
    for(i in 1:p){
      k <- k*as.numeric(abs(x[i]/h[i])<=1)/2/h[i]
    }
    return(as.numeric(k))
  }
  K_Epanechnikov <- function(x, h){
    x <- as.vector(x)
    p <- length(x)
    k <- 1
    for(i in 1:p){
      k <- k*(abs(x[i]/h[i])<=1)*3/4*(1-(x[i]/h[i])^2)/h[i]
    }
    return(as.numeric(k))
  }
  if(substr(kernel, 1, 1) == 'g'){
    K <- K_Gaussian
  } else if(substr(kernel, 1, 1) == 'u'){
    K <- K_Uniform
  } else if(substr(kernel, 1, 1) == 'e'){
    K <- K_Epanechnikov
  } else{
    stop('Valid Kernels are Gaussian, Uniform and Epanechnikov')
  }
  
  # choose bandwidth by cross-validation
  if(is.na(sum(bw))){
    hs <- matrix(0, p, 20)
    for(l in 1:p){
      hs[l, ] <- exp(seq(from = log(n^(-1/(1+p))*(max(x[, l])-min(x[, l]))/10), 
                         to = log(5*n^(-1/(1+p))*(max(x[, l])-min(x[, l]))), 
                         length.out = 20))
    }
    cv <- array(0, 20^p)
    for(k in 0:(20^p-1)){
      h <- array(0, p)
      for(l in 1:p){
        kl <- floor((k %% (20^l)) / (20^(l-1))) + 1
        h[l] <- hs[l, kl]
      }
      for(j in 1:n){
        a <- x[j, ]
        if(p>1){
          mu1 <- rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi-a, h)*(xi-a)))
          mu2 <- matrix(rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi-a, h)*((xi-a) %*% t(xi-a)))), ncol=p)
        }
        else{
          mu1 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi-a, h)*(xi-a)))
          mu2 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi-a, h)*((xi-a) %*% t(xi-a))))
        }
        skip <- FALSE
        tryCatch(solve(mu2), error = function(e) skip <<- TRUE)
        if(skip){
          cv[k+1] <- Inf
          break
        }
        wc <- t(mu1)%*%solve(mu2)# 1 by p
        w <- apply(as.matrix(x[-j, ]), 1, function(xi){
          K(xi-a, h)*(1-wc%*%(xi-a))
        })# weight
        if(substr(metric, 1, 1)=="f"){
          qNew <- apply(yVec[-j, ], 2, weighted.mean, w)# m^2
          fitj <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
          fitj <- (fitj+t(fitj))/2# symmetrize
          if(!is.na(digits)) fitj <- round(fitj, digits = digits)# round
          cv[k+1] <- cv[k+1] + sum((y[[j]]-fitj)^2)/n
        }
        else if(substr(metric, 1, 1)=="p"){
          bAlpha <- matrix(apply(yAlphaVec[-j, ], 2, weighted.mean, w), ncol = m)# m by m
          eigenDecom <- eigen(bAlpha)
          Lambda <- pmax(Re(eigenDecom$values), 0)# projection to M_m
          U <- eigenDecom$vectors
          qNew <- as.vector(U%*%diag(Lambda^(1/alpha))%*%t(U))# inverse power
          fitj <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
          fitj <- (fitj+t(fitj))/2# symmetrize
          if(!is.na(digits)) fitj <- round(fitj, digits = digits)# round
          cv[k+1] <- cv[k+1] + sum((y[[j]]-fitj)^2)/n
        }
      }
    }
    bwi <- which.min(cv)
    bw <- array(0, p)
    for(l in 1:p){
      kl <- floor((bwi %% (20^l)) / (20^(l-1))) + 1
      bw[l] <- hs[l, kl]
    }
  }
  
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  if(substr(metric, 1, 1)=="f"){
    for(i in 1:n){
      a <- x[i, ]
      if(p>1){
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a)))), ncol=p)
      }
      else{
        mu1 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a))))
      }
      wc <- t(mu1)%*%solve(mu2)# 1 by p
      w <- apply(x, 1, function(xi){
        K(xi-a, bw)*(1-wc%*%(xi-a))
      })# weight
      qNew <- apply(yVec, 2, weighted.mean, w)# m^2
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp+t(temp))/2# symmetrize
      if(!is.na(digits)) temp <- round(temp, digits = digits)# round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]]-temp)^2))
    }
    if(nOut>0){
      predict <- vector(mode = "list", length = nOut)
      for(i in 1:nOut){
        a <- xOut[i, ]
        if(p>1){
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a)))), ncol=p)
        }
        else{
          mu1 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a))))
        }
        wc <- t(mu1)%*%solve(mu2)# 1 by p
        w <- apply(x, 1, function(xi){
          K(xi-a, bw)*(1-wc%*%(xi-a))
        })# weight
        qNew <- apply(yVec, 2, weighted.mean, w)# m^2
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp+t(temp))/2# symmetrize
        if(!is.na(digits)) temp <- round(temp, digits = digits)# round
        predict[[i]] <- temp
      }
      rst <- list(fit = fit, predict = predict, residuals = residuals, y = y, x = x, xOut = xOut, optns = optns)
    }
    else{
      rst <- list(fit = fit, residuals = residuals, y = y, x = x, optns = optns)
    }
  }
  else if(substr(metric, 1, 1)=="p"){
    for(i in 1:n){
      a <- x[i, ]
      if(p>1){
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a)))), ncol=p)
      }
      else{
        mu1 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a))))
      }
      wc <- t(mu1)%*%solve(mu2)# 1 by p
      w <- apply(x, 1, function(xi){
        K(xi-a, bw)*(1-wc%*%(xi-a))
      })# weight
      bAlpha <- matrix(apply(yAlphaVec, 2, weighted.mean, w), ncol = m)# m by m
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0)# projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U%*%diag(Lambda^(1/alpha))%*%t(U))# inverse power
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp+t(temp))/2# symmetrize
      if(!is.na(digits)) temp <- round(temp, digits = digits)# round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]]-temp)^2))
    }
    if(nOut>0){
      predict <- vector(mode = "list", length = nOut)
      for(i in 1:nOut){
        a <- xOut[i, ]
        if(p>1){
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a)))), ncol=p)
        }
        else{
          mu1 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*(xi-a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi-a, bw)*((xi-a) %*% t(xi-a))))
        }
        wc <- t(mu1)%*%solve(mu2)# 1 by p
        w <- apply(x, 1, function(xi){
          K(xi-a, bw)*(1-wc%*%(xi-a))
        })# weight
        bAlpha <- matrix(apply(yAlphaVec, 2, weighted.mean, w), ncol = m)# m^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0)# projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U%*%diag(Lambda^(1/alpha))%*%t(U))# inverse power
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp+t(temp))/2# symmetrize
        if(!is.na(digits)) temp <- round(temp, digits = digits)# round
        predict[[i]] <- temp
      }
      rst <- list(fit = fit, predict = predict, residuals = residuals, y = y, x = x, xOut = xOut, optns = optns)
    }
    else{
      rst <- list(fit = fit, residuals = residuals, y = y, x = x, optns = optns)
    }
  }
  class(rst) <- 'corReg'
  rst
}