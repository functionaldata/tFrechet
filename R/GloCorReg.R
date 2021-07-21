#'@title Global Fr$\'{e}$chet regression for correlation matrices
#'@description Global Fr$\'{e}$chet regression for correlation matrices with Euclidean predictors.
#'@param x an n by p matrix of predictors.
#'@param y a list (length n) of $m$ by $m$ correlation matrix.
#'@param xOut an nOut by p matrix of output predictor levels.
#'@param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'@details Available control options are
#'\describe{
#'\item{metric}{choice of metric. 'frobenius' and 'power' are supported, which corresponds to Frobenius metric and Euclidean power metric, repectively. Default is Frobenius metric}
#'\item{alpha}{the power for Euclidean power metric. Default is $1$ which corresponds to Frobenius metric.}
#'\item{digits}{the integer indicating the number of decimal places (round) to be kept in the output. Default is NULL, which means no round operation.}
#'}
#'@return A \code{corReg} object --- a list containing the follwing fields:
#'\item{fit}{a list of estimated correlation matrices at \code{x}.}
#'\item{predict}{a list of estimated correlation matrices at \code{xOut}. Included if \code{xOut} is not \code{NULL}.}
#'\item{RSquare}{Fr\'{e}chet coefficient of determination.}
#'\item{AdjRSquare}{adjusted Fr\'{e}chet coefficient of determination.}
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
#'  yVec <- rbeta(d, shape1 = x[i], shape2 = 1-x[i])
#'  y[[i]] <- matrix(0, nrow = m, ncol = m)
#'  y[[i]][lower.tri(y[[i]])] <- yVec
#'  y[[i]] <- y[[i]] + t(y[[i]])
#'  diag(y[[i]]) <- 1
#'}
#'# Frobenius metric
#'fit1 <- GloCorReg(x, y, xOut, 
#'                  optns = list(metric = 'frobenius', digits = 5))
#'# Euclidean power metric
#'fit2 <- GloCorReg(x, y, xOut, 
#'                  optns = list(metric = 'power', alpha = .5))
#'@references
#'\itemize{
#'\item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#'}
#'@export

GloCorReg <- function(x, y, xOut = NULL, optns = list()){
  if(is.null(optns$metric)){
    metric <- "frobenius"
  } else {
    metric <- optns$metric
  }
  if(!(metric %in% c("frobenius", "power"))){
    stop("metric choice not supported.")
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
  xMean <- colMeans(x)
  invVa <- solve(var(x)*(n-1)/n)
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  wc <- t(apply(x, 1, function(xi) t(xi-xMean)%*%invVa))# n by p
  totVa <- sum((scale(yVec, scale = FALSE))^2)
  if(nrow(wc)!=n) wc <- t(wc)# for p=1
  if(substr(metric, 1, 1)=="f"){
    for(i in 1:n){
      w <- apply(wc, 1, function(wci) 1+t(wci)%*%(x[i, ]-xMean))
      qNew <- apply(yVec, 2, weighted.mean, w)# m^2
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp+t(temp))/2# symmetrize
      if(!is.na(digits)) temp <- round(temp, digits = digits)# round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]]-temp)^2))
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa/totVa
    AdjRSquare <- RSquare-(1-RSquare)*p/(n-p-1)
    if(nOut>0){
      predict <- vector(mode = "list", length = nOut)
      for(i in 1:nOut){
        w <- apply(wc, 1, function(wci) 1+t(wci)%*%(xOut[i, ]-xMean))
        qNew <- apply(yVec, 2, weighted.mean, w)# m^2
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = m), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp+t(temp))/2# symmetrize
        if(!is.na(digits)) temp <- round(temp, digits = digits)# round
        predict[[i]] <- temp
      }
      rst <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, y = y, x = x, xOut = xOut, optns = optns)
    }
    else{
      rst <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, y = y, x = x, optns = optns)
    }
  }
  else if(substr(metric, 1, 1)=="p"){
    for(i in 1:n){
      w <- apply(wc, 1, function(wci) 1+t(wci)%*%(x[i, ]-xMean))
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
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa/totVa
    AdjRSquare <- RSquare-(1-RSquare)*p/(n-p-1)
    if(nOut>0){
      predict <- vector(mode = "list", length = nOut)
      for(i in 1:nOut){
        w <- apply(wc, 1, function(wci) 1+t(wci)%*%(xOut[i, ]-xMean))
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
      rst <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, y = y, x = x, xOut = xOut, optns = optns)
    }
    else{
      rst <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, y = y, x = x, optns = optns)
    }
  }
  class(rst) <- 'corReg'
  rst
}