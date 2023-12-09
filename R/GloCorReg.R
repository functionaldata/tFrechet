#' @title Global Fréchet regression for correlation matrices
#' @description Global Fréchet regression for correlation matrices
#'   with Euclidean predictors.
#' @param x an n by p matrix or data frame of predictors.
#' @param M a q by q by n array (resp. a list of q by q matrices) where
#'   \code{M[, , i]} (resp. \code{M[[i]]}) contains the i-th correlation matrix
#'   of dimension q by q.
#' @param xOut an m by p matrix or data frame of output predictor levels.
#'   It can be a vector of length p if m = 1.
#' @param optns A list of options control parameters specified by
#'   \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{metric}{choice of metric. \code{'frobenius'} and \code{'power'} are supported,
#'   which corresponds to Frobenius metric and Euclidean power metric,
#'   respectively. Default is Frobenius metric.}
#' \item{alpha}{the power for Euclidean power metric.
#'   Default is 1 which corresponds to Frobenius metric.}
#' \item{digits}{the integer indicating the number of decimal places (round)
#'   to be kept in the output. Default is NULL, which means no round operation.}
#' }
#' @return A \code{corReg} object --- a list containing the following fields:
#' \item{fit}{a list of estimated correlation matrices at \code{x}.}
#' \item{predict}{a list of estimated correlation matrices at \code{xOut}.
#'   Included if \code{xOut} is not \code{NULL}.}
#' \item{RSquare}{Fréchet coefficient of determination.}
#' \item{AdjRSquare}{adjusted Fréchet coefficient of determination.}
#' \item{residuals}{Frobenius distance between the true and
#'   fitted correlation matrices.}
#' \item{xOut}{the output predictor level used.}
#' \item{optns}{the control options used.}
#' @examples
#' # Generate simulation data
#' n <- 100
#' q <- 10
#' d <- q * (q - 1) / 2
#' xOut <- seq(0.1, 0.9, length.out = 9)
#' x <- runif(n, min = 0, max = 1)
#' y <- list()
#' for (i in 1:n) {
#'   yVec <- rbeta(d, shape1 = x[i], shape2 = 1 - x[i])
#'   y[[i]] <- matrix(0, nrow = q, ncol = q)
#'   y[[i]][lower.tri(y[[i]])] <- yVec
#'   y[[i]] <- y[[i]] + t(y[[i]])
#'   diag(y[[i]]) <- 1
#' }
#' # Frobenius metric
#' fit1 <- GloCorReg(x, y, xOut,
#'   optns = list(metric = "frobenius", digits = 5)
#' )
#' # Euclidean power metric
#' fit2 <- GloCorReg(x, y, xOut,
#'   optns = list(metric = "power", alpha = .5)
#' )
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @importFrom Matrix nearPD
#' @export

GloCorReg <- function(x, M, xOut = NULL, optns = list()) {
  if (is.null(optns$metric)) {
    metric <- "frobenius"
  } else {
    metric <- optns$metric
  }
  if (!(metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    alpha <- 1
  } else {
    alpha <- optns$alpha
  }
  if (alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$digits)) {
    digits <- NA
  } else {
    digits <- optns$digits
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame")
    }
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (ncol(x) == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame")
      }
    }
    if (ncol(x) != ncol(xOut)) {
      stop("x and xOut must have the same number of columns")
    }
    m <- nrow(xOut) # number of predictions
  } else {
    m <- 0
  }
  y <- M
  if (!is.list(y)) {
    if (is.array(y)) {
      y <- lapply(seq(dim(y)[3]), function(i) y[, , i])
    } else {
      stop("y must be a list or an array")
    }
  }
  if (nrow(x) != length(y)) {
    stop("the number of rows in x must be the same as the number of correlation matrices in y")
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  q <- ncol(y[[1]]) # dimension of the correlation matrix
  yVec <- matrix(unlist(y), ncol = q^2, byrow = TRUE) # n by q^2
  if (substr(metric, 1, 1) == "p") {
    yAlpha <- lapply(y, function(yi) {
      eigenDecom <- eigen(yi)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^alpha) %*% t(U)
    })
    yAlphaVec <- matrix(unlist(yAlpha), ncol = q^2, byrow = TRUE) # n by q^2
  }
  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  wc <- t(apply(x, 1, function(xi) t(xi - xMean) %*% invVa)) # n by p
  totVa <- sum((scale(yVec, scale = FALSE))^2)
  if (nrow(wc) != n) wc <- t(wc) # for p=1
  if (substr(metric, 1, 1) == "f") {
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      qNew <- apply(yVec, 2, weighted.mean, w) # q^2
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(digits)) temp <- round(temp, digits = digits) # round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]] - temp)^2))
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (m > 0) {
      predict <- vector(mode = "list", length = m)
      for (i in 1:m) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        qNew <- apply(yVec, 2, weighted.mean, w) # q^2
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(digits)) temp <- round(temp, digits = digits) # round
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, optns = optns)
    }
  } else if (substr(metric, 1, 1) == "p") {
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      bAlpha <- matrix(apply(yAlphaVec, 2, weighted.mean, w), ncol = q) # q by q
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / alpha)) %*% t(U)) # inverse power
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(digits)) temp <- round(temp, digits = digits) # round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]] - temp)^2))
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (m > 0) {
      predict <- vector(mode = "list", length = m)
      for (i in 1:m) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        bAlpha <- matrix(apply(yAlphaVec, 2, weighted.mean, w), ncol = q) # q^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / alpha)) %*% t(U)) # inverse power
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(digits)) temp <- round(temp, digits = digits) # round
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, optns = optns)
    }
  }
  class(res) <- "corReg"
  res
}
