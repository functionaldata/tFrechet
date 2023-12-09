#' @title Local Fréchet regression for correlation matrices
#' @description Local Fréchet regression for correlation matrices
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
#' \item{alpha}{the power for Euclidean power metric. Default is 1 which
#'   corresponds to Frobenius metric.}
#' \item{kernel}{Name of the kernel function to be chosen from \code{'gauss'},
#'   \code{'rect'}, \code{'epan'}, \code{'gausvar'} and \code{'quar'}. Default is \code{'gauss'}.}
#' \item{bw}{bandwidth for local Fréchet regression, if not entered
#'   it would be chosen from cross validation.}
#' \item{digits}{the integer indicating the number of decimal places (round)
#'   to be kept in the output. Default is NULL, which means no round operation.}
#' }
#' @return A \code{corReg} object --- a list containing the following fields:
#' \item{fit}{a list of estimated correlation matrices at \code{x}.}
#' \item{predict}{a list of estimated correlation matrices at \code{xOut}.
#'   Included if \code{xOut} is not \code{NULL}.}
#' \item{residuals}{Frobenius distance between the true and
#'   fitted correlation matrices.}
#' \item{xOut}{the output predictor level used.}
#' \item{optns}{the control options used.}
#' @examples
#' # Generate simulation data
#' \donttest{
#' n <- 100
#' q <- 10
#' d <- q * (q - 1) / 2
#' xOut <- seq(0.1, 0.9, length.out = 9)
#' x <- runif(n, min = 0, max = 1)
#' y <- list()
#' for (i in 1:n) {
#'   yVec <- rbeta(d, shape1 = sin(pi * x[i]), shape2 = 1 - sin(pi * x[i]))
#'   y[[i]] <- matrix(0, nrow = q, ncol = q)
#'   y[[i]][lower.tri(y[[i]])] <- yVec
#'   y[[i]] <- y[[i]] + t(y[[i]])
#'   diag(y[[i]]) <- 1
#' }
#' # Frobenius metric
#' fit1 <- LocCorReg(x, y, xOut,
#'   optns = list(metric = "frobenius", digits = 2)
#' )
#' # Euclidean power metric
#' fit2 <- LocCorReg(x, y, xOut,
#'   optns = list(
#'     metric = "power", alpha = .5,
#'     kernel = "epan", bw = 0.08
#'   )
#' )
#' }
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @importFrom Matrix nearPD
#' @export

LocCorReg <- function(x, M, xOut = NULL, optns = list()) {
  if (is.null(optns$metric)) {
    metric <- "frobenius"
  } else {
    metric <- optns$metric
  }
  if (!(metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$kernel)) {
    kernel <- "gauss"
  } else {
    kernel <- optns$kernel
  }
  if (is.null(optns$bw)) {
    bw <- NA
  } else {
    bw <- optns$bw
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
  if (p > 2) {
    stop("local method is designed to work in low dimensional case (p is either 1 or 2)")
  }
  if (!is.na(sum(bw))) {
    if (sum(bw <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(bw) != p) {
      stop("dimension of bandwidth does not agree with x")
    }
  }
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
  # define different kernels
  Kern <- kerFctn(kernel)
  K <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * Kern(x[i] / h[i])
    }
    return(as.numeric(k))
  }

  # choose bandwidth by cross-validation
  if (is.na(sum(bw))) {
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
        if (substr(metric, 1, 1) == "f") {
          qNew <- apply(yVec[-j, ], 2, weighted.mean, w) # q^2
          fitj <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(digits)) fitj <- round(fitj, digits = digits) # round
          cv[k + 1] <- cv[k + 1] + sum((y[[j]] - fitj)^2) / n
        } else if (substr(metric, 1, 1) == "p") {
          bAlpha <- matrix(apply(yAlphaVec[-j, ], 2, weighted.mean, w), ncol = q) # q by q
          eigenDecom <- eigen(bAlpha)
          Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
          U <- eigenDecom$vectors
          qNew <- as.vector(U %*% diag(Lambda^(1 / alpha)) %*% t(U)) # inverse power
          fitj <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(digits)) fitj <- round(fitj, digits = digits) # round
          cv[k + 1] <- cv[k + 1] + sum((y[[j]] - fitj)^2) / n
        }
      }
    }
    bwi <- which.min(cv)
    bw <- array(0, p)
    for (l in 1:p) {
      kl <- floor((bwi %% (20^l)) / (20^(l - 1))) + 1
      bw[l] <- hs[l, kl]
    }
  }

  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  if (substr(metric, 1, 1) == "f") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, bw) * (1 - wc %*% (xi - a))
      }) # weight
      qNew <- apply(yVec, 2, weighted.mean, w) # q^2
      temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(digits)) temp <- round(temp, digits = digits) # round
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((y[[i]] - temp)^2))
    }
    if (m > 0) {
      predict <- vector(mode = "list", length = m)
      for (i in 1:m) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, bw) * (1 - wc %*% (xi - a))
        }) # weight
        qNew <- apply(yVec, 2, weighted.mean, w) # q^2
        temp <- as.matrix(Matrix::nearPD(matrix(qNew, ncol = q), corr = TRUE, maxit = 1000)$mat)
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(digits)) temp <- round(temp, digits = digits) # round
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, optns = optns)
    }
  } else if (substr(metric, 1, 1) == "p") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, bw) * (1 - wc %*% (xi - a))
      }) # weight
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
    if (m > 0) {
      predict <- vector(mode = "list", length = m)
      for (i in 1:m) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, bw) * (1 - wc %*% (xi - a))
        }) # weight
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
      res <- list(fit = fit, predict = predict, residuals = residuals, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, optns = optns)
    }
  }
  class(res) <- "corReg"
  res
}
