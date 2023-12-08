#' @title Glocal Wasserstein Regression
#' @noRd
#' @description  Glocal Frechet regression with respect to the Wasserstein distance.
#'
#' @param xin An n by p matrix or a vector of length n (if p=1) with input measurements of the predictors.
#' @param qin An n by m matrix with values of quantile functions of which each row holds the quantile function values on an equispaced grid on [0, 1].
#' @param xout A k by p matrix or a vector of length k (if p=1) with output measurements of the predictors.
#' @param optns A list of control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{lower}{A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.}
#' \item{upper}{A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.}
#' \item{Rsquared}{A logical variable indicating whether R squared would be returned. Default is \code{FALSE}.}
#' \item{qSup}{A numerical vector of length m holding the probability grid on [0, 1] at which the input quantile functions take values. If \code{optns$Rsquared} is TRUE, \code{qSup} is needed. Default is \code{seq(1,2*m,2)/2/m}.}
#' }
#' @importFrom osqp solve_osqp osqpSettings
#' @importFrom pracma trapz

GloWassReg <- function(xin, qin, xout, optns=list()){
  if (is.null(optns$Rsquared)) optns$Rsquared <- FALSE

  if(is.vector(xin)){
    xin <- as.matrix(xin)
  }
  if(is.vector(xout)){
    xout <- as.matrix(xout)
  }
  if(nrow(xin)!=nrow(qin))
    stop("The numbers of observations in xin and qin are not the same.")
  if(ncol(xin)!=ncol(xout))
    stop("The numbers of variables in xin and xout are not the same.")
  if(optns$Rsquared & is.null(optns$qSup)){
    warning("optns$qSup is missing and taking the default value.")
  }

  k <- nrow(xout)
  n <- nrow(xin)
  m <- ncol(qin)
  xbar <- colMeans(xin)
  Sigma <- cov(xin) * (n-1) / n
  invSigma <- solve(Sigma)

  # if lower & upper are neither NULL
  A <- cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  if(!is.null(optns$upper) & !is.null(optns$lower)){
    b0 <- c(optns$lower, rep(0,m-1), -optns$upper)
  }else if(!is.null(optns$upper)){
    A <- A[,-1]
    b0 <- c(rep(0,m-1), -optns$upper)
  }else if(!is.null(optns$lower)){
    A <- A[,-ncol(A)]
    b0 <- c(optns$lower,rep(0,m-1))
  }else{
    A <- A[,-c(1,ncol(A))]
    b0 <- rep(0,m-1)
  }
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")

  qout <- sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    gx <- colMeans(qin * s)

    #res <- do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
    #return(sort(res$solution)) #return(res$solution)


    res <- do.call(osqp::solve_osqp,
                   list(P=Pmat, q= -gx, A=Amat, l=b0, pars = osqp::osqpSettings(verbose = FALSE)))
    return(sort(res$x))
  })
  qout <- t(qout)

  if (!optns$Rsquared) {
    return(list(qout=qout))
  } else {
    qMean <- colMeans(qin)
    if (k == n) {
      if (sum(abs(xout-xin)) > 1e-10*length(xout))
        qin.est <- qout
    } else {
      qin.est <- sapply(1:n, function(j){
        s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xin[j,] - xbar)
        s <- as.vector(s)
        gx <- colMeans(qin * s)

        #res <- do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
        #return(sort(res$solution)) #return(res$solution)

        res <- do.call(osqp::solve_osqp,
                       list(P=Pmat, q= -gx, A=Amat, l=b0, pars = osqp::osqpSettings(verbose = FALSE)))
        return(sort(res$x))
      })
      qin.est <- t(qin.est)
    }
    Rsq <- ifelse(
      is.null(optns$qSup),
      1 - sum(t(qin - qin.est)^2) / sum((t(qin) - qMean)^2),
      1 - pracma::trapz(x=optns$qSup, y=colSums((qin - qin.est)^2)) /
        pracma::trapz(x=optns$qSup, y=rowSums((t(qin) - qMean)^2))
    )
    if(Rsq < 0) Rsq <- 0
    return(list(qout=qout, R.squared=Rsq))
  }
}
