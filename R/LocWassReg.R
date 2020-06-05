#' @title Local Wasserstein Regression
#' @noRd
#' @description  Local Fréchet regression with respect to the Wasserstein distance.
#'
#' @param xin A vector of length n with input measurement points.
#' @param qin An n by m matrix with values of quantile functions of which each row holds the quantile function values on an equispaced grid on [0, 1] of length m.
#' @param xout A vector of length k with output measurement points.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bw}{A scalar used as the bandwidth.}
#' \item{ker}{A character holding the type of kernel functions.}
#' \item{lower}{A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.}
#' \item{upper}{A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.}
#' }

LocWassReg = function(xin, qin, xout, optns = list()){

  if (is.null(optns$bw))
    stop ("optns$bw has no default values and must be input by users.")
  bw <- optns$bw
  if(!is.numeric(optns$bw) | (length(optns$bw)>1))
    stop("optns$bw should be a numerical vector of length 1.")
  if(!is.vector(xin) | !is.numeric(xin))
    stop("xin should be a numerical vector.")
  if(!is.matrix(qin) | !is.numeric(qin))
    stop("qin should be a numerical matrix.")
  if(!is.vector(xout) | !is.numeric(xout))
    stop("xout should be a numerical vector.")
  if(length(xin)!=nrow(qin))
    stop("The length of xin should be the same as the number of rows of qin.")

  #if (is.null(optns$ker)) optns$ker <- 'gauss'
  ker <- kerFctn(optns$ker)

  k <- length(xout)
  n <- length(xin)
  m <- ncol(qin)

  # if lower & upper are neither NULL
  A <- cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  if (!is.null(optns$upper) & !is.null(optns$lower)) {
    b0 <- c(optns$lower, rep(0,m-1), -optns$upper)
  } else if(!is.null(optns$upper)) {
    A <- A[,-1]
    b0 <- c(rep(0,m-1), -optns$upper)
  } else if(!is.null(optns$lower)) {
    A <- A[,-ncol(A)]
    b0 <- c(optns$lower,rep(0,m-1))
  } else {
    A <- A[,-c(1,ncol(A))]
    b0 <- rep(0,m-1)
  }
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")

  qout <- sapply(1:k, function(j){
    mu0 <- mean(ker((xout[j] - xin) / bw))
    mu1 <- mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
    mu2 <- mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
    s <- ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
      (mu0 * mu2 - mu1^2)
    gx <- colMeans(qin * s)

    #res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
    #return(sort(res$solution))

    res <- do.call(rosqp::solve_osqp,
                   list(P=Pmat, q= -gx, A=Amat, l=b0, pars = rosqp::osqpSettings(verbose = FALSE)))
    return(sort(res$x))
  })
  qout = t(qout)
  return(qout)
}
