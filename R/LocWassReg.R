#' @title Local Wasserstein Regression
#' @noRd
#' @description  Local Fréchet regression for probability distributions with respect to the Wasserstein distance.
#'
#' @param xin A n by p matrix holding the n observations of the predictor.
#' @param qin An n by m matrix with values of quantile functions of which each row holds the quantile function values on an equispaced grid on [0, 1] of length m.
#' @param xout A k by p matrix holding the k output predictor values.
#' @param optns A list of control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bw}{A vector of length p containing the bandwidths of each predictor dimension.}
#' \item{ker}{A character holding the type of kernel functions.}
#' \item{lower}{A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.}
#' \item{upper}{A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.}
#' }

LocWassReg = function(xin, qin, xout, optns = list()){

  if(!is.matrix(xin)&!is.vector(xin)){
    stop('xin must be a matrix or vector')
  }
  if(is.vector(xin)){
    xin<- matrix(xin,length(xin))
  }
  if (is.null(xout)){
    xout <- xin
  }
  if(!is.matrix(xout)&!is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  if(is.vector(xout)){
    xout<- matrix(xout,length(xout))
  }
  if(ncol(xin) != ncol(xout)){
    stop('xin and xout must have the same number of columns')
  }
  if(is.null(optns$bw)){
    stop ("optns$bw has no default values and must be input by user.")
  }
  if(!is.numeric(optns$bw) | (length(optns$bw)!=ncol(xin))){
    stop("optns$bw should be a numerical vector of length p.")
  }
  if(!is.matrix(qin) | !is.numeric(qin)){
    stop("qin should be a numerical matrix.")
  }
  if(nrow(xin)!=nrow(qin)){
    stop("The number of rows of xin should be the same as the number of rows of qin.")
  }
  if(is.null(optns$ker)){
    optns$ker <- 'gauss'
  }
  ker <- kerFctn(optns$ker)
  
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*ker(x[,i]/h[i])
    }
    return(as.numeric(k))
  }
  
  k <- nrow(xout)
  n <- nrow(xin)
  m <- ncol(qin)
  p <- ncol(xin)
  
  getLFRweights=function(x0){
    #x0 is a vector in R^p that corresponds to the covariate value for which we want to predict
    aux=K(xin-matrix(t(x0),nrow=n,ncol=length(x0),byrow=TRUE),optns$bw)
    mu0 = mean(aux)
    mu1 = colMeans(aux*(xin - matrix(t(x0),nrow=n,ncol=length(x0),byrow=TRUE)))
    mu2=0
    for(i in 1:n){
      mu2 = mu2 + aux[i]*(xin[i,]-x0) %*% t(xin[i,]-x0)/n
    }
    sL = array(0,n)
    for(i in 1:n){
      sL[i] =aux[i]*(1-t(mu1)%*%solve(mu2)%*%(xin[i,]-x0))
    }
    s = sum(sL)
    return(sL/s)
  }
  
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
    s=getLFRweights(xout[j,])
    s=as.vector(s)
    gx <- colMeans(qin * s)*n
    #res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
    #return(sort(res$solution))
    res <- do.call(osqp::solve_osqp,
                   list(P=Pmat, q= -gx, A=Amat, l=b0, pars = osqp::osqpSettings(verbose = FALSE)))
    return(sort(res$x))
  })
  qout = t(qout)
  return(qout)
}
