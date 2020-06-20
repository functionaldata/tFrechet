#'@title Fréchet mean of covariance matrices
#'@description Fréchet mean computation for covariance matrices.
#'@param M A q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{metric}{Metric type choice, \code{"frobenius"}, \code{"power"}, \code{"log_cholesky"}, \code{"cholesky"} - default: \code{"frobenius"} which corresponds to the power metric with \code{alpha} equal to 2.}
#' \item{alpha}{The power parameter for the power metric, which can be any non-negative integer. Default is 2 which corresponds to Frobenius metric.}
#' }
#' @return A list containing the following fields:
#' \item{Mout}{A list containing the Fréchet mean of the covariance matrices in \code{M}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#'#Example M input
#'n=10 #sample size
#'m=5 # dimension of covariance matrices
#'M <- array(0,c(m,m,n))
#'for (i in 1:n){
#'  y0=rnorm(m)
#'  aux<-diag(m)+y0%*%t(y0)
#'  M[,,i]<-aux
#'}
#' Fmean=CovFMean(M=M,optns=list(metric="frobenius"))
#'
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \item \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}
#' \item \cite{Lin, Z. (2019). Riemannian geometry of symmetric positive definite matrices via Cholesky decomposition. Siam. J. Matrix. Anal, A. 40, 1353--1370.}
#' }
#' @export

CovFMean= function(M=NULL, optns = list()){
  if(is.list(M)){
    n=length(M)
  } else {
    if(!is.array(M)){
      stop('M must be an array or a list')
    }
    n=dim(M)[3]
  }
  if(n==1){
    stop("Sample size n should be at least 2")
  }
  x=matrix(1:n,nrow=n,ncol=1)
  xout=matrix((n+1)/2) #mean of x

  if (is.null(optns$metric)){
    metric="frobenius"
  } else {
    metric=optns$metric
  }
  if(!metric%in%c("frobenius","power","cholesky","log_cholesky")){
    stop("metric choice not supported.")
  }

  if(metric=="frobenius"){
    res <- list(Mout=GFRCov(x=x, y=NULL,M=M,xout=xout,optns = optns)$Mout,optns=optns)
  } else if(metric=="power"){
    res <- list(Mout=GFRCovPower(x=x, y=NULL,M=M,xout=xout,optns = optns)$Mout,optns=optns)
  } else {
    if (is.null(M))
      stop("M must be input for Cholesky and log-Cholesky metrics; y does not apply.")
    res <- list(Mout=GFRCovCholesky(x=x, M=M, xout=xout, optns = optns)$Mout,optns=optns)
  }
  class(res) <- "covReg"
  return(res)
}
