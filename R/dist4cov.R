#' @title Distance between covariance matrices 
#'  
#' @description Distance computation between two covariance matrices
#' @param A an p by p matrix
#' @param B an p by p matrix
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{metric}{Metric type choice, \code{"frobenius"}, \code{"power"}, \code{"log_cholesky"} and \code{"cholesky"} - default: \code{"frobenius"}, which corresponds to the power metric with \code{alpha} equal to 1.}
#' \item{alpha}{The power parameter for the power metric, which can be any non-negative number. Default is 1 which corresponds to Frobenius metric.}
#' }
#' @return A list containing the following fields:
#' \item{dist}{the distance between covariance matrices \code{A} and \code{B}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples 
#'# M input as array
#' m <- 5 # dimension of covariance matrices
#' M <- array(0,c(m,m,2))
#' for (i in 1:2) {
#'  y0 <- rnorm(m)
#'  aux <- diag(m) + y0 %*% t(y0)
#'  M[,,i] <- aux
#' }
#' A <- M[,,1]
#' B <- M[,,2]
#' frobDist <- dist4cov(A=A, B=B, optns=list(metric="frobenius"))
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2016). Fréchet integration and adaptive metric selection for interpretable covariances of multivariate functional data. Biometrika, 103, 103--120.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \item \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}
#' }
#' @export

dist4cov= function(A=NULL,B=NULL, optns = list()){
  if(!is.matrix(A) | !is.matrix(B) ){
    stop('A and B must be of matrix class')
  }
  if(nrow(A)!=ncol(A) | nrow(B)!=ncol(B)){
    stop('Both A and B must be square matrices')
  }
  if(sum(dim(A)!=dim(B))>0){
    stop('Both A and B must have the same dimension')
  }
  
  if (is.null(optns$metric)){
    metric="frobenius"
  } else {
    metric=optns$metric
    if(metric%in%c("power")){
      if(is.null(optns$alpha)){
        optns$alpha=1
      }
    }
  }
  if(!metric%in%c("frobenius","power","log_cholesky","cholesky")){
    stop("metric choice is not supported.")
  }
  if(metric%in%c("frobenius","power")){
    res=CovFPowerDist(A=A,B=B,optns=optns)
  }
  else{
    res=DistCholesky(A=A,B=B,optns=optns)
  }
  return(res)
}
  
  
  
  
  
  