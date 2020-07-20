#'@title Power distance for covariance matrices
#'@description Power distance computation for covariance matrices.
#'@param M A q by q by 2 array (resp. a list of two q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{metric}{Metric type choice, \code{"frobenius"}, \code{"power"} - default: \code{"frobenius"} which corresponds to the power metric with \code{alpha} equal to 1.}
#' \item{alpha}{The power parameter for the power metric, which can be any non-negative number. Default is 1 which corresponds to frobenius metric.}
#' }
#' @return A list containing the following fields:
#' \item{dist}{The distance between the two covariance matrices in \code{M}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#'#Example M input as array
#'m=5 # dimension of covariance matrices
#'M <- array(0,c(m,m,2))
#'for (i in 1:2){
#'  y0=rnorm(m)
#'  aux<-diag(m)+y0%*%t(y0)
#'  M[,,i]<-aux
#'}
#' covDistance=CovFPowerDist(M=M,optns=list(metric="frobenius"))
#'
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2016). Fréchet integration and adaptive metric selection for interpretable covariances of multivariate functional data. Biometrika, 103, 103--120.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \item \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}
#' }
#' @export

CovFPowerDist= function(M=NULL, optns = list()){
  
  if(is.null(M)){
    stop("M must be provided")
  }
  if(is.list(M)){
    n=length(M)
    M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
  } else{
    if(!is.array(M)){
      stop('M must be an array or a list')
    }
    n=dim(M)[3]
  }
  if(n>2){
    stop('Only two covariance matrices should be provided')
  }
  
  if (is.null(optns$metric)){
    metric="frobenius"
  } else {
    metric=optns$metric
    if(is.null(optns$alpha)){
      optns$alpha=1
      alpha=1
    }
    else{
      alpha=optns$alpha
    }
  }
  if(!metric%in%c("frobenius","power")){
    stop("metric choice not supported.")
  }
  
  if(metric=="frobenius"){
    dist <- sqrt(sum(diag((M[,,1]-M[,,2])%*%(M[,,1]-M[,,2]))))
  } else if(metric=="power"){
    if(alpha==1){
      dist <- sqrt(sum(diag((M[,,1]-M[,,2])%*%(M[,,1]-M[,,2]))))
    }
    else{
      if(alpha>0){
        #Transform M1
        P=eigen(M[,,1])$vectors
        Lambd_alpha=diag(pmax(0,eigen(M[,,1])$values)**alpha)
        M_alpha_1=P%*%Lambd_alpha%*%t(P)
        #Transform M2
        P=eigen(M[,,2])$vectors
        Lambd_alpha=diag(pmax(0,eigen(M[,,2])$values)**alpha)
        M_alpha_2=P%*%Lambd_alpha%*%t(P)       
        dist <- sqrt(sum(diag((M_alpha_1-M_alpha_2)%*%(M_alpha_1-M_alpha_2))))/alpha
      }
      else{
        if(alpha<0){
          stop('alpha has to be non-negative')
        }
        #Transform M1
        P=eigen(M[,,1])$vectors
        Lambd_alpha=diag(log(pmax(1e-30,eigen(M[,,1])$values)))
        M_alpha_1=P%*%Lambd_alpha%*%t(P)
        #Transform M2
        P=eigen(M[,,2])$vectors
        Lambd_alpha=diag(log(pmax(1e-30,eigen(M[,,2])$values)))
        M_alpha_2=P%*%Lambd_alpha%*%t(P)
        dist <- sqrt(sum(diag((M_alpha_1-M_alpha_2)%*%(M_alpha_1-M_alpha_2))))
      }
    }
  }
  res <- list(dist=dist,optns=optns)
  return(res)
}
