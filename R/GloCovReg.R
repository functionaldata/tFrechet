#'@title Global Fréchet regression of covariance matrices
#'@description Global Fréchet regression of covariance matrices with Euclidean predictors.
#'@param x An n by p matrix of predictors.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed. See 'metric' option in 'Details' for more details.
#'@param M A q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.  See 'metric' option in 'Details' for more details.
#'@param xout An m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{metric}{Metric type choice, \code{"frobenius"}, \code{"power"}, \code{"log_cholesky"}, \code{"cholesky"} - default: \code{"frobenius"} which corresponds to the power metric with \code{alpha} equal to 2.
#' For power (and Frobenius) metrics, either \code{y} or \code{M} must be input; \code{y} would override \code{M}. For Cholesky and log-Cholesky metrics, \code{M} must be input and \code{y} does not apply.}
#' \item{alpha}{The power parameter for the power metric, which can be any non-negative integer. Default is 2 which corresponds to Frobenius metric.}
#' }
#' @return A \code{covReg} object --- a list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance or correlation matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#'#Example y input
#'n=50             # sample size
#'t=seq(0,1,length.out=100)       # length of data
#'x = matrix(runif(n),n)
#'theta1 = theta2 = array(0,n)
#'for(i in 1:n){
#'  theta1[i] = rnorm(1,x[i],x[i]^2)
#'  theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
#'}
#'y = matrix(0,n,length(t))
#'phi1 = sqrt(3)*t
#'phi2 = sqrt(6/5)*(1-t/2)
#'y = theta1%*%t(phi1) + theta2 %*% t(phi2)
#'xout = matrix(c(0.25,0.5,0.75),3)
#'Cov_est=GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3))
#'#Example M input
#'n=10 #sample size
#'m=5 # dimension of covariance matrices
#'M <- array(0,c(m,m,n))
#'for (i in 1:n){
#'  y0=rnorm(m)
#'  aux<-diag(m)+y0%*%t(y0)
#'  M[,,i]<-aux
#'}
#'x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
#'xout=cbind(runif(3),runif(3)) #output predictor levels
#'Cov_est=GloCovReg(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,metric="power",alpha=3))
#'
#' @references
#' \itemize{
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \item \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}
#' \item \cite{Lin, Z. (2019). Riemannian geometry of symmetric positive definite matrices via Cholesky decomposition. Siam. J. Matrix. Anal, A. 40, 1353--1370.}
#' }
#' @export

GloCovReg= function(x,y=NULL,M=NULL,xout,optns = list()){
  if (is.null(optns$metric)){
    metric="frobenius"
  } else {
    metric=optns$metric
  }
  if(!metric%in%c("frobenius","power","cholesky","log_cholesky")){
    stop("metric choice not supported.")
  }
  if(metric=="frobenius"){
    res <- GFRCov(x=x, y=y,M=M,xout=xout,optns = optns)
  } else if(metric=="power"){
    res <- GFRCovPower(x=x, y=y,M=M,xout=xout,optns = optns)
  } else {
    if (is.null(M))
      stop("M must be input for Cholesky and log-Cholesky metrics; y does not apply.")
    res <- GFRCovCholesky(x=x, M=M, xout=xout, optns = optns)
  }
  class(res) <- "covReg"
  return(res)
}
