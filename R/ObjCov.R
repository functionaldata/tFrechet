#' @title Object Covariance
#' @description Calculating covariance for time varying object data
#' @param tgrid Time grid for the time varying object data and covariance function
#' @param I A four dimension array of \code{n} x \code{n} matrix of squared distances between the time point u of the ith process and process and the time point v of the jth object process, 
#' e.g.: \eqn{I[i,j,u,v] = d^2(X_i(u) X_j(v))}
#' @param K Numbers of principal components 
#' @param smooth Logical indicating if the smoothing is enabled when calculating the eigenvalues and eigenfunctions
#' @return A list of the following:
#' \item{C}{Estimated object covariance (non-smooth) on the 2D grid of dimension \code{length(tgrid)} X \code{length(tgrid)}}
#' \item{sC}{Estimated object covariance (smooth) on the 2D grid of dimension \code{length(tgrid)} X \code{length(tgrid)}}
#' \item{tgrid}{Time grid for the time varying object data and covariance function}
#' \item{K}{Numbers of principal components}
#' \item{phi}{Matrix of smooth eigenfunctions (dimension: \code{length(tgrid)} X \code{K})}
#' \item{lambda}{Vector of eigenvalues of dimension \code{K} }
#' @examples 
#' \donttest{
#' ### functional covariate
#' phi1 <- function(x) -cos(pi*x/10)/sqrt(5)
#' phi2 <- function(x)  sin(pi*x/10)/sqrt(5)
#' 
#' lambdaX <- c(4,2)
#' # training set
#' n <- 100
#' N <- 50
#' tgrid <- seq(0,10,length.out = N)
#' 
#' Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
#' CovX <- lambdaX[1] * phi1(tgrid) %*% t(phi1(tgrid)) + lambdaX[2] * phi2(tgrid) %*% t(phi2(tgrid))
#' comp1 = lambdaX[1]^(1/2) * Xi[,1] %*% t(phi1(tgrid))
#' comp2 = lambdaX[2]^(1/2) * Xi[,2] %*% t(phi2(tgrid))
#' SampleX <- comp1 + comp2
#' 
#' I <- array(0, c(n,n,N,N))
#' for (u in 1:N){
#'   for (v in 1:N){
#'     temp1 <- SampleX[,u]
#'    temp2 <- SampleX[,v]
#'     I[,,u,v] <- outer(temp1, temp2, function(v1,v2){
#'      (v1 - v2)^2
#'    })
#'  }
#' }
#' 
#' result_cov <- ObjCov(tgrid, I, 2)
#' result_cov$lambda #4 2
#' 
#' sC <- result_cov$sC
#' sum((sC-CovX)^2) / sum(sC^2)
#' sum((phi1(tgrid)-result_cov$phi[,1])^2)/sum(phi1(tgrid)^2)
#' }
#' @references 
#' \cite{Dubey, P., & Müller, H. G. (2020). Functional models for time‐varying random objects. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(2), 275-327.}
#' @export
#' @import pracma
#' @import fdapace
#' @importFrom utils getFromNamespace


ObjCov <- function(tgrid, I, K, smooth=TRUE){
  if(length(dim(I))!=4){
    stop("I must be a four dimensional array")
  }
  n = dim(I)[1]
  N = length(tgrid)
  if(dim(I)[2]!=n){
    stop("length of the first and second arguments of array I are inconsistent")
  }
  if(dim(I)[3]!=dim(I)[4]){
    stop("length of the third and fourth arguments of array I are inconsistent")
  }
  if(dim(I)[3]!=N){
    stop("length of tgrid and the length of third arguments of array I are inconsistent")
  }
  t <- length(tgrid)
  C <- matrix(0, nrow = t, ncol = t)
  # calculate empirical metric autocovariance matrix C 
  for(u in 1:(t-1)){
    for(v in (u+1):t){
      C[u,v] <-  ((sum(I[,,u,v]) - sum(diag(I[,,u,v]))) - (n-1)*sum(diag(I[,,u,v]))) / (2*n*(n-1))
      C[v,u] <- C[u,v]
    } 
    C[u,u] <- sum(sum(I[,,u,u])) / (2*n*(n-1))
  }
  C[t,t] <- sum(sum(I[,,t,t])) / (2*n*(n-1))
  #obtain smooth covariance
  GetEigenAnalysisResults <- utils::getFromNamespace("GetEigenAnalysisResults", "fdapace")
  if(smooth == TRUE){
    sC <- getSmoothCov(C, tgrid, "GMeanAndGCV", "gauss",n)
    #derivat the eigenfunctions and eigenvalues
    eigResult <- GetEigenAnalysisResults(sC, tgrid, optns = list(maxK = K, verbose = FALSE, FVEthreshold = 0.9999))
  }else{
    sC <- getSmoothCov(C, tgrid, "GMeanAndGCV", "gauss",n)
    #derivat the eigenfunctions and eigenvalues
    eigResult <- GetEigenAnalysisResults(C, tgrid, optns = list(maxK = K, verbose = FALSE, FVEthreshold = 0.9999))
  }
  return(list(C=C, sC = sC, tgrid = tgrid, K = K, phi = eigResult$phi, lambda = eigResult$lambda))
}

