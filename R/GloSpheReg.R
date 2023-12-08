#' @title Global Fréchet Regression for Spherical Data
#' 
#' @description  Global Fréchet regression for spherical data with respect to the geodesic distance.
#' 
#' @param xin A vector of length \eqn{n} or an \eqn{n}-by-\eqn{p} matrix with input measurement points.
#' @param yin An \eqn{n}-by-\eqn{m} matrix holding the spherical data, of which the sum of squares of elements within each row is 1.
#' @param xout A vector of length \eqn{k} or an \eqn{k}-by-\eqn{p}  with output measurement points; Default: the same grid as given in \code{xin}.
#' @return A list containing the following components:
#' \item{xout}{Input \code{xout}.}
#' \item{yout}{A \eqn{k}-by-\eqn{m} matrix holding the fitted responses, of which each row is a spherical vector, corresponding to each element in \code{xout}.}
#' \item{xin}{Input \code{xin}.}
#' \item{yin}{Input \code{yin}.}
#' 
#' @examples
#' \donttest{
#' n <- 101
#' xin <- seq(-1,1,length.out = n)
#' theta_true <- rep(pi/2,n)
#' phi_true <- (xin + 1) * pi / 4
#' ytrue <- apply( cbind( 1, phi_true, theta_true ), 1, pol2car )
#' yin <- t( ytrue )
#' xout <- xin
#' res <- GloSpheReg(xin=xin, yin=yin, xout=xout)
#' }
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export 

GloSpheReg <- function(xin=NULL, yin=NULL, xout=NULL){
  
  if (is.null(xin))
    stop ("xin has no default and must be input by users.")
  if (is.null(yin))
    stop ("yin has no default and must be input by users.")
  if (is.null(xout))
    xout <- xin
  if (!is.numeric(xin))
    stop("xin should be a numerical vector or matrix.")
  if (!is.matrix(yin) | !is.numeric(yin))
    stop("yin should be a numerical matrix.")
  if (!is.numeric(xout))
    stop("xout should be a numerical vector or matrix.")
  if(is.vector(xin)){
    xin <- as.matrix(xin)
  }
  if(is.vector(xout)){
    xout <- as.matrix(xout)
  }
  if (length(xin)!=nrow(yin))
    stop("The length of xin should be the same as the number of rows in yin.")
  if (sum(abs(rowSums(yin^2) - rep(1,nrow(yin))) > 1e-6)){
    yin = yin / sqrt(rowSums(yin^2))
    warning("Each row of yin has been standardized to enforce sum of squares equal to 1.")
  }
  
  yout <- GloSpheGeoReg(xin = xin, yin = yin, xout = xout)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin)
  class(res) <- "spheReg"
  return(res)
}
