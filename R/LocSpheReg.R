#' @title Local Fréchet Regression for Spherical Data
#' 
#' @description  Local Fréchet regression for spherical data with respect to the geodesic distance.
#' 
#' @param xin A vector of length n with input measurement points.
#' @param yin An n by m matrix holding the spherical data, of which the sum of squares of elements within each row is 1.
#' @param xout A vector of length k with output measurement points; Default: \code{xout = xin}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bw}{A scalar used as the bandwidth or \code{"CV"} (default).}
#' \item{kernel}{A character holding the type of kernel functions for local Fréchet regression for densities; \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' }
#' @return A list containing the following components:
#' \item{xout}{Input \code{xout}.}
#' \item{yout}{A k by m matrix holding the fitted responses, of which each row is a spherical vector, corresponding to each element in \code{xout}.}
#' \item{xin}{Input \code{xin}.}
#' \item{yin}{Input \code{yin}.}
#' \item{optns}{A list of control options used.}
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' # simulate the data according to the simulation in Petersen & Müller (2019)
#' xin <- runif(n)
#' err_sd <- 0.2
#' xout <- seq(0,1,length.out = 51)
#' 
#' phi_true <- acos(xin)
#' theta_true <- pi * xin
#' ytrue <- cbind(
#'   sin(phi_true) * cos(theta_true),
#'   sin(phi_true) * sin(theta_true),
#'   cos(phi_true)
#' )
#' basis <- list(
#'   b1 = cbind(
#'     cos(phi_true) * cos(theta_true),
#'     cos(phi_true) * sin(theta_true),
#'     -sin(phi_true)
#'   ),
#'   b2 = cbind(
#'     sin(theta_true),
#'     -cos(theta_true),
#'     0
#'   )
#' )
#' yin_tg <- basis$b1 * rnorm(n, mean = 0, sd = err_sd) + 
#'   basis$b2 * rnorm(n, mean = 0, sd = err_sd)
#' yin <- t(sapply(seq_len(n), function(i) {
#'   tgNorm <- sqrt(sum(yin_tg[i,]^2))
#'   if (tgNorm < 1e-10) {
#'     return(ytrue[i,])
#'   } else {
#'     return(sin(tgNorm) * yin_tg[i,] / tgNorm + cos(tgNorm) * ytrue[i,])
#'   }
#' }))
#' 
#' res <- LocSpheReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.15, kernel = "epan"))
#' }
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export 

LocSpheReg <- function(xin=NULL, yin=NULL, xout=NULL, optns=list()){
  
  if (is.null(xin))
    stop ("xin has no default and must be input by users.")
  if (is.null(yin))
    stop ("yin has no default and must be input by users.")
  if (is.null(xout))
    xout <- xin
  if (!is.vector(xin) | !is.numeric(xin))
    stop("xin should be a numerical vector.")
  if (!is.matrix(yin) | !is.numeric(yin))
    stop("yin should be a numerical matrix.")
  if (!is.vector(xout) | !is.numeric(xout))
    stop("xout should be a numerical vector.")
  if (length(xin)!=nrow(yin))
    stop("The length of xin should be the same as the number of rows in yin.")
  if (sum(abs(rowSums(yin^2) - rep(1,nrow(yin))) > 1e-6)){
    yin = yin / sqrt(rowSums(yin^2))
    warning("Each row of yin has been standardized to enforce sum of squares equal to 1.")
  }
  
  if (is.null(optns$bw)){
    optns$bw <- "CV" #max(sort(xin)[-1] - sort(xin)[-length(xin)]) * 1.2
  }
  if (is.character(optns$bw)) {
    if (optns$bw != "CV") {
      warning("Incorrect input for optns$bw.")
    }
  } else if (!is.numeric(optns$bw)) {
    stop("Mis-specified optns$bw.")
  }
  if (length(optns$bw) > 1)
    stop("bw should be of length 1.")
  
  if (is.null(optns$kernel))
    optns$kernel <- "gauss"
  
  if (is.numeric(optns$bw)) {
    bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$kernel)
    if (optns$bw < bwRange$min | optns$bw > bwRange$max) {
      optns$bw <- "CV"
      warning("optns$bw is too small or too large; reset to be chosen by CV.")
    }
  } 
  if (optns$bw == "CV") {
    optns$bw <- bwCV_sphe(xin = xin, yin = yin, xout = xout, optns = optns)
  }
  yout <- LocSpheGeoReg(xin = xin, yin = yin, xout = xout, optns = optns)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin, optns = optns)
  class(res) <- "spheReg"
  return(res)
}
