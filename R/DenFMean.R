#' @title Fréchet means of densities.
#' @description Obtain Fréchet means of densities with respect to \eqn{L^2}-Wasserstein distance.
#' @param yin A matrix or list holding the sample of observations of the response. If \code{yin} is a matrix, each row holds the observations of the response corresponding to a row in \code{xin}.
#' @param hin A list holding the histograms of the response corresponding to each row in \code{xin}.
#' @param qin A matrix or list holding the quantile functions of the response. If \code{qin} is a matrix, each row holds the quantile function of the response taking values on \code{optns$qSup} corresponding to a row in \code{xin}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}.
#' @details Available control options are \code{qSup}, \code{nqSup}, \code{lower}, \code{upper}, \code{bwDen}, \code{nRegGrid}, \code{delta}, \code{kernelDen}, \code{infSupport}, \code{outputGrid}. See \code{\link{LocDenReg}} for details.
#' @return A list containing the following components:
#' \item{dout}{A matrix or list holding the output densities corresponding to \code{xout}. If \code{dout} is a matrix, each row gives a density and the domain grid is given in \code{dSup}. If \code{dout} is a list, each element is a list of two components, \code{x} and \code{y}, giving the domain grid and density function values, respectively.}
#' \item{dSup}{A numeric vector giving the domain grid of \code{dout} when it is a matrix.}
#' \item{qout}{A matrix holding the quantile functions of the output densities. Each row corresponds to a value in \code{xout}.}
#' \item{qSup}{A numeric vector giving the domain grid of \code{qout}.}
#' \item{optns}{A list of control options used.}
#' @examples
#' xin = seq(0,1,0.05)
#' yin = lapply(xin, function(x) {
#'   rnorm(100, rnorm(1,x + x^2,0.005), 0.05)
#' })
#' res <- DenFMean(yin=yin)
#' plot(res)
#'
#' @export

DenFMean <- function(yin=NULL, hin=NULL, qin=NULL, optns=list()) {
  if (is.list(yin)) {
    xin <- rep(1, length(yin))
  } else if (is.matrix(yin)) {
    xin <- rep(1, nrow(yin))
  }
  res <- GloDenReg(xin = xin, yin = yin, hin = hin, qin = qin, optns = optns)
  if (is.list(res$dout)) {
    res$dSup <- res$dout[[1]]$x
    res$dout <- matrix(res$dout[[1]]$y, nrow = 1)
  }
  #with(res, list(dout = dout, dSup = dSup, qout = c(qout), qSup = qSup, optns = optns))
  class(res) <- "denReg"
  res
}
