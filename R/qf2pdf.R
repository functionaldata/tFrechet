#' @title Converting a quantile function to a density
#' @noRd
#' @param qf a numerical vector holding the values of a quantile function at a probability grid.
#' @param prob a numerical vector holding the probability grid at which the quantile function takes values. By default, \code{prob} is an equidistant sequence with the same length as \code{qf}.
#' @param breaks a numerical vector holding the breakpoints between histogram bins. By default, \code{breaks} is an equidistant sequence from the minimum to the maximum values of \code{qf}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{userBwMu}{The bandwidth value for the smoothed mean function; positive numeric - default: determine automatically based on the data-driven bandwidth selector proposed by Sheather and Jones (1991)}
#' \item{nRegGrid}{The number of support points the KDE; numeric - default: 101.}
#' \item{delta}{The size of the bin to be used; numeric - default: determine automatically as "max(c(diff(range(y))/1000, min(diff(sort(unique(y))))))". It only works when the raw sample is available.}
#' \item{kernel}{smoothing kernel choice, "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"}
#' \item{infSupport}{logical if we expect the distribution to have infinite support or not; logical - default: TRUE}
#' \item{outputGrid}{User defined output grid for the support of the KDE, it overrides nRegGrid; numeric - default: NULL}
#' }
#'
#' @return A list containing the following fields:
#' \item{bw}{Variance for measure error.The bandwidth used by smoothing.}
#' \item{x}{A vector of length \emph{nGridReg} with the values of the KDE's support points.}
#' \item{y}{A vector of length \emph{nGridReg} with the values of the KDE at the support points.}

qf2pdf <- function(qf=NULL, prob=NULL, breaks=NULL, optns=list()){
  hist = qf2hist(qf=qf, prob=prob, breaks=breaks)
  return(CreateDensity(histogram = hist, optns=optns))
}
