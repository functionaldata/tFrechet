#'@title \eqn{L^2} Wasserstein distance between two distributions.
#'@param d1,d2 Two lists holding the density functions of the two distributions, each consisting of two numeric vectors \code{x} and \code{y} of the same length,
#' where \code{x} holds the support grid and \code{y} holds the values of the density function.
#' Note that \code{d1$x} and \code{d2$x} can have different lengths, and \code{d1} and \code{d2} overwrite \code{q1} and \code{q2}.
#'@param q1,q2 Two numeric vectors holding the values of the quantile functions of the two distributions evaluated on the same grid given in \code{qSup}.
#'@param qSup A numeric vector holding the support grid of \code{q1} and \code{q2}; default: \code{seq(0,1,length(q1))}.
#'@return A scalar holding the \eqn{L^2} Wasserstein distance between \code{d1} and \code{d2} or between \code{q1} and \code{q2}.
#'@examples
#'d <- 3
#'y1 <- rnorm(d)
#'y1 <- y1 / sqrt(sum(y1^2))
#'y2 <- rnorm(d)
#'y2 <- y2 / sqrt(sum(y2^2))
#'dist <- SpheGeoDist(y1,y2)
#'@export

DenDist <- function(d1=NULL,d2=NULL,q1=NULL,q2=NULL,qSup=NULL) {
  if (!is.null(d1)) {
    tol <- 1e-5
    if (is.null(d2)) {
      stop("Requires the input of both d1 and d2, or the input of both q1 and q2.")
    }
    if (!is.list(d1) | !is.list(d2)) {
      stop("d1 and d2 should be lists.")
    }
    if (!all(c("x", "y") %in% names(d1))) {
      stop("d1 should consist of two elements x and y.")
    }
    if (!all(c("x", "y") %in% names(d2))) {
      stop("d2 should consist of two elements x and y.")
    }
    if (abs(length(d1$x) - length(d1$y)) > 0) {
      stop("d1$x and d1$y should have the same length.")
    }
    if (abs(length(d2$x) - length(d2$y)) > 0) {
      stop("d2$x and d2$y should have the same length.")
    }
    if (any(d1$y < 0) | abs(pracma::trapz(d1$x,d1$y) - 1) > tol) {
      stop("d1 should be a density function, i.e., it does not integrate to 1 with tolerance of ",tol,", or d1$y is not all non-negative.")
    }
    if (any(d2$y < 0) | abs(pracma::trapz(d2$x,d2$y) - 1) > tol) {
      stop("d2 should be a density function, i.e., it does not integrate to 1 with tolerance of ",tol,", or d2$y is not all non-negative.")
    }
    qSup <- seq(0,1,length.out = 201)
    q1 <- fdadensity::dens2quantile(d1$y, dSup = d1$x, qSup = qSup)
    q2 <- fdadensity::dens2quantile(d2$y, dSup = d2$x, qSup = qSup)
  } else {
    if (is.null(q1) | is.null(q2)) {
      stop("Requires the input of both d1 and d2, or the input of both q1 and q2.")
    }
    if (abs(length(q1) - length(q2)) > 0) {
      stop("q1 and q2 should q1 and q2 should have the same length")
    }
    if (is.null(qSup)) {
      qSup <- seq(0,1,length.out = length(q1))
    } else {
      if (abs(length(qSup) - length(q1)) > 0) {
        stop("q1, q2, and qSup should have the same length.")
      }
    }
  }
  sqrt(pracma::trapz(qSup, (q1 - q2)^2))
}
