#'@title \eqn{L^2} Wasserstein distance between two distributions.
#'@param d1,d2 Lists holding the density functions or quantile functions of the two distributions.
#' Each list consists of two numeric vectors \code{x} and \code{y} of the same length,
#' where \code{x} holds the support grid and \code{y} holds the values of the function.
#' Note that the type of functions representing the distributions in \code{d1} and \code{d2}
#' should be the same---either both are density functions, or both are quantile functions. 
#' If both are quantile functions, all elements in \code{d1$x} and \code{d2$x} must be between 0 and 1.
#' \code{d1$x} and \code{d2$x} may have different lengths. 
#'@param fctn_type Character vector of length 1 holding the function type in \code{d1} and \code{d2} 
#' representing the distributions: \code{"density"} (default), \code{"quantile"}.
#'@param optns A list of control parameters specified by \code{list(name=value)}.
#' @details Available control options are:
#' \describe{
#' \item{nqSup}{A scalar giving the length of the support grid of quantile functions based on which the \eqn{L^2} Wasserstein distance (i.e., the \eqn{L^2} distance between the quantile functions) is computed. Default is 201.}
#' }
#'@return A scalar holding the \eqn{L^2} Wasserstein distance between \code{d1} and \code{d2}.
#'@examples
#' d1 <- list(x = seq(-6,6,0.01))
#' d1$y <- dnorm(d1$x)
#' d2 <- list(x = d1$x + 1)
#' d2$y <- dnorm(d2$x, mean = 1)
#' dist <- dist4den(d1 = d1,d2 = d2)
#'@export
#'@importFrom fdadensity dens2quantile
#'@importFrom pracma trapz

dist4den <- function(d1 = NULL, d2 = NULL, fctn_type = NULL, optns = list()) {
  tol <- 1e-5
  if (is.null(d1) | is.null(d2)) {
    stop("Requires the input of both d1 and d2.")
  }
  if (is.null(fctn_type)) {
    fctn_type <- "density"
  }
  if (length(fctn_type) > 1) {
    fctn_type <- fctn_type[1]
    warning("fctn_type has length greater than 1---only the first element is used.")
  }
  if (!fctn_type %in% c("density","quantile")) {
    stop("Unrecognized value of fctn_type.")
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
  
  if (is.null(optns$nqSup)) {
    optns$nqSup <- 201
  }
  nqSup <- optns$nqSup
  
  if (fctn_type == "density") {
    if (any(d1$y < 0) | abs(pracma::trapz(d1$x,d1$y) - 1) > tol) {
      stop("d1 should be a density function, i.e., it does not integrate to 1 with tolerance of ",tol,", or d1$y is not all non-negative.")
    }
    if (any(d2$y < 0) | abs(pracma::trapz(d2$x,d2$y) - 1) > tol) {
      stop("d2 should be a density function, i.e., it does not integrate to 1 with tolerance of ",tol,", or d2$y is not all non-negative.")
    }
    qSup <- seq(0,1,length.out = nqSup)
    q1 <- fdadensity::dens2quantile(d1$y, dSup = d1$x, qSup = qSup)
    q2 <- fdadensity::dens2quantile(d2$y, dSup = d2$x, qSup = qSup)
  } else if (fctn_type == "quantile") {
    if (any(d1$x < 0 | d1$x > 1)) {
      stop("Some elements in d1$x do not lie in [0,1].")
    }
    if (any(d2$x < 0 | d2$x > 1)) {
      stop("Some elements in d2$x do not lie in [0,1].")
    }
    if (is.unsorted(d1$x)) {
      d1$y <- d1$y[order(d1$x)]
      d1$x <- sort(d1$x)
    }
    if (is.unsorted(d2$x)) {
      d2$y <- d2$y[order(d2$x)]
      d2$x <- sort(d2$x)
    }
    if (is.unsorted(d1$y) | is.unsorted(d2$y)) {
      stop("Quantile functions given in d1 and d2 are not monotonic.")
    } else {
      if (any(diff(d1$y)<0)) {
        len <- length(d1$x)
        d1$x <- d1$x[len:1]
        d1$y <- d1$y[len:1]
      }
      if (any(diff(d2$y)<0)) {
        len <- length(d2$x)
        d2$x <- d2$x[len:1]
        d2$y <- d2$y[len:1]
      }
    }
    
    diffSupp <- FALSE
    if (abs(length(d1$x) - length(d2$x)) > 0) {
      diffSupp <- TRUE
    } else if (sum(abs(d1$x - d2$x)) > 0) {
      diffSupp <- TRUE
    }
    if (diffSupp) {
      qSup <- seq(0,1,length.out = nqSup)
      q1 <- approx(x = d1$x, y = d1$y, xout = qSup)$y
      q2 <- approx(x = d2$x, y = d2$y, xout = qSup)$y
    } else {
      qSup <- d1$x
      q1 <- d1$y
      q2 <- d2$y
    }
  }
  sqrt(pracma::trapz(qSup, (q1 - q2)^2))
}
