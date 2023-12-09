#' Compute gradient w.r.t. y of the geodesic distance \eqn{\arccos(x^\top y)} on a unit hypersphere
#' @param x,y Two unit vectors.
#' @return A vector holding gradient w.r.t. \code{y} of the geodesic distance between \code{x} and \code{y}.
#' @export
SpheGeoGrad <- function(x,y) { 
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}