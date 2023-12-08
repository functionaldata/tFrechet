#' Hessian \eqn{\partial^2/\partial y \partial y^\top} of the geodesic distance \eqn{\arccos(x^\top y)} on a unit hypersphere
#' @param x,y Two unit vectors.
#' @return A Hessian matrix.
#' @export
SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}