# gradient w.r.t. y of the geodesic distance \eqn{\arccos(x^\top y)} on a unit hypersphere
SpheGeoGrad <- function(x,y) { #, tol = 1e-10){
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}