#'@title Geodesic distance on spheres.
#'@param y1,y2 Two unit vectors, i.e., with \eqn{L^2} norm equal to 1, of the same length.
#'@return A scalar holding the geodesic distance between \code{y1} and \code{y2}.
#'@examples
#'d <- 3
#'y1 <- rnorm(d)
#'y1 <- y1 / sqrt(sum(y1^2))
#'y2 <- rnorm(d)
#'y2 <- y2 / sqrt(sum(y2^2))
#'dist <- SpheGeoDist(y1,y2)
#'@export

SpheGeoDist <- function(y1,y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
    stop("y1 is not a unit vector.")
  }
  if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}
