#' Transform Cartesian to polar coordinates
#' @param x A vector holding the Cartesian coordinates.
#' @return A vector holding the polar coordinates. See \code{\link{pol2car}} for details.
#' @examples 
#' car2pol(c(-1,0,0)) # should equal c(1, pi, 0)
#' car2pol(c(1,0,1)/sqrt(2)) # should equal c(1, 0, pi/4)
#' @noRd
car2pol <- function(x, tol = 1e-10) {
  d <- length(x)
  if (d <= 1) stop("x must have at least 2 elements.")
  r <- sqrt(sum(x^2))
  if (abs(r) < tol) {
    message("Input x has norm 0.")
    return(rep(0,d)) # angles can be any feasible values; return all zeros here
  }
  existsZeroCos <- FALSE
  theta <- asin(x[d] / r)
  if (d > 2) {
    for (i in 1:(d-2)) {
      if (prod(cos(theta)) == 0) {
        existsZeroCos <- TRUE
        theta[(i+1):(d-1)] <- 0 # can be any feasible values; return all zeros here
        break
      } else {
        tmp_sintheta <- x[d-i] / (r * prod(cos(theta)))
        if (tmp_sintheta > 1) {
          theta[i+1] <- pi/2
        } else if (tmp_sintheta < -1) {
          theta[i+1] <- -pi/2
        } else {
          theta[i+1] <- asin(tmp_sintheta)
        }
      }
    }
    theta <- theta[(d-1):1]
    if (!existsZeroCos) {
      theta[1] <- abs(theta[1])
      costheta1 <- x[1] / (r * prod(cos(theta[-1])))
      sintheta1 <- x[2] / (r * prod(cos(theta[-1])))
    }
  } else {
    costheta1 <- x[1] / r
    sintheta1 <- x[2] / r
  }
  if (sintheta1 < 0) {
    if (costheta1 > 0) {
      theta[1] <- 2 * pi - theta[1]
    } else {
      theta[1] <- theta[1] + pi
    }
  } else if (costheta1 < 0) {
    theta[1] <- pi - theta[1]
  }
  return(c(r,theta))
}
