#' Generate a "natural" frame (orthonormal basis)
#' @description Generate a "natural" frame (orthonormal basis) for the tangent space at \code{x} on the unit sphere.
#' @param x A unit vector of length \eqn{d}.
#' @return A \eqn{d}-by-\eqn{(d-1)} matrix where columns hold the orthonormal basis of the tangent space at \code{x} on the unit sphere.
#' @details The first \eqn{(i+1)} elements of the \eqn{i}th basis vector are given by 
#' \eqn{\sin\theta_i\prod_{j=1}^{i-1}\cos\theta_j}, \eqn{\sin\theta_i\sin\theta_1 \prod_{j=2}^{i-1}\cos\theta_j},
#' \eqn{\sin\theta_i\sin\theta_2 \prod_{j=3}^{i-1}\cos\theta_j}, \eqn{\dots}, \eqn{\sin\theta_i\sin\theta_{i-1}}, \eqn{-\cos\theta_i}, respectively.
#' The rest elements (if any) of the \eqn{i}th basis vector are all zero.
#' @examples
#' frameSphere(c(1,0,0,0))
#' @export
frameSphere <- function(x) {
  theta <- car2pol(x)[-1]
  d <- length(x)
  frm <- matrix(numeric(d), ncol = 1)
  frm[1:2,1] <- c(sin(theta[1]), -cos(theta[1]))
  if (d == 2) return(frm)
  
  # product of (sin(theta[1]), cos(theta[2]), ..., cos(theta[l-1]), sin(theta[l]))
  sinfirstlast <- function(theta) {
    len <- length(theta)
    if (len < 2) stop("theta must have at least two elements.")
    res <- prod(sin(theta[c(1,len)]))
    if (len > 2) res <- res * prod(cos(theta[-c(1,len)]))
    return(res)
  }
  
  frm <- cbind(
    frm,
    sapply(2:(d-1), function(j) {
      tmp <- rep(0,d)
      tmp[1] <- sin(theta[j]) * prod(cos(theta[1:(j-1)]))
      tmp[2:j] <- sapply(1:(j-1), function(k) {
        sinfirstlast(theta[k:j])
      })
      tmp[j+1] <- -cos(theta[j])
      return(tmp)
    })
  )
  return(frm)
}
