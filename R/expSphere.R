#' Compute an exponential map for a unit hypersphere.
#' @param base A unit vector of length \eqn{m} holding the base point of the tangent space.
#' @param tg A vector of length \eqn{m} of which the exponential map is taken.
#' @return A unit vector of length \eqn{m}.
#' @export
expSphere <- function(base,tg) {
  tgNorm <- l2norm(tg)
  if (!is.na(tgNorm) & isTRUE(all.equal(tgNorm,0)) ) {
    base
  } else {
    sin(tgNorm) * tg / tgNorm + cos(tgNorm) * base
  }
}