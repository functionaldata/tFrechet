#' Compute a log map for a unit hypersphere.
#' @param base A unit vector of length \eqn{m} holding the base point of the tangent space.
#' @param x A unit vector of length \eqn{m} which the log map is taken.
#' @return A tangent vector of length \eqn{m}.
#' @export
logSphere <- function(base, x) {
  tg <- (x - sum(x * base) * base)
  tgNorm <- l2norm(tg)
  if ( !is.na(tgNorm) & isTRUE(all.equal(tgNorm, 0)) ) {
    rep( 0,length(base) )
  } else {
    tg / l2norm(tg) * SpheGeoDist(base,x)
  }
}