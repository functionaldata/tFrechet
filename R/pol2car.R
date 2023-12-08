#' Transform polar to Cartesian coordinates
#' @param p A vector of length \eqn{d} \eqn{(d\ge 2)} with the first element being the radius and the others being the angles,
#' where \code{p[2]} takes values in \eqn{[0,2\pi]} and \code{p[i]} takes values in \eqn{[-\pi/2,\pi/2]}, for all \eqn{i>2} if any.
#' @return A vector of length \eqn{d} holding the corresponding Cartesian coordinates 
#' \deqn{\left(r\prod_{i=1}^{d-1}\cos\theta_i, r\sin\theta_1\prod_{i=2}^{d-1}\cos\theta_i, r\sin\theta_2\prod_{i=3}^{d-1}\cos\theta_i,\dots, r\sin\theta_{d-2}\cos\theta_{d-1}, r\sin\theta_{d-1}\right),}
#' where \eqn{r} is given by \code{p[1]} and \eqn{\theta_i} is given by \code{p[i+1]} for \eqn{i=1,\dots,d-1}.
#' @examples
#' pol2car(c(1, 0, pi/4)) # should equal c(1,0,1)/sqrt(2)
#' pol2car(c(1, pi, 0)) # should equal c(-1,0,0)
#' @export
pol2car <- function(p) {
  d <- length(p)
  if (d < 2) stop("p must have at least 2 elements.")
  r <- p[1]
  theta <- p[-1]
  x <- r * sin(theta[d-1])
  if (d == 2) {
    x <- c(r * cos(theta), x)
  } else {
    for (i in 2:(d-1)) {
      x[i] <- r * prod(cos(theta[(d-i+1):(d-1)])) * sin(theta[d-i])
    }
    x[d] <- r * prod(cos(theta))
    x <- x[d:1]
  }
  return(x)
}
