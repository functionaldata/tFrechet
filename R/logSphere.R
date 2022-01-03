# log map of a unit hypersphere.
logSphere <- function(base, x, tol = 1e-10) {
  tg <- (x - sum(x * base) * base)
  tgNorm <- l2norm(tg)
  if (!is.na(tgNorm) & tgNorm <= tol) {
    rep(0,length(base))
  } else {
    tg / l2norm(tg) * SpheGeoDist(base,x)
  }
}