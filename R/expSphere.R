# exp map of a unit hypersphere.
expSphere <- function(base,tg, tol = 1e-10) {
  tgNorm <- l2norm(tg)
  if (!is.na(tgNorm) & tgNorm <= tol) {
    base
  } else {
    sin(tgNorm) * tg / tgNorm + cos(tgNorm) * base
  }
}