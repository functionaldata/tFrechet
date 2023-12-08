# using trust package and perturbation for initial value
#' @noRd
#' @import trust trust

GloSpheGeoReg <- function(xin, yin, xout) {
  k = length(xout)
  n = length(xin)
  m = ncol(yin)
  
  xbar <- colMeans(xin)
  Sigma <- cov(xin) * (n-1) / n
  invSigma <- solve(Sigma)
  
  yout = sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}
