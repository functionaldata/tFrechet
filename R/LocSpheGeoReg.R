# using trust package and perturbation for initial value
#' @noRd
#' @import trust trust
#' 
LocSpheGeoReg <- function(xin, yin, xout, optns = list()) {
  k = length(xout)
  n = length(xin)
  m = ncol(yin)
  
  bw <- optns$bw
  ker <- kerFctn(optns$kernel)
  
  yout = sapply(1:k, function(j){
    mu0 = mean(ker((xout[j] - xin) / bw))
    mu1 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
    mu2 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
    s = ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
      (mu0 * mu2 - mu1^2)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if (sum(sapply(1:n, function(i) sum(yin[i,]*y0))[ker((xout[j] - xin) / bw)>0] > 1-1e-8)){
    #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) +
    #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      # return(y0)
      y0 = y0 + rnorm(3) * 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      # y <- y / l2norm(y)
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
    res = trust::trust(objFctn, y0, 0.1, 1e5)
    # res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument / l2norm(res$argument))
  })
  return(t(yout))
}
