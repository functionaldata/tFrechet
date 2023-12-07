#'@noRd
#'@import pracma
#'@import fdapace
#'@importFrom utils getFromNamespace

getSmoothCov <- function(C, tgrid, method, kern, n){
  t = length(tgrid)
  t1 = tgrid
  cyy = t(as.vector(C))
  cxxn = cyy
  xygrid = pracma::meshgrid(t1)
  xx = xygrid$X
  yy = xygrid$Y
  tpairn = cbind(as.vector(xx), as.vector(yy))
  win = pracma::ones(1, length(cxxn))
  indx = c()
  count = c()
  rcov = list(tpairn = tpairn, 
              cxxn = cxxn,
              indx, indx, 
              win = win,
              cyy = cyy,
              count = count)
  GCVLwls2DV2 <- utils::getFromNamespace("GCVLwls2DV2", "fdapace")
  if (method %in% c('GCV', 'GMeanAndGCV')){
    gcvObj = GCVLwls2DV2(tgrid, tgrid, kern = kern, rcov = rcov, t = lapply(1:n, function(o) tgrid))
    bwCov <- gcvObj$h
    if (method == 'GMeanAndGCV') {
      bwCov <- sqrt(bwCov * gcvObj$minBW)
    } 
  }
  sC = fdapace::Lwls2D(bwCov, "gauss", xin = cbind(rep(tgrid, times = t), rep(tgrid, each = t)), yin = as.vector(C), xout1 = tgrid, xout2 = tgrid)
  
  return(sC)
}
