#'@title Local Fréchet regression of conditional mean
#'@noRd
#'@description Local Fréchet regression of conditional mean function with functional response and vector predictors. Input must be regular.
#'@param x An n by p matrix of predictors, p should be at most 3.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed.
#'@param xout An m by p matrix of output predictor levels
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bwMean}{A vector of length p holding the bandwidths for conditional mean estimation if \code{y} is provided. If \code{bwMean} is not provided, it is chosen by cross validation.}
#' \item{kernel}{Name of the kernel function to be chosen from 'gauss', 'rect', 'epan', 'gausvar' and 'quar'. Default is 'gauss'.}
#' }
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{mean_out}{A list of estimated conditional mean vectors at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#' \donttest{
#'# Example
#'n=200             # sample size
#'t=seq(0,1,length.out=100)       # length of data
#'x = matrix(runif(n),n)
#'theta1 = theta2 = array(0,n)
#'for(i in 1:n){
#'  theta1[i] = rnorm(1,x[i],x[i]^2)
#'  theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
#'}
#'y = matrix(0,n,length(t))
#'phi1 = sqrt(3)*t
#'phi2 = sqrt(6/5)*(1-t/2)
#'y = theta1%*%t(phi1) + theta2 %*% t(phi2)
#'xout = matrix(c(0.25,0.5,0.75),3)
#'Mean_rst = mean4LocCovReg(x=x,y=y,xout=xout,list(bwMean=0.1))
#'}
#' @references
#' \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}

mean4LocCovReg = function(x, y, xout,optns = list()){
  if(is.null(optns$kernel)){
    kernel='gauss'
  } else{
    kernel=optns$kernel
  }
  if(is.null(optns$bwMean)){
    bwMean=NA
  } else{
    bwMean=optns$bwMean
  }
  if(!is.matrix(x)){
    stop('x must be a matrix')
  }
  if(!is.matrix(y)){
    stop('y must be a matrix')
  }
  if(!is.matrix(xout)){
    stop('xout must be a matrix')
  }
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have the same number of columns')
  }
  if(nrow(x) != nrow(y)){
    stop('x and y must have the same number of rows')
  }

  Kern=kerFctn(kernel)
  K = function(x,h){
    x=matrix(x,nrow=1)
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }

  n = nrow(y)
  p = ncol(x)
  m = nrow(xout)
  nGrid = ncol(y)
  cmh = matrix(0,m,nGrid)

  if(p > 3){
    warning('Local method is designed to work in low dimensional case, the result might be unstable.')
  }

  if(is.na(sum(bwMean))){
      hs = matrix(0,p,20)
      for(l in 1:p){
        hs[l,] = exp(seq(log(n^(-1/(1+p))*(max(x[,l])-min(x[,l]))/10),log(5*n^(-1/(1+p))*(max(x[,l])-min(x[,l]))),length.out =  20))
      }
      cv = array(0,20^p)
      for(k in 0:(20^p-1)){
        h = array(0,p)
        for(l in 1:p){
          kl = floor((k %% (20^l)) / (20^(l-1))) + 1
          h[l] = hs[l,kl]
        }
        yfit = matrix(0,n,nGrid)
        for(j in 1:n){
          a = x[j,]
          sL = array(0,n)
          mu0 = 0
          mu1 = 0
          mu2 = 0
          for(i in (1:n)[-j]){
            mu0 = mu0 + K(x[i,]-a,h) /n
            mu1 = mu1 + K(x[i,]-a,h)*(x[i,]-a) /n
            mu2 = mu2 + K(x[i,]-a,h)*((x[i,]-a) %*% t(x[i,]-a))/n
          }
          for(i in (1:n)[-j]){
            sL[i] =K(x[i,]-a,h)*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-a))
          }
          s = sum(sL)
          for(i in (1:n)[-j]){
            yfit[j,] = yfit[j,]+sL[i]*y[i,]/s
          }
        }
        cv[k+1] = sum((y - yfit)^2)/n/nGrid
      }
      bwi = which.min(cv)
      bwMean = array(0,p)
      for(l in 1:p){
        kl = floor((bwi %% (20^l)) / (20^(l-1))) + 1
        bwMean[l] = hs[l,kl]
      }
  }

  for(j in 1:m){
    if(length(bwMean) != p){
      stop('Dimension of bandwidth does not agree with that of Euclidean predictor X')
    }
    a = xout[j,]
    sL = array(0,n)
    mu0 = 0
    mu1 = 0
    mu2 = 0
    for(i in 1:n){
      mu0 = mu0 + K(x[i,]-a,bwMean) /n
      mu1 = mu1 + K(x[i,]-a,bwMean)*(x[i,]-a) /n
      mu2 = mu2 + K(x[i,]-a,bwMean)*(x[i,]-a) %*% t(x[i,]-a)/n
    }
    for(i in 1:n){
      sL[i] =K(x[i,]-a,bwMean)*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-a))
    }
    s = sum(sL)
    if(s == 0){
      stop('Bandwidth too small')
    }
    for(i in 1:n){
      cmh[j,] = cmh[j,]+sL[i]*y[i,]/s
    }
  }
  optns$kernel=kernel
  optns$bwMean=bwMean
  return(list(xout=xout,mean_out = cmh,optns=optns))
}








