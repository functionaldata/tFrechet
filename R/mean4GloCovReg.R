#'@title Global Fréchet regression of conditional mean
#'@noRd
#'@description Global Fréchet regression of conditional mean function with functional response and vector predictors, input must be regular.
#'@param x An n by p matrix of predictors.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed.
#'@param xout An m by p matrix of output predictor levels
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{mean_out}{An m by l matrix of estimated conditional means at \code{xout}.}
#' @examples
#' \donttest{
#'#Example
#'n=200             # sample size
#'t=seq(0,1,length.out=100)       # length of data
#'x = matrix(runif(n),n)
#'theta1 = theta2 = array(0,n)
#'for(i in 1:n){
#'   theta1[i] = rnorm(1,x[i],x[i]^2)
#'   theta2[i] = rnorm(1,x[i]/2,(1-x[i])^2)
#'}
#'y = matrix(0,n,length(t))
#'phi1 = sqrt(3)*t
#'phi2 = sqrt(6/5)*(1-t/2)
#'y = theta1%*%t(phi1) + theta2 %*% t(phi2)
#'xout = matrix(c(0.25,0.5,0.75),3)
#'Mean_rst = mean4GloCovReg(x,y,xout)
#'}
#' @references
#' \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}

mean4GloCovReg  = function(x, y, xout){
  if(!is.matrix(x)){
    stop('x must be a matrix')
  }
  if(!is.matrix(y)){
    stop('y must be a matrix')
  }
  if(!is.matrix(xout)){
    stop('y must be a matrix')
  }
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have the same number of columns')
  }
  if(nrow(x) != nrow(y)){
    stop('x and y must have the same number of columns')
  }
  n = nrow(y)
  invVa = solve(var(x))
  p = ncol(invVa)
  mx = apply(x,2,mean)
  nGrid = ncol(y)
  m = nrow(xout)
  cm = matrix(0,m,nGrid)

  for(j in 1:m){
    a=xout[j,]
    sL = array(0,n)
    for(i in 1:n){
      sL[i] = 1-(x[i,]-mx)%*%invVa %*% (mx-a)
    }
    for(i in 1:n){
      cm[j,] = cm[j,]+sL[i]*y[i,]/n
    }
  }
  return(list(xout=xout,mean_out = cm))
}
