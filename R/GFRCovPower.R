#'@title Global Fréchet regression of covariance matrices with power metric
#'@noRd
#'@description Global Fréchet regression of covariance matrices with Euclidean predictors and power metric.
#'@param x An n by p matrix of predictors.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed.
#'@param M A q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#'@param xout An m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{alpha}{Non-negative parameter from the power metric. Default is 1 which corresponds to Frobenius metric.}
#' }
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance or correlation matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#' \donttest{
#'#Example y input
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
#'Cov_est=GFRCovPower(x=x,y=y,xout=xout,optns=list(alpha=3,corrOut=FALSE))
#'#Example M input
#'n=10 #sample size
#'m=5 # dimension of covariance matrices
#'M <- array(0,c(m,m,n))
#'for (i in 1:n){
#'  y0=rnorm(m)
#'  aux<-diag(m)+y0%*%t(y0)
#'  M[,,i]<-aux
#'}
#'x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
#'xout=cbind(runif(3),runif(3)) #output predictor levels
#'Cov_est=GFRCovPower(x=x,M=M,xout=xout,optns=list(alpha=3,corrOut=FALSE))
#'}
#' @references
#' \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}


GFRCovPower  = function(x,y=NULL,M=NULL,xout,optns = list()){
  if(is.null(optns$corrOut)){
    corrOut=FALSE
  } else{
    corrOut=optns$corrOut
  }
  if(is.null(optns$alpha)){
    alpha=1
  } else{
    alpha=optns$alpha
  }

  if(!is.matrix(x)){
    stop('x must be a matrix')
  }
  if(!is.matrix(xout)){
    stop('xout must be a matrix')
  }
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have the same number of columns')
  }
  if(alpha<0){
    stop("alpha must be non-negative")
  }
  invVa = solve(var(x))
  mx = apply(x,2,mean)

  if(!is.null(y)){
    if(!is.matrix(y)){
      stop('y must be a matrix')
    }
    if(nrow(x) != nrow(y)){
      stop('x and y must have the same number of rows')
    }
    n = nrow(y)
    cm = mean4GloCovReg(x=x,y=y,xout=x)$mean_out
    #conditional covariance
    M=array(0,c(dim(y)[2], dim(y)[2], dim(y)[1]))
    for(i in 1:n){
      M[,,i] = (y[i,] - cm[i,]) %*% t(y[i,] - cm[i,])
    }
  } else{
    if(is.null(M)){
      stop("y or M must be provided.")
    }
    if(is.list(M)){
      M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
    }else{
      if(!is.array(M)){
        stop('M must be an array or a list')
      }
    }
    if(nrow(x)!=dim(M)[3]){
      stop("The number of rows of x must be the same as the number of covariance matrices in M")
    }
  }
  M_hat=array(0,c(dim(M)[1],dim(M)[1],nrow(xout)))
  n=dim(x)[1]

  if(alpha>0){
    for(j in 1:nrow(xout)){
      s = array(0,n)
      for(i in 1:n){
        s[i] = 1+(x[i,]-mx)%*%invVa%*%(xout[j,]-mx)
      }
      for(i in 1:n){
        P=eigen(M[,,i])$vectors
        Lambd_alpha=diag(eigen(M[,,i])$values**alpha)
        M_alpha=P%*%Lambd_alpha%*%t(P)
        M_hat[,,j]=M_hat[,,j]+s[i]*M_alpha/n
      }
      M_hat[,,j]=as.matrix(Matrix::nearPD(M_hat[,,j],corr = FALSE)$mat)
      P=eigen(M_hat[,,j])$vectors
      Lambd_alpha=diag(eigen(M_hat[,,j])$values**(1/alpha))
      M_hat[,,j]=P%*%Lambd_alpha%*%t(P)
      M_hat[,,j]=as.matrix(Matrix::forceSymmetric(M_hat[,,j]))
    }
  } else{
    for(j in 1:nrow(xout)){
      s = array(0,n)
      for(i in 1:n){
        s[i] = 1+(x[i,]-mx)%*%invVa%*%(xout[j,]-mx)
      }
      for(i in 1:n){
        P=eigen(M[,,i])$vectors
        Lambd_alpha=diag(log(eigen(M[,,i])$values))
        M_alpha=P%*%Lambd_alpha%*%t(P)
        M_hat[,,j]=M_hat[,,j]+s[i]*M_alpha/n
      }
      M_hat[,,j]=as.matrix(Matrix::nearPD(M_hat[,,j],corr = FALSE)$mat)
      P=eigen(M_hat[,,j])$vectors
      Lambd_alpha=diag(exp(eigen(M_hat[,,j])$values))
      M_hat[,,j]=P%*%Lambd_alpha%*%t(P)
      M_hat[,,j]=as.matrix(Matrix::forceSymmetric(M_hat[,,j]))
    }
  }
  if(corrOut){
    for(j in 1:nrow(xout)){
      D=diag(1/sqrt(diag(M_hat[,,j])))
      M_hat[,,j]=D%*%M_hat[,,j]%*%D
      M_hat[,,j]=as.matrix(Matrix::forceSymmetric(M_hat[,,j]))
    }
  }
  Mout=list()
  for(j in 1:nrow(xout)){
    Mout=c(Mout,list(M_hat[,,j]))
  }
  optns$corrOut=corrOut
  optns$alpha=alpha
  return(list(xout=xout, Mout=Mout, optns=optns))
}






