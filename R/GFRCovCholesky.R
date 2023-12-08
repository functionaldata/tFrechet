#' @title Global Fréchet regression of covariance matrices with Log-Cholesky and Cholesky metric
#' @noRd
#' @description Global Fréchet regression of covariance matrices with Euclidean predictors.
#' 
#' @param x an n by p matrix of predictors.
#' @param M an q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#' @param xout an m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{metric}{Metric type choice, "log_cholesky", "cholesky" - default: \code{log_cholesky} for log Cholesky metric}
#' }
#' 
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' 
#' @examples
#' n=30 #sample size
#' m=5  #dimension of covariance matrices
#' x=cbind(matrix(rnorm(n),n),matrix(rnorm(n),n)) #vector of predictor values
#' M <- array(0,c(m,m,n))
#' a = rnorm(m); b = rnorm(m)
#' A = diag(m)+a%*%t(a);
#' B = diag(m)+3*b%*%t(b);
#' for (i in 1:n){
#'   aux <- x[i,1]*A + x[i,2]**2*B
#'   M[,,i] <- aux %*% t(aux)
#' }
#' xout=cbind(runif(5),runif(5)) #output predictor levels
#' Covlist <- GFRCovCholesky(x=x, M=M, xout=xout)
#' 
#' @references
#' \cite{A Petersen and HG Müller (2019). "Fréchet regression for random objects with Euclidean predictors." An. Stat. 47, 691-719.}
#' \cite{Z Lin (2019). " Riemannian Geometry of Symmetric Positive Definite Matrices via Cholesky Decomposition." Siam. J. Matrix. Anal, A. 40, 1353–1370.}
#' @importFrom Matrix forceSymmetric

GFRCovCholesky <- function(x, M, xout, optns = list()){
  
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(!is.matrix(xout)&!is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  if(is.vector(x)){x<- matrix(x,length(x)) }
  
  if(is.vector(xout)){xout<- matrix(xout,length(xout)) }
  
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have same number of columns')
  }

  
  if(is.null(M)){
    stop("M must be provided")
  }
  if(is.list(M)){
    M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
  }else{
    if(!is.array(M)){
      stop('M must be an array or a list')
    } else if (length(dim(M))!=3) {
      stop('M must be an array or a list')
    }
  }
  if(nrow(x)!=dim(M)[3]){
    stop("the number of rows of x must be the same as the number of covariance matrices in M")
  }
  
  if(is.null(optns$corrOut)){
    corrOut = FALSE
  } else {
    corrOut = optns$corrOut
  }
  
  if(is.null(optns$metric)){
    metric = 'log_cholesky'
  } else {
    metric =  optns$metric
  }

 
  
  n = nrow(x)
  p = ncol(x)
  nout = nrow(xout)
  
  invVa = solve(var(x))
  mx = apply(x,2,mean)
  
  MM = list()
  if(is.array(M)){
    for (i in 1:n) {
      MM[[i]] = M[,,i]
    }
  } else {MM = M}
  
  M = lapply(MM, function(X) (X+t(X))/2)
  Mout = list()
  if(metric == 'log_cholesky'){
    LL = lapply(M, chol)
    L = lapply(LL, function(X) X - diag(diag(X)))
    D = lapply(LL, function(X) diag(X))
    
    for (j in 1:nout) {
      ss = 0
      U = 0
      E = 0
      s = array(0,n)
      for (i in 1:n) {
        s[i] = 1+(x[i,]-mx)%*%invVa%*%(xout[j,]-mx)
        ss = ss + s[i]
        U = U + s[i]*L[[i]]
        E = E + s[i]*log(D[[i]])
      }
      SS = U/ss + diag(exp(E/ss))
      Mout[[j]] = t(SS)%*%SS
    }
  } else{
    L = lapply(M, chol)
    for (j in 1:nout) {
      ss = 0
      U = 0
      s = array(0,n)
      for (i in 1:n) {
        s[i] = 1+(x[i,]-mx)%*%invVa%*%(xout[j,]-mx)
        ss = ss + s[i]
        U = U + s[i]*L[[i]]
      }
      Mout[[j]] = t(U/ss) %*% (U/ss)
    }
  }
  
  if(corrOut){
    for(j in 1:nrow(xout)){
      D=diag(1/sqrt(diag(Mout[[j]])))
      Mout[[j]]=D%*%Mout[[j]]%*%D
      Mout[[j]]=as.matrix(Matrix::forceSymmetric(Mout[[j]]))
    }
  }
  out = list(xout=xout, Mout=Mout, optns=list(corrOut=corrOut,metric=metric))
  return(out)
}
  
  



