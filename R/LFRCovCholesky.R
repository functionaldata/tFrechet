#'@title Local Fréchet regression of conditional covariance matrices with Log-Cholesky and Cholesky metric
#'@noRd
#'@description Local Fréchet regression of covariance matrices with Euclidean predictors.
#'
#'@param x an n by p matrix of predictors.
#'@param M an q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#'@param xout an m by p matrix of output predictor levels.
#'@param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'@details Available control options are 
#'\describe{
#'\item{bwCov}{bandwidth for conditional covariance estimation. If \code{bwCov} is not provided, it is chosen by cross validation.}
#'\item{kernel}{Name of the kernel function to be chosen from 'gauss', 'rect', 'epan', 'gausvar' and 'quar'. Default is 'gauss'.}
#'\item{corrOut}{Boolean indicating if Mout is shown as correlation or covariance matrix. Default: \code{FALSE} for only a covariance matrix.}
#'\item{metric}{Metric type choice, "log_cholesky", "cholesky" - default: \code{log_cholesky} for log Cholesky metric}
#' }
#' 
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance matrices at \code{xout}.}
#' \item{opts}{A list containing the \code{opts} parameters utilized.}
#' 
#' @examples
#' n=30 #sample size
#' m=5 # dimension of covariance matrices
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
#' Covlist = LFRCovCholesky(x,M,xout)
#'
#' @references
#' \cite{A Petersen and HG Müller (2019). "Fréchet regression for random objects with Euclidean predictors." An. Stat. 47, 691-719.}
#' \cite{Z Lin (2019). " Riemannian Geometry of Symmetric Positive Definite Matrices via Cholesky Decomposition." Siam. J. Matrix. Anal, A. 40, 1353–1370.}
#' @importFrom Matrix forceSymmetric

LFRCovCholesky <- function(x, M, xout, optns=list()){
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(!is.matrix(xout)&!is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  
  if(is.vector(x)){x<- matrix(x,length(x)) }
  if(is.vector(xout)){xout<- matrix(xout,length(xout))}
  
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have same number of columns')
  }

  
  if(is.null(optns$bwCov)){
    bwCov = NA
  } else {
    bwCov = optns$bwCov
    if(min(bwCov)<=0){
      stop("bandwidth must be positive")
    }
  }

  if(is.null(optns$kernel)){
    kernel= 'gauss'
  } else {
    kernel = optns$kernel
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
  
  p = ncol(x)
  if(p>2){
    stop("The number of dimensions of the predictor x must be at most 2")
  }
  m = nrow(xout)
  n = nrow(x)

   
  Kern=kerFctn(kernel)
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }
  if(is.null(M)){
    stop("M must be provided")
  }
  
  if(is.array(M)){
    if (length(dim(M))!=3) {
      stop('M must be an array or a list')
    }
    MM = list()
    if(is.array(M)){
      for (i in 1:dim(M)[3]) {
        MM[[i]] = M[,,i]
      }
    }
    M = lapply(MM, function(X) (X+t(X))/2)
  }else{
    if(!is.list(M)){
      stop('M must be an array or a list')
    }
    M = lapply(M, function(X) (X+t(X))/2)
  }
  
  if(nrow(x)!= length(M)){
    stop("the number of rows of x must be the same as the number of covariance matrices in M")
  }
  
  computeLFRSPD=function(idx,x0,bw2){
    #idx: index for x
    #x0 m-by-p matrix,
    #bw2 are in b-by-p
    x=as.matrix(x[idx,])
    aux=K(x-matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE),bw2)
    mu0 = mean(aux)
    mu1 = colMeans(aux*(x - matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE)))
    mu2=0
    for(i in 1:length(idx)){
      mu2 = mu2 + aux[i]*(x[i,]-x0) %*% t(x[i,]-x0)/length(idx)
    }
    sL = array(0,length(idx))
    for(i in 1:length(idx)){
      sL[i] =aux[i]*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-x0))
    }
    s = sum(sL)
    
    Mout = list()
    MM = M[idx]
    n = length(idx)
    if(metric == 'log_cholesky'){
      LL = lapply(MM, chol)
      L = lapply(LL, function(X) X - diag(diag(X)))
      D = lapply(LL, function(X) diag(X))

      U = 0
      E = 0
      for (i in 1:n) {
        U = U + sL[i]*L[[i]]
        E = E + sL[i]*log(D[[i]])
      }
      SS = U/s + diag(exp(E/s))
      Mout = t(SS)%*%SS
      
    } else {
      L = lapply(MM, chol)
      U = 0
      for (i in 1:n) {
        U = U + sL[i]*L[[i]]
      }
      Mout = t(U/s) %*% (U/s)
    }

    return(Mout)
  }
  
  distance <- function(M1, M2){
    if(metric == 'log_cholesky'){
      LL1 = chol(M1); LL2 = chol(M2)
      L1 = LL1 - diag(diag(LL1)); L2 = LL2 - diag(diag(LL2))
      D1 = diag(LL1); D2 = diag(LL2)
      L = L1 - L2; D = log(D1) - log(D2)
      res = sqrt(sum(sum(L^2))+sum(D^2))
    }else{ 
      L1 = chol(M1); L2 = chol(M2)
      L = L1 - L2;
      res = sqrt(sum(sum(L^2)))
    }
    return(res)
  }
  
  #CV for bwCov selection
  if(is.na(sum(bwCov))){
    if(p==1){
      bw_choice=SetBwRange(as.vector(x), as.vector(xout), kernel)
      objF=matrix(0,nrow=20,ncol=1)
      aux1=as.matrix(seq(bw_choice$min,bw_choice$max,length.out=nrow(objF)))
      for(i in 1:nrow(objF)){
        #Try-catch statement in case bandwidth is too small and produces numerical issues
        objF[i] = tryCatch({
          sum(sapply(1:dim(x)[1],function(j){
            distance(computeLFRSPD(setdiff(1:dim(x)[1],j),x[j],aux1[i]), M[[j]])
          }))
        }, error = function(e) {
          return(NA)
        })
      }
      if(sum(is.na(objF))==dim(objF)[1]*dim(objF)[2]){
        stop("Bandwidth too small in cross-validation search")
      }else{
        ind=which(objF==min(objF,na.rm=TRUE))[1]
        bwCV=aux1[ind]
      }
    }
    if(p==2){
      bw_choice1=SetBwRange(as.vector(x[,1]), as.vector(xout[,1]), kernel)
      bw_choice2=SetBwRange(as.vector(x[,2]), as.vector(xout[,2]), kernel)
      if(n<=30){
        objF=matrix(0,nrow=6,ncol=6)
        aux1=seq(bw_choice1$min,bw_choice1$max,length.out=nrow(objF))
        aux2=seq(bw_choice2$min,bw_choice2$max,length.out=ncol(objF))
        for(i1 in 1:nrow(objF)){
          for(i2 in 1:ncol(objF)){
            #Try-catch statement in case bandwidth is too small and produces numerical issues
            objF[i1,i2] = tryCatch({
              sum(sapply(1:dim(x)[1],function(j){
                distance(computeLFRSPD(setdiff(1:dim(x)[1],j),x[j,],c(aux1[i1],aux2[i2])), M[[j]])
              }))
            }, error = function(e) {
              return(NA)
            })
          }
        }
        if(sum(is.na(objF))==dim(objF)[1]*dim(objF)[2]){
          stop("Bandwidth too small in cross-validation search")
        }else{
          ind=which(objF==min(objF,na.rm=TRUE),arr.ind = TRUE)
          bwCV=c(aux1[ind[1]],aux2[ind[2]])
        }
      }else{
        randIndices=sample(dim(x)[1])
        groupIndices=cut(seq(1,dim(x)[1]),breaks=10,labels=FALSE)
        cv10fold_compute=function(v,leaveIn){
          distance(computeLFRSPD(leaveIn,x[v,],c(aux1[i1],aux2[i2])),M[[v]])
        }
        objF=matrix(0,nrow=6,ncol=6)
        aux1=seq(bw_choice1$min,bw_choice1$max,length.out=nrow(objF))
        aux2=seq(bw_choice2$min,bw_choice2$max,length.out=ncol(objF))
        for(i1 in 1:nrow(objF)){
          for(i2 in 1:ncol(objF)){
            #Try-catch statement in case bandwidth is too small and produces numerical issues
            objF[i1,i2] = tryCatch({
              sum(sapply(1:10,function(j){
                leaveIn=setdiff(1:(dim(x)[1]),randIndices[groupIndices==j])
                sum(sapply(randIndices[groupIndices==j],function(v){cv10fold_compute(v,leaveIn)}))
              }))
            }, error = function(e) {
              return(NA)
            })
          }
        }
        if(sum(is.na(objF))==dim(objF)[1]*dim(objF)[2]){
          stop("Bandwidth too small in cross-validation search")
        }else{
          ind=which(objF==min(objF,na.rm=TRUE),arr.ind = TRUE)
          bwCV=c(aux1[ind[1]],aux2[ind[2]])
        }
      }
    }
    bwCov=bwCV
  }
  
  
  Mout = list()
  for (j in 1:nrow(xout)) {
    Mout[[j]] = computeLFRSPD(1:dim(x)[1], xout[j,], bwCov)
  }
  
  if(corrOut){
    for(j in 1:nrow(xout)){
      D=diag(1/sqrt(diag(Mout[[j]])))
      Mout[[j]]=D%*%Mout[[j]]%*%D
      Mout[[j]]=as.matrix(Matrix::forceSymmetric(Mout[[j]]))
    }
  }
  out = list(xout=xout,Mout=Mout,optns=list(bwCov =bwCov,kernel=kernel,corrOut=corrOut,metric=metric))
  return(out)
}




