#'@title Local Fréchet regression of conditional covariance matrices with Frobenius metric
#'@noRd
#'@description Local Fréchet regression of covariance matrices with Euclidean predictors and Frobenius metric.
#'@param x An n by p matrix of predictors.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed.
#'@param M A q by q by n array (resp. list) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#'@param xout An m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{bwMean}{A vector of length p holding the bandwidths for conditional mean estimation if \code{y} is provided. If \code{bwMean} is not provided, it is chosen by cross validation.}
#' \item{bwCov}{A vector of length p holding the bandwidths for conditional covariance estimation. If \code{bwCov} is not provided, it is chosen by cross validation.}
#' \item{kernel}{Name of the kernel function to be chosen from 'gauss', 'rect', 'epan', 'gausvar' and 'quar'. Default is 'gauss'.}
#' }
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance or correlation matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#' \donttest{
#'### Example y input
#'n=120             # sample size
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
#'xout=matrix(c(0.25,0.5,0.75),3)
#'Cov_rst = LFRCov(x=x,y=y,xout=xout,optns=list(corrOut=FALSE,bwMean=1.5))
#'### Example M input
#'n=30 #sample size
#'m=30 # dimension of covariance matrices
#'M <- array(0,c(m,m,n))
#'for (i in 1:n){
#'  y0=rnorm(m)
#'  aux<-15*diag(m)+y0%*%t(y0)
#'  M[,,i]<-aux
#'}
#'x=matrix(rnorm(n),n)
#'xout = matrix(c(0.25,0.5,0.75),3) #output predictor levels
#'Cov_rst=LFRCov(x=x,M=M,xout=xout,optns=list(corrOut=FALSE,bwCov=2))
#'}
#' @references
#' \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' \cite{Petersen, A., Deoni, S. and Müller, H.-G. (2019). Fréchet estimation of time-varying covariance matrices from sparse data, with application to the regional co-evolution of myelination in the developing brain. The Annals of Applied Statistics, 13(1), 393--419.}
#' @importFrom Matrix nearPD forceSymmetric
#' 
LFRCov  = function(x, y=NULL,M=NULL, xout,optns = list()){
  if(is.null(optns$corrOut)){
    corrOut=FALSE
  } else{
    corrOut=optns$corrOut
  }
  if(is.null(optns$kernel)){
    kernel = 'gauss'
  } else{
    kernel=optns$kernel
  }
  if(is.null(optns$bwMean)){
    bwMean = NA
  } else{
    bwMean=optns$bwMean
  }
  if(is.null(optns$bwCov)){
    bwCov=NA
  } else{
    bwCov=optns$bwCov
  }
  bw2=bwCov
  
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(is.vector(x)){
    x<- matrix(x,length(x))
  }
  if(is.vector(xout)){
    xout<- matrix(xout,length(xout))
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
  p = ncol(x)
  if(p > 2){
    stop("The number of dimensions of the predictor x is greater than 2.")
  }
  m = nrow(xout)

  Kern=kerFctn(kernel)
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }

  computeLFR=function(idx,x0,bw2){
    #x0 and bw2 are in R^p
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

    M_aux=array(0,c(dim(M)[1],dim(M)[1],1))
    for(i in 1:length(idx)){
      M_aux[,,1]=M_aux[,,1]+sL[i]*M[,,idx[i]]/s
    }
    M_aux[,,1]
  }

  if(!is.null(y)){
    if(!is.matrix(y)){
      stop('y must be a matrix')
    }
    if(nrow(x) != nrow(y)){
      stop('x and y must have the same number of rows')
    }
    n = nrow(y)
    nGrid = ncol(y)
    cm = mean4LocCovReg(x=x,y=y,xout=x,list(bwMean = bwMean))
    bwMean = cm$optns$bwMean
    cmh = cm$mean_out

    M=array(0,c(dim(y)[2], dim(y)[2], dim(y)[1]))
    for(i in 1:n){
      M[,,i] = (y[i,] - cmh[i,]) %*% t(y[i,] - cmh[i,])
    }
    
  } else{
    if(!is.null(M)){
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
        stop("The number of rows of x must be the same as the number of covariance matrices in M")
      }
      
    } else{
      stop("y or M must be provided.")
    }
  }
  
  #CV for bandwidth bw2 selection
  if(is.na(sum(bw2))){
    if(p==1){
      bw_choice=SetBwRange(as.vector(x), as.vector(xout), kernel)
      objF=matrix(0,nrow=20,ncol=1)
      aux1=as.matrix(seq(bw_choice$min,bw_choice$max,length.out=nrow(objF)))
      for(i in 1:length(aux1)){
        #Try-catch statement in case bandwidth is too small and produces numerical issues
        objF[i] = tryCatch({
          sum(sapply(1:dim(x)[1],function(j){
            aux=as.matrix(Matrix::nearPD(computeLFR(setdiff(1:dim(x)[1],j),x[j],aux1[i]),corr = FALSE)$mat)-M[,,j]
            sum(diag(aux%*%t(aux)))
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
      objF=matrix(0,nrow=5,ncol=5)
      aux1=seq(bw_choice1$min,bw_choice1$max,length.out=nrow(objF))
      aux2=seq(bw_choice2$min,bw_choice2$max,length.out=ncol(objF))
      for(i1 in 1:length(aux1)){
        for(i2 in 1:length(aux2)){
          #Try-catch statement in case bandwidth is too small and produces numerical issues
          objF[i1,i2] = tryCatch({
            sum(sapply(1:dim(x)[1],function(j){
              aux=as.matrix(Matrix::nearPD(computeLFR(setdiff(1:dim(x)[1],j),x[j,],c(aux1[i1],aux2[i2])),corr = FALSE)$mat)-M[,,j]
              sum(diag(aux%*%t(aux)))
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
    bw2=bwCV
  }

  Mout = list()
  if(corrOut){
    for(j in 1:m){
      x0 = xout[j,]
      aux=as.matrix(Matrix::nearPD(computeLFR(1:dim(x)[1],x0,bw2),corr=FALSE)$mat)
      D=diag(1/sqrt(diag(aux)))
      aux=D%*%aux%*%D
      aux=as.matrix(Matrix::forceSymmetric(aux))
      Mout = c(Mout,list(aux))
    }
  } else{
    for(j in 1:m){
      x0 = xout[j,]
      Mout = c(Mout,list(as.matrix(Matrix::nearPD(computeLFR(1:dim(x)[1],x0,bw2),corr = FALSE)$mat)))
    }
  }
  optns$corrOut=corrOut
  optns$kernel=kernel
  optns$bwMean=bwMean
  optns$bwCov=bw2
  return(list(xout=xout, Mout=Mout, optns=optns))
}






