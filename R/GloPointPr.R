#' @title Global Cox point process regression.
#' @description Global Fréchet regression for replicated Cox point processes with respect to \eqn{L^2}-Wasserstein distance on shape space and Euclidean 2-norm on intensity factor space.
#' @param xin An n by p matrix with input measurements of the predictors.
#' @param tin A list holding the sample of event times of each replicated point process, where the ith element of the list \code{tin} holds the event times of the point process corresponding to the ith row of \code{xin}.
#' @param T0 A positive scalar that defines the time window [0,\code{T0}] where the replicated Cox point processes are observed.
#' @param xout A k by p matrix with output measurements of the predictors. Default is \code{xin}.
#' @param optns A list of control parameters specified by \code{list(name=value)}.
#' @details Available control options are \code{bwDen} (see \code{\link{LocDenReg}} for this option description) and
#' \describe{
#' \item{L}{Upper Lipschitz constant for quantile space; numeric -default: 1e10.}
#' \item{M}{Lower Lipschitz constant for quantile space; numeric -default: 1e-10.}
#' \item{dSup}{User defined output grid for the support of kernel density estimation used in \code{CreateDensity()} for mapping from quantile space to shape space. This grid must be in [0,\code{T0}]. Default is an equidistant grid with \code{nqSup}+2 points.}
#' \item{nqSup}{A scalar with the number of equidistant points in (0,1) used to obtain the empirical quantile function from each point process. Default: 500.}
#' }
#' @return A list containing the following components:
#' \item{xout}{Input \code{xout}.}
#' \item{dSup}{Support of each estimated (up to a constant) conditional intensity regression function in the columns of \code{intensityReg}.}
#' \item{intensityReg}{A matrix of dimension \code{length(dSup)} by \code{nrow(xout)} holding the estimated intensity regression functions up to a constant over the support grid \code{dSup}, where each column corresponds to a predictor level in the corresponding row of \code{xout}.}
#' \item{xin}{Input \code{xin}.}
#' \item{optns}{A list of control options used.}
#' @examples
#' \donttest{
#' n=100
#' alpha_n=sqrt(n)
#' alpha1=2.0
#' beta1=1.0
#' gridQ=seq(0,1,length.out=500+2)[2:(500+1)]
#' X=runif(n,0,1)#p=1
#' tau=matrix(0,nrow=n,ncol=1)
#' for(i in 1:n){
#'   tau[i]=alpha1+beta1*X[i]+truncnorm::rtruncnorm(1, a=-0.3, b=0.3, mean = 0, sd = 1.0)
#' }
#' Ni_n=matrix(0,nrow=n,ncol=1)
#' u0=0.4
#' u1=0.5
#' u2=0.05
#' u3=-0.01
#' tin=list()
#' for(i in 1:n){
#'   Ni_n[i]=rpois(1,alpha_n*tau[i])
#'   mu_x=u0+u1*X[i]+truncnorm::rtruncnorm(1,a=-0.1,b=0.1,mean=0,sd=1)
#'   sd_x=u2+u3*X[i]+truncnorm::rtruncnorm(1,a=-0.02,b=0.02,mean=0,sd=0.5)
#'   if(Ni_n[i]==0){
#'     tin[[i]]=c()
#'   }else{
#'     tin[[i]]=truncnorm::rtruncnorm(Ni_n[i],a=0,b=1,mean=mu_x,sd=sd_x) #Sample from truncated normal
#'   }
#' }
#' res=GloPointPrReg(
#'   xin=matrix(X,ncol=1),tin=tin,
#'   T0=1,xout=matrix(seq(0,1,length.out=10),ncol=1),
#'   optns=list(bwDen=0.1)
#' )
#' }
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' 
#' \cite{Gajardo, Á. and Müller, H.-G. (2022). "Cox Point Process Regression." IEEE Transactions on Information Theory, 68(2), 1133-1156.}
#' @importFrom quadprog solve.QP
#' @export

GloPointPrReg <- function(xin=NULL, tin=NULL,T0=NULL, xout=NULL, optns=list()) {
  
  if(is.null(xin)){
    stop("xin has no default and must be input")
  }
  if(!is.matrix(xin)){
    stop("xin must be a matrix")
  }
  if(is.null(xout)){
    xout=xin
  }
  if(!is.matrix(xout)){
    stop("xout must be a matrix")
  }
  if(ncol(xin)!=ncol(xout)){
    stop("xin and xout must have same number of columns")
  }
  if(is.null(tin)){
    stop("tin has no default and must be input")
  }
  if(is.null(T0)){
    stop("T0 has no default and must be input")
  }
  if(T0<max(unlist(tin))){
    stop("T0 cannot be smaller than any event time")
  }
  if(is.null(optns$L)){
    optns$L=1e10
  }
  if(is.null(optns$M)){
    optns$M=1e-10
  }
  if(optns$L<0 | optns$M<0 | optns$L<optns$M){
    stop("L and M must be positive with L>M")
  }
  if(is.null(xout)){
    xout <- xin
  }
  if(is.vector(xin)){
    xin = as.matrix(xin)
  }
  
  if(!is.null(optns$bwDen)){
    if(sum(!is.numeric(optns$bwDen) | !length(optns$bwDen)==1 | !(optns$bwDen>0))>0){
      stop("bwDen option must be a positive scalar")
    }
    optns$userBwMu=optns$bwDen
  }
  if(is.null(optns$nqSup)){
    optns$nqSup=500
  }else{
    if(sum(!is.numeric(optns$nqSup) | !length(optns$nqSup)==1 | !(optns$nqSup>0))>0){
      stop("nqSup must be a positive integer")
    }
    if(optns$nqSup<100){
      warning("nqSup option may be too small")
    }
  }
  
  if(!is.null(optns$dSup)){
    if(optns$dSup[1]!=0 | optns$dSup[length(optns$dSup)]!=T0){
      stop("dSup must be a vector with endpoints 0 and T0")
    }
  }
  
  n = nrow(xin)
  
  getNplus=function(x_predict){
    cov_X=var(xin)*(n-1)/n
    aux=xin-mean(xin)
    N_plus=max((t(Ni_n)%*%(1+aux%*%solve(cov_X)%*%(x_predict-mean(xin)))/n)/mean(Ni_n),0)
    return(N_plus)
  }
  
  gridQ=seq(0,1,length.out=optns$nqSup+2)[2:(optns$nqSup+1)]
  Qi_hat=matrix(0,nrow=n,ncol=length(gridQ))
  Ni_n=matrix(0,nrow=n,ncol=1)
  for(i in 1:n){
    Ni_n[i]=length(tin[[i]])
    if(Ni_n[i]==0){
      Qi_hat[i,]=gridQ #uniform quantile distribution function if there are no points
    }else{
      Qi_hat[i,]=as.vector(quantile(tin[[i]],gridQ))
    }
  }
  
  k = nrow(xout)
  m = ncol(Qi_hat)
  p=ncol(xin)
  xbar = colMeans(xin)
  Sigma = cov(xin) * (n-1) / n
  invSigma = solve(Sigma)
  
  A = cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  A=A[,c(1,ncol(A),2:(ncol(A)-1))]
  A=cbind(A,-A[,3:ncol(A)])
  A=cbind(c(-1,rep(0,m-1)),c(rep(0,m-1),1),c(rep(0,m-1),-1),A)
  b0 = c(-optns$L/(m+1),1-optns$L/(m+1),optns$M/(m+1)-1,optns$M/(m+1),-T0,rep(optns$M/(m+1),m-1))
  b0=c(b0,-rep(optns$L/(m+1),m-1))
  
  quantile_GFR = t(sapply(1:k, function(j){
    s = 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s = as.vector(s)/n
    gx = (s %*% Qi_hat)
    res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
    return(sort(res$solution))
  })) #each row contains the GFR quantile regression function on shape space at xout[j,]
  
  
  if(is.null(optns$dSup)){
    density_grid=as.vector(c(0,gridQ,T0))
  }else{
    density_grid=optns$dSup
    optns$dSup=NULL
  }
  
  intensityReg=sapply(1:k,function(j){
    qf2pdf(qf=as.vector(c(0,quantile_GFR[j,],T0)),prob=as.vector(c(0,gridQ,1)),optns=list(outputGrid=density_grid,infSupport=FALSE,userBwMu=optns$bwDen))$y*getNplus(xout[j,])
  })
  
  res <- list(xout = xout, dSup = density_grid, intensityReg = intensityReg, xin=xin, optns=optns)
  
  class(res) <- "PointPrReg"
  return(res)
}

