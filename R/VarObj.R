#' @title Fréchet Variance Trajectory for densities
#' @description Modeling time varying density objects with respect to $L^2$-Wasserstein distance by Fréchet variance trajectory
#' @param tgrid Time grid vector for the time varying object data.
#' @param din A three dimension array of size \code{n} x \code{length(tgrid)} x \code{\length(optns$dSup)} holding the observed densities, such that \code{din[i,j,]} holds the observed density function taking values on \code{optns$dSup} corresponding to the ith sample at the jth time grid.
#' @param qin A three dimension array of size \code{n} x \code{length(tgrid)} x \code{\length(optns$qSup)} holding the observed quantiles, such that \code{din[i,j,]} holds the observed density function taking values on \code{optns$qSup} corresponding to the ith sample at the jth time grid.
#' Note that only one of \code{din} and \code{qin} needs to be input. If more than one of them are specified, \code{qin} overwrites \code{din}.
#' where each row holds the observations for one subject on the common grid \code{tGrid}.
#' @param optns a list of options control parameters specified by \code{list(name=value)}.
#' @details Available control options are \code{qSup}, \code{nqSup}, \code{dSup} and other options in FPCA of fdapace.
#' @return A list of the following:
#' \item{tgridout}{Time grid vector for the output time varying object data.}
#' \item{K}{Numbers of principal components.}
#' \item{nu}{A vector of dimension \code{length(tgridout)} giving the mean function support on \code{tgridout} of the Fréchet variance function.}
#' \item{lambda}{A vector of dimension \code{K} containing eigenvalues.}
#' \item{phi}{A \code{length(tgridout)} X \code{K} matrix containing eigenfunctions support on \code{tgridout} of the Fréchet variance function.}
#' \item{xiEst}{A \code{n} X \code{K} matrix containing the FPC estimates.}
#' \item{cumFVE}{A vector of dimension \code{K} with the fraction of the cumulative total variance explained with each additional FPC.}
#' \item{FPCAObj}{FPCA Object of Fréchet variance function.}
#' \item{tgridin}{Input tgrid.}
#' \item{qSup}{A vector of dimension \code{length(tgridin)} giving the domain grid of quantile functions \code{qout}.}
#' \item{qout}{A three dimension array of dimension \code{n} x \code{length(tgridin)} x \code{length(qSup)} holding the observed quantiles, such that \code{qout[i,j,]} holds the observed density function taking values on \code{qSup} corresponding to the ith sample at the jth time grid.}
#' \item{qmean}{A \code{length(tgridin)} X \code{length(qSup)} matrix containing the time varying Fréchet mean function.}
#' \item{VarTraj}{A \code{n} X \code{length(tgridin)} matrix containing the variance trajectory.}

#' @examples
#' # use din
#' tgrid = seq(1, 50, length.out = 50)
#' dSup = seq(-10, 60, length.out = 100)
#' din = array(dim=c(30, 50, 100))
#' for(i in 1:30){
#'   din[i,,] = t(sapply(tgrid, function(t){
#'     dnorm(dSup, mean = rnorm(1, mean = t, sd = 1/t))
#'   }))
#' }
#' result = VarObj(tgrid, din = din, optns=list(dSup = dSup))
#' 
#' # use qin
#' tgrid = seq(1, 50, length.out = 50)
#' qSup = seq(0.00001,1-0.00001,length.out = 100)
#' qin = array(dim=c(30, 50, 100))
#' for(i in 1:30){
#'   qin[i,,] = t(sapply(tgrid, function(t){
#'     qnorm(qSup, mean = rnorm(1, mean = t, sd = 1/t))
#'   }))
#' }
#' result = VarObj(tgrid, qin = qin, optns=list(qSup = round(qSup, 4)))
#' 
#' 
#' @references
#' \cite{Dubey, P., & Müller, H. G. (2021). Modeling Time-Varying Random Objects and Dynamic Networks. Journal of the American Statistical Association, 1-33.}
#' @export

VarObj <- function(tgrid, din = NULL, qin = NULL, optns=list()){
  require(fdadensity)
  require(fdapace)
  if(is.null(tgrid)){
    stop ("tgrid has no default and must be input by users.")
  }
  if (is.null(din) & is.null(qin))
    stop ("One of the two arguments, din and qin, should be input by users.")
  if (!is.null(din)){ #check din
    if (!is.array(din)){
      stop ("din must be a three dimensional array ")
    }
    if (length(dim(din))!=3 || dim(din)[2] != length(tgrid)){
      stop ("din must be a three dimensional array with the dimension of second argument consisting with the length of tgrid")
    }
  }else{ #check qin
    if (!is.array(qin)){
      stop ("qin must be a three dimensional array")
    }
    if (length(dim(qin))!=3 || dim(qin)[2] != length(tgrid)){
      stop ("din must be a three dimensional array with the dimension of second argument consisting with the length of tgrid")
    }
  }

  
  #Specify qSup #modified from GloDenReg.R
  if (!is.null(optns$qSup)) {  #optns$qSup available
    if (min(optns$qSup) != 0 | max(optns$qSup) - 1 != 0)
      stop ("optns$qSup must have minimum 0 and maximum 1.")
    if (sum(duplicated(optns$qSup)) > 0) {
      optns$qSup <- unique(optns$qSup)
      warning ("optns$qSup has duplicated elements which has been removed.")
    }
    if (is.unsorted(optns$qSup)) {
      optns$qSup <- sort(optns$qSup)
      warning ("optns$qSup has been reordered to be increasing.")
    }
  } else {  #optns$qSup not available, create one!
    if (!is.null(qin)){ #qin available
      optns$qSup <- seq(0,1,length.out = dim(qin)[3])
      warning ("optns$qSup is missing and is set by default as an equidistant grid on [0,1] with length equal to the number of columns in matrix qin.")
    } else{ #din available, check optns$nqSup
      if(is.null(optns$nqSup)) {
        optns$nqSup <- 201
      }
      optns$qSup <- seq(0,1,length.out = optns$nqSup)
    }
  }
  qSup <- optns$qSup
  
  #generate qin based on din 
  nt = length(tgrid)
  nq = length(qSup)
  if(!is.null(din)){
    if(!is.null(qin)){
      warning ("qin overwrites din")
    }else{ #use input din
      n = dim(din)[1]
      if(is.null(optns$dSup)){
        stop ("dSup need to be specified when using din")
      }
      dSup = optns$dSup
      if (length(dSup)!=dim(din)[3]){
        stop ('dimension of the third arguement of din is not consistent with length of dSup')
      }
      qin = array(dim = c(n, nt, length(qSup)))
      for (i in 1:n){
        for (j in 1: nt){
          qin[i,j,] = fdadensity::dens2quantile(din[i,j,], dSup = dSup, qSup = qSup)
        }
      }
    }
  }else{ #use input qin
    n = dim(qin)[1]
  }
  
  
  #compute time varying fréchet mean based on qin 
  qmean = sapply(1:nt, function(t){
    apply(qin[,t,], 2, mean)
  }) 
  qmean = t(qmean) # dim nt x nq

  #compute variance trajectory based on qin, qmean
  VarTraj = t(sapply(1:n, function(i){
    apply(qin[i,,] - qmean, 1, function(v){
      fdapace::trapzRcpp(qSup, v^2)
    })
  })) # dim of n x nt
  
  #FPCA based on VarTraj
  Lt = lapply(1:n, function(i) tgrid)
  Ly = lapply(1:n, function(i) VarTraj[i,])
  optnsFPCA = optns[!(names(optns) %in% c("qSup", "nqSup", "dSup"))]
  VarTrajFPCA = fdapace::FPCA(Ly = Ly, Lt = Lt, optns = optnsFPCA)
  
  return(list(tgridout = VarTrajFPCA$workGrid,
              nu = VarTrajFPCA$mu,
              lambda = VarTrajFPCA$lambda,
              phi = VarTrajFPCA$phi,
              xiEst = VarTrajFPCA$xiEst,
              K = VarTrajFPCA$selectK,
              cumFVE = VarTrajFPCA$cumFVE,
              FPCAObj = VarTrajFPCA,
              tgridin = tgrid,
              qSup = qSup,
              qout = qin,
              qmean = qmean,
              VarTraj = VarTraj
              ))
}

