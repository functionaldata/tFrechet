#' @title Fréchet Variance Trajectory for densities
#' @description Modeling time varying density objects with respect to $L^2$-Wasserstein distance by Fréchet variance trajectory
#' @param tgrid Time grid vector for the time varying object data.
#' @param yin An array or list of lists holding the samples of observations. If \code{yin} is an array, it has size \code{n} x \code{length(tgrid)} x numbers of samples holding the observation, such that \code{yin[i,j,]} holds the observations to the ith sample at the jth time grid. If \code{yin} is a list of lists, \code{yin[[i]][[j]]} holds the observations to the ith sample at the jth time grid.     
#' @param hin A list of lists holding the histogram for each subject. \code{hin[[i]][[j]]} holds the histogram to the ith sample at the jth time grid.     
#' @param din A three dimension array of size \code{n} x \code{length(tgrid)} x \code{length(optns$dSup)} holding the observed densities, such that \code{din[i,j,]} holds the observed density function taking values on \code{optns$dSup} corresponding to the ith sample at the jth time grid.
#' @param qin A three dimension array of size \code{n} x \code{length(tgrid)} x \code{length(optns$qSup)} holding the observed quantiles, such that \code{din[i,j,]} holds the observed density function taking values on \code{optns$qSup} corresponding to the ith sample at the jth time grid.
#' Note that only one of \code{yin}, \code{hin}, \code{din} and \code{qin} needs to be input. If more than one of them are specified, \code{yin} overwrites \code{hin}, \code{hin} overwrites \code{din} and \code{din} overwrites \code{qin}.
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
#' \item{tgridin}{Input \code{tgrid}.}
#' \item{qSup}{A vector of dimension \code{length(tgridin)} giving the domain grid of quantile functions \code{qout}.}
#' \item{qout}{A three dimension array of dimension \code{n} x \code{length(tgridin)} x \code{length(qSup)} holding the observed quantiles, such that \code{qout[i,j,]} holds the observed density function taking values on \code{qSup} corresponding to the ith sample at the jth time grid.}
#' \item{qmean}{A \code{length(tgridin)} X \code{length(qSup)} matrix containing the time varying Fréchet mean function.}
#' \item{VarTraj}{A \code{n} X \code{length(tgridin)} matrix containing the variance trajectory.}
#' @examples
#' \donttest{
#' set.seed(1)
#' #use yin 
#' tgrid = seq(1, 50, length.out = 50)
#' dSup = seq(-10, 60, length.out = 100)
#' yin = array(dim=c(30, 50, 100))
#' for(i in 1:30){
#'   yin[i,,] = t(sapply(tgrid, function(t){
#'     rnorm(100, mean = rnorm(1, mean = 1, sd = 1/t))
#'   }))
#' }
#' result1 = VarObj(tgrid, yin = yin)
#' plot(result1$phi[,1])
#' plot(result1$phi[,2])
#' yin2 = replicate(30, vector("list", 50), simplify = FALSE)
#' for(i in 1:30){
#'   for(j in 1:50){
#'     yin2[[i]][[j]] = yin[i,j,]
#'   }}
#' result1 = VarObj(tgrid, yin = yin2)
#' 
#' # use hin
#' tgrid = seq(1, 50, length.out = 50)
#' dSup = seq(-10, 60, length.out = 100)
#' hin =  replicate(30, vector("list", 50), simplify = FALSE)
#' for(i in 1:30){
#'   for (j in 1:50){
#'     hin[[i]][[j]] = hist(yin[i,j,])
#'   }
#' }
#' result2 = VarObj(tgrid, hin = hin)
#' 
#' # use din
#' tgrid = seq(1, 50, length.out = 50)
#' dSup = seq(-10, 60, length.out = 100)
#' din = array(dim=c(30, 50, 100))
#' for(i in 1:30){
#'   din[i,,] = t(sapply(tgrid, function(t){
#'     dnorm(dSup, mean = rnorm(1, mean = t, sd = 1/t))
#'   }))
#' }
#' result3 = VarObj(tgrid, din = din, optns=list(dSup = dSup))
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
#' result4 = VarObj(tgrid, qin = qin, optns=list(qSup = round(qSup, 4)))
#' }
#' @references
#' \cite{Dubey, P., & Müller, H. G. (2021). Modeling Time-Varying Random Objects and Dynamic Networks. Journal of the American Statistical Association, 1-33.}
#' @export
#' @import fdadensity
#' @import fdapace

VarObj <- function(tgrid, yin = NULL, hin = NULL, din = NULL, qin = NULL, optns=list()){
  
  ### Check input ###
  if(is.null(tgrid)){
    stop ("tgrid has no default and must be input by users.")
  }
  if (is.null(din) & is.null(qin) & is.null(yin) & is.null(hin))
    stop ("One of the four arguments, yin, hin, din and qin, should be input by users.")
  
  #check din
  if (!is.null(din)){ 
    if (!is.array(din)){
      stop ("din must be a three dimensional array ")
    }
    if (length(dim(din))!=3 || dim(din)[2] != length(tgrid)){
      stop ("din must be a three dimensional array with the dimension of second argument consisting with the length of tgrid")
    }
  }
  
  #check qin
  if (!is.null(qin)){ 
    if (!is.array(qin)){
      stop ("qin must be a three dimensional array")
    }
    if (length(dim(qin))!=3 || dim(qin)[2] != length(tgrid)){
      stop ("qin must be a three dimensional array with the dimension of second argument consisting with the length of tgrid")
    }
  }
  
  #check yin
  if (!is.null(yin)){
    if (is.list(yin)){
      length_y <- sapply(yin, function(yi) length(yi))
      if (sum(length_y != length(tgrid))){
        stop ("length of each subject of yin shoule be consistent with the length of tgrid")
      }
    }else if(is.array(yin)){
      if (length(dim(yin)) != 3 || dim(yin)[2] != length(tgrid) ){
        stop ("yin must be a three dimensional array with the dimension of second argument consisting with the length of tgrid")
      }
    }else{
      stop("input for yin should be either a three dimensional array of a list of a list ")
    }
    
  }
  
  #check hin
  if (!is.null(hin)){
    length_h <- sapply(hin, function(hi) length(hi))
    if (sum(length_h != length(tgrid))){
      stop ("length of each subject of hin shoule be consistent with the length of tgrid")
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
  }else {  #optns$qSup not available, create one!
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
  
  ### Create qin for different types of input ###
  nt <- length(tgrid)
  nq <- length(qSup)
  #create qin based on yin or hin 
  optnsRegIdx <- match(c("bwReg","kernelReg","lower","upper","qSup","nqSup","bwRange"), names(optns))
  optnsRegIdx <- optnsRegIdx[!is.na(optnsRegIdx)]
  optnsReg <- optns[optnsRegIdx]
  if (is.null(optnsReg$kernelReg))
    optnsReg$kernelReg <- "gauss"
  names(optnsReg)[which(names(optnsReg) == "kernelReg")] <- "ker"
  if (!is.null(optnsReg$bwReg))
    names(optnsReg)[which(names(optnsReg) == "bwReg")] <- "bw"
  
  optnsDen <- optns[-optnsRegIdx]
  if (!is.null(optnsDen$kernelDen))
    names(optnsDen)[which(names(optnsDen) == "kernelDen")] <- "kernel"
  if (!is.null(optnsDen$bwDen))
    names(optnsDen)[which(names(optnsDen) == "bwDen")] <- "userBwMu"
  
  if (!is.null(yin)){ 
    # use yin 
    if (!is.null(hin) || !is.null(din) || !is.null(qin)){
      warning("yin overwrites the other inputs")
    }
    # create qin based on yin
    if (is.list(yin)){
      n <- length(yin)
    }else{
      n <- dim(yin)[1]
    }
    qin <- array(dim=c(n, nt, nq))
    for(i in 1:n){
      for(j in 1:nt){
        if (is.list(yin)){
          dij <- CreateDensity(yin[[i]][[j]], optns = optnsDen)
        }else{
          dij <- CreateDensity(yin[i,j,], optns = optnsDen)
        }
        qin[i,j,] <- fdadensity::dens2quantile(dens = dij$y, dSup = dij$x, qSup = qSup)
      }
    }
  } else if (!is.null(hin)){
    if (!is.null(din) || !is.null(qin)){
      warning("hin overwrites the other inputs")
    }
    n <- length(hin)
    # create qin based on hin 
    qin <- array(dim=c(length(hin), nt, nq))
    for(i in 1:n){
      for(j in 1:nt){
        dij <- CreateDensity(histogram = hin[[i]][[j]], optns = optnsDen)
        qin[i,j,] <- fdadensity::dens2quantile(dens = dij$y, dSup = dij$x, qSup = qSup)
      }
    }
  } else if (!is.null(din)){
    if (!is.null(qin)){
      warning("din overwrites the other inputs")
    }
    # create qin based on din
    if(is.null(optns$dSup)){
      stop ("dSup need to be specified when using din")
    }
    dSup <- optns$dSup
    if (length(dSup)!=dim(din)[3]){
      stop ('dimension of the third arguement of din is not consistent with length of dSup')
    }
    n = dim(din)[1]
    qin <- array(dim = c(n, nt, nq))
    for (i in 1:n){
      for (j in 1: nt){
        qin[i,j,] <- fdadensity::dens2quantile(din[i,j,], dSup = dSup, qSup = qSup)
      }
    }
  } 
  n <- dim(qin)[1]
  

  
  ### FPCA ###
  #compute time varying fréchet mean based on qin 
  qmean <- sapply(1:nt, function(t){
    apply(qin[,t,], 2, mean)
  }) 
  qmean <- t(qmean) # dim nt x nq

  #compute variance trajectory based on qin, qmean
  VarTraj <- t(sapply(1:n, function(i){
    apply(qin[i,,] - qmean, 1, function(v){
      fdapace::trapzRcpp(qSup, v^2)
    })
  })) # dim of n x nt
  
  #FPCA based on VarTraj
  Lt <- lapply(1:n, function(i) tgrid)
  Ly <- lapply(1:n, function(i) VarTraj[i,])
  optnsFPCA <- optns[!(names(optns) %in% c("qSup", "nqSup", "dSup"))]
  VarTrajFPCA <- fdapace::FPCA(Ly = Ly, Lt = Lt, optns = optnsFPCA)
  
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

