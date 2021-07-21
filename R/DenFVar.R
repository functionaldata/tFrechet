#'@title Fr\'{e}chet variance of densities.
#'@description Obtain Fr\'{e}chet variance of densities with respect to \eqn{L^2}-Wasserstein distance.
#'@param Ly a list of $n$ vectors containing the observed values for each density/quantile function. See `Details`.
#'@param Lx a list of $n$ vectors containing the support grid for each density/quantile function corresponding to \code{Ly}.
#'@param optns a list of control parameters specified by \code{list(name=value)}. See `Details`.
#'@details Note that the type of functions representing the distributions in \code{Ly}
#'should be the same---either all are density functions, or all are quantile functions. 
#'
#'If all support grids in \code{Lx} are the same, \code{qSup} will be set as \code{Lx[[1]]}, no matter the value of \code{nqSup}.
#'
#'Available control options are:
#'\describe{
#'\item{fctn_type}{function type in \code{Ly} representing the distributions: \code{"density"} (default), \code{"quantile"}.}
#'\item{nqSup}{a scalar giving the length of the support grid of quantile functions based on which the \eqn{L^2} Wasserstein distance (i.e., the \eqn{L^2} distance between the quantile functions) is computed. Default is 201.}
#'}
#'@return A list containing the following fields:
#'\item{DenFVar}{a scalar holding the Fr\'{e}chet variance.}
#'\item{qSup}{a numeric vector holding the support grid used.}
#'\item{optns}{the control options used.}
#'@examples
#'n <- 100
#'mu <- rnorm(n, mean = 0, sd = 0.5)
#'Ly <- lapply(1:n, function(i) qnorm(qSup, mu[i], sd = 1))
#'qSup <- seq(0.01, 0.99, (0.99-0.01)/50)
#'Lx <- rep(list(qSup), n)
#'res <- DenFVar(Ly, Lx, optns = list(fctn_type = 'quantile'))
#'res$DenFVar
#'@export

DenFVar <- function(Ly, Lx, optns = list()){
  if(is.null(Ly) | is.null(Lx)){
    stop("Requires the input of both Ly and Lx.")
  }
  if(!is.list(Ly) | !is.list(Lx)){
    stop("Ly and Lx should be lists.")
  }
  if(length(Ly)!=length(Lx)){
    stop('Ly and Lx should have the same length.')
  }
  if(sum(sapply(Ly, length)-sapply(Lx, length))){
    stop('Each vector in Ly and its corresponding vector in Lx should have the same length.')
  }
  if(is.null(optns$fctn_type)){
    fctn_type <- 'density'
  }
  else{
    fctn_type <- optns$fctn_type
  }
  if(length(fctn_type) > 1){
    fctn_type <- fctn_type[1]
    warning("fctn_type has length greater than 1---only the first element is used.")
  }
  if(!fctn_type %in% c("density","quantile")){
    stop("Unrecognized value of fctn_type.")
  }
  if(is.null(optns$nqSup)){
    nqSup <- 201
  }
  else{
    nqSup <- optns$nqSup
  }
  n <- length(Ly)
  tol <- 1e-5
  if(fctn_type == "density"){
    if(any(sapply(1:n, function(i) {
      abs(pracma::trapz(Lx[[i]], Ly[[i]]) - 1) > tol
    }))){
      stop("Each element of Ly should be a density function (integrates to $1$ with tolerance of ",tol,").")
    }
    if(any(sapply(1:n, function(i) {
      any(Ly[[i]] < 0)
    }))){
      stop("Each vector in Ly should be all non-negative.")
    }
    qSup <- seq(0, 1, length.out = nqSup)
    Ly <- lapply(1:n, function(i) fdadensity::dens2quantile(Ly[[i]], dSup = Lx[[i]], qSup = qSup))
  }
  else{
    if(any(sapply(Lx, function(Lxi) {
      any(Lxi < 0 | Lxi > 1)
    }))){
      stop("Each vector in Lx should lie in [0, 1].")
    }
    for(i in 1:n){
      if (is.unsorted(Lx[[i]])) {
        Ly[[i]] <- Ly[[i]][order(Lx[[i]])]
        Lx[[i]] <- sort(Lx[[i]])
      }
    }
    if(any(sapply(1:n, function(i) is.unsorted(Ly[[i]])))){
      stop("Quantile functions given in Ly are not monotonic.")
    }
    diffSupp <- TRUE
    if(length(unique(sapply(Lx, length)))==1){
      if(sum(diff(matrix(unlist(Lx), nrow = n, byrow = TRUE)))==0) diffSupp <- FALSE
    }
    if(diffSupp){
      qSup <- seq(0, 1, length.out = nqSup)
      Ly <- lapply(1:n, function(i) approx(x = Lx[[i]], y = Ly[[i]], xout = qSup)$y)
    }
    else{
      qSup <- Lx[[1]]
    }
  }
  DenFMean <- rowMeans(matrix(unlist(Ly), nrow = length(Ly[[1]]), ncol = n))
  DenFVar <- mean(sapply(Ly, function(Lyi) {
    pracma::trapz(qSup, (Lyi-DenFMean)^2)
  }))
  res <- list(DenFVar = DenFVar, qSup = qSup, optns = optns)
  res
}