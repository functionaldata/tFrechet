#'@title Plots for Fréchet regression for covariance matrices.
#' @param x A \code{covReg} object obtained from \code{\link{CovFMean}}, \code{\link{GloCovReg}} or \code{\link{LocCovReg}}.
#' @param optns A list of control options specified by \code{list(name=value)}. See 'Details'.
#' @details Available control options are
#' \describe{
#' \item{ind.xout}{A vector holding the indices of elements in \code{x$Mout} at which the plots will be made. Default is \itemize{
#' \item \code{1:length(x$Mout)} when \code{x$Mout} is of length no more than 3;
#' \item \code{c(1,round(length(x$Mout)/2),length(x$Mout))} when \code{x$Mout} is of length greater than 3.
#' }}
#' \item{nrow}{An integer --- default: 1; subsequent figures will be drawn in an \code{optns$nrow}-by-\cr
#' \code{ceiling(length(ind.xout)/optns$nrow)} array.}
#' \item{plot.type}{Character with two choices, "continuous" and "categorical".
#' The former plots the correlations in a continuous scale of colors by magnitude
#' while the latter categorizes the positive and negative entries into two different colors.
#' Default is "continuous"}
#' \item{plot.clust}{Character, the ordering method of the correlation matrix. 
#' \code{"original"} for original order (default);
#' \code{"AOE"} for the angular order of the eigenvectors;
#' \code{"FPC"} for the first principal component order;
#' \code{"hclust"} for the hierarchical clustering order, drawing 4 rectangles on the graph according to the hierarchical cluster;
#' \code{"alphabet"} for alphabetical order.}
#' \item{plot.method}{Character, the visualization method of correlation matrix to be used.
#' Currently, it supports seven methods, named "circle" (default), "square", "ellipse", "number", "pie", "shade" and "color". }
#' \item{CorrOut}{Logical, indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{plot.display}{Character, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.}
#'}
#'@return No return value.
#'@examples
#'#Example y input
#'n=20             # sample size
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
#'Cov_est=GloCovReg(x=x,y=y,xout=xout,optns=list(corrOut = FALSE, metric="power",alpha=3))
#'CreateCovRegPlot(Cov_est, optns = list(ind.xout = 2, plot.method = "shade"))
#'\donttest{
#'CreateCovRegPlot(Cov_est, optns = list(plot.method = "color"))
#'}
#'@export

CreateCovRegPlot <- function(x, optns = list()) {
  if (is.null(x$xout)) {
    ind.xout <- 1
  } else {
    if (is.null(optns$ind.xout)) {
      lMout <- length(x$Mout)
      if (lMout > 3) {
        ind.xout <- c(1,round(lMout/2),lMout)
      } else ind.xout <- seq_len(lMout)
    } else {
      if (!all(optns$ind.xout %in% (1:length(x$Mout))))
        stop("Each element in optns$ind.xout must be integers between 1 and length of x$Mout.")
      ind.xout <- optns$ind.xout
    }
  }

  if (is.null(optns$nrow)) {
    nrow <- 1
  } else {
    if (optns$nrow < 1)
      optns$nrow <- 1
    optns$nrow <- round(optns$nrow)
    nrow <- optns$nrow
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(nrow, ceiling(length(ind.xout)/nrow)))
  for (ind in 1:length(ind.xout)) {
    mout = x$Mout[[ind]]
    covplot(mout, optns = optns)
  }
}


