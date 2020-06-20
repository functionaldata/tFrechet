#' @title Global density regression.
#' @description Global Fréchet regression for densities with respect to \eqn{L^2}-Wasserstein distance.
#' @param xin An n by p matrix or a vector of length n (if p=1) with input measurements of the predictors.
#' @param yin A matrix or list holding the sample of observations of the response. If \code{yin} is a matrix, each row holds the observations of the response corresponding to a row in \code{xin}.
#' @param hin A list holding the histograms of the response corresponding to each row in \code{xin}.
#' @param qin A matrix or list holding the quantile functions of the response. If \code{qin} is a matrix, each row holds the quantile function of the response taking values on \code{optns$qSup} corresponding to a row in \code{xin}.
#' @param xout A k by p matrix or a vector of length k (if p=1) with output measurements of the predictors. Default is \code{xin}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}.
#' @details Available control options are \code{qSup}, \code{nqSup}, \code{lower}, \code{upper}, \code{Rsquared}, \code{bwDen}, \code{nRegGrid}, \code{delta}, \code{kernelDen}, \code{infSupport}, \code{outputGrid}. \code{Rsquared} is explained as follows and see \code{\link{LocDenReg}} for the other options.
#' \describe{
#' \item{Rsquared}{A logical variable indicating whether R squared would be returned. Default is \code{FALSE}.}
#' }
#' @return A list containing the following components:
#' \item{xout}{Input \code{xout}.}
#' \item{dout}{A matrix or list holding the output densities corresponding to \code{xout}. If \code{dout} is a matrix, each row gives a density and the domain grid is given in \code{dSup}. If \code{dout} is a list, each element is a list of two components, \code{x} and \code{y}, giving the domain grid and density function values, respectively.}
#' \item{dSup}{A numeric vector giving the domain grid of \code{dout} when it is a matrix.}
#' \item{qout}{A matrix holding the quantile functions of the output densities. Each row corresponds to a value in \code{xout}.}
#' \item{qSup}{A numeric vector giving the domain grid of \code{qout}.}
#' \item{xin}{Input \code{xin}.}
#' \item{din}{Densities corresponding to the input \code{yin}, \code{hin} or \code{qin}.}
#' \item{qin}{Quantile functions corresponding to the input \code{yin}, \code{hin} or \code{qin}.}
#' \item{Rsq}{A scalar giving the R squared value if \code{optns$Rsquared = TRUE}.}
#' \item{optns}{A list of control options used.}
#'
#' @examples
#' xin = seq(0,1,0.05)
#' yin = lapply(xin, function(x) {
#'   rnorm(100, rnorm(1,x,0.005), 0.05)
#' })
#' qSup = seq(0,1,0.02)
#' xout = seq(0,1,0.25)
#' res1 <- GloDenReg(xin=xin, yin=yin, xout=xout, optns = list(qSup = qSup))
#' plot(res1)
#'\donttest{
#' hin = lapply(yin, function(y) hist(y, breaks = 50, plot=FALSE))
#' res2 <- GloDenReg(xin=xin, hin=hin, xout=xout, optns = list(qSup = qSup))
#' plot(res2)
#'}
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export

GloDenReg <- function(xin=NULL, yin=NULL, hin=NULL, qin=NULL, xout=NULL, optns=list()) {
  if (is.null(optns$Rsquared)) optns$Rsquared <- FALSE
  if (is.null(xin))
    stop ("xin has no default and must be input by users.")
  if (is.null(yin) & is.null(qin) & is.null(hin))
    stop ("One of the three arguments, yin, hin and qin, should be input by users.")
  if (is.null(xout))
    xout <- xin
  if (!is.null(optns$qSup)) {
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
  } else {
    if (!(is.null(yin) & is.null(hin))) {
      if(is.null(optns$nqSup)) {
        optns$nqSup <- 201
      }
      optns$qSup <- seq(0,1,length.out = optns$nqSup)
    } else {
      if (is.matrix(qin)) {
        optns$qSup <- seq(0,1,length.out = ncol(qin))
        warning ("optns$qSup is missing and is set by default as an equidistant grid on [0,1] with length equal to the number of columns in matrix qin.")
      } else {
        if(is.null(optns$nqSup)) {
          optns$nqSup <- 201
        }
        optns$qSup <- seq(0,1,length.out = optns$nqSup)
      }
    }
  }
  qSup <- optns$qSup

  optnsRegIdx <- match(c("Rsquared","lower","upper","qSup","nqSup"), names(optns))
  optnsRegIdx <- optnsRegIdx[!is.na(optnsRegIdx)]
  optnsReg <- optns[optnsRegIdx]

  optnsDen <- optns[-optnsRegIdx]
  if (!is.null(optnsDen$kernelDen))
    names(optnsDen)[which(names(optnsDen) == "kernelDen")] <- "kernel"
  if (!is.null(optnsDen$bwDen))
    names(optnsDen)[which(names(optnsDen) == "bwDen")] <- "userBwMu"
  # moved to just the last step transforming output quantile to densities
  # don't want a common support before F reg
  #if (!is.null(optnsDen$ndSup))
  #  names(optnsDen)[which(names(optnsDen) == "ndSup")] <- "nRegGrid"
  #if (!is.null(optnsDen$dSup))
  #  names(optnsDen)[which(names(optnsDen) == "dSup")] <- "outputGrid"

  if (!(is.null(yin) & is.null(hin))) {
    #require(fdadensity)

    if (!is.null(yin)) {
      if (!is.null(hin) | !is.null(qin))
        warning ("hin and qin are redundant when yin is available.")
      if (is.matrix(yin))
        yin <- as.data.frame(t(yin))
      if (!is.list(yin))
        stop ("yin must be a matrix or list.")
      den <- lapply(yin, CreateDensity, optns = optnsDen)
    } else if (!is.null(hin)) {
      if (!is.null(qin))
        warning ("qin is redundant when hin is available.")
      if (!is.list(hin) | length(hin) != length(xin))
        stop ("hin must be a list of the same length as xin.")
      for (histogram in hin) {
        if (!is.list(histogram))
          stop ("Each element of hin must be a list.")
        if (is.null(histogram$breaks) & is.null(histogram$mids))
          stop ("Each element of hin must be a list with at least one of the components breaks or mids.")
        if (is.null(histogram$counts))
          stop ("Each element of hin must be a list with component counts.")
      }
      den <- lapply(hin, function(histogram) {
        CreateDensity(histogram = histogram, optns = optnsDen)
      })
    }
    qin <- sapply(den, function(deni) {
      fdadensity::dens2quantile(dens = deni$y, dSup = deni$x, qSup = qSup)
    })
    qin <- t(qin)
  } else {
    #if (!is.matrix(qin))
    #  stop ("qin must be a matrix, of which each row holding the values of a quantile function evaluated on a common grid from 0 to 1.")
    if (!is.matrix(qin)) {
      if (!is.list(qin))
        stop ("qin must be a matrix or list.")
      for (qt in qin) {
        if (!is.list(qt)) {
          stop ("If qin is a list, each element must also be a list with two components, x and y.")
        } else if (is.null(qt$x) | is.null(qt$y)) {
          stop ("If qin is a list, each element must also be a list with two components, x and y.")
        }
      }
      qin <- sapply(qin, function(q) {
        approx(x = q$x, y = q$y, xout = qSup, rule = 2)$y
      })
      qin <- t(qin)
    }
    den <- apply(qin, 1, function(q) qf2pdf(qf = sort(q),prob = qSup))
  }

  if (is.null(optns$denLowerThreshold)) {
    optns$denLowerThreshold <- 0.001 * mean(qin[,ncol(qin)] - qin[,1])
  } else if (optns$denLowerThreshold) {
    if(!is.numeric(optns$denLowerThreshold) | optns$denLowerThreshold < 0)
      optns$denLowerThreshold <- 0.001 * mean(qin[,ncol(qin)] - qin[,1])
  }

  if (optns$denLowerThreshold) {
    # density thresholding from below
    if (sum(sapply(den, function(d) sum(d$y < optns$denLowerThreshold/diff(range(d$x))))) > 0) {
      den <- lapply(den, function(d) {
        lower <- optns$denLowerThreshold/diff(range(d$x))
        if (sum(d$y < lower) > 0) {
          d$y[d$y < lower] <- lower
          d$y <- d$y / pracma::trapz(d$x,d$y)
        }
        list(x=d$x, y=d$y)
      })
      qin <- sapply(den, function(deni) {
        fdadensity::dens2quantile(dens = deni$y, dSup = deni$x, qSup = qSup)
      })
      qin <- t(qin)
    }
  }

  if (sum(abs(xin - 1)) == 0) {
    # compute the Fréchet mean
    qout <- t(colMeans(qin))
  } else {
    regRes <- GloWassReg(xin = xin, qin = qin, xout = xout, optns = optnsReg)
    qout <- regRes$qout
  }

  if (!is.null(optnsDen$ndSup))
    names(optnsDen)[which(names(optnsDen) == "ndSup")] <- "nRegGrid"
  if (!is.null(optnsDen$dSup))
    names(optnsDen)[which(names(optnsDen) == "dSup")] <- "outputGrid"

  if (is.null(optnsDen$outputGrid)) {
    dout <- apply(qout, 1, qf2pdf, prob = qSup, optns = optnsDen)
    dout <- lapply(dout, function(d) d[c("x","y")])
    res <- list(xout = xout, dout = dout, qout = qout, qSup = qSup, xin=xin, din=den, qin=qin, optns=optns)
  } else {
    dSup <- optnsDen$outputGrid
    dout <- apply(qout, 1, function(q) qf2pdf(q, prob = qSup, optns = optnsDen)$y)
    #dout <- apply(qout, 1, qnt2dens, qSup = qSup, dSup = optnsDen$outputGrid)
    dout <- t(dout)
    res <- list(xout = xout, dout = dout, dSup = dSup, qout = qout, qSup = qSup, xin=xin, din=den, qin=qin, optns=optns)
  }
  if (optns$Rsquared & sum(abs(xin - 1)) > 0) res$Rsq <- regRes$R.squared
  class(res) <- "denReg"
  return(res)
}
