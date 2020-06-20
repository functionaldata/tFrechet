#' @title Local density regression.
#' @description Local Fréchet regression for densities with respect to \eqn{L^2}-Wasserstein distance.
#' @param xin A numeric vector holding the observations of the predictor.
#' @param yin A matrix or list holding the sample of observations of the response. If \code{yin} is a matrix, each row holds the observations of the response corresponding to a predictor value in \code{xin}.
#' @param hin A list holding the histograms of the response corresponding to each predictor value in \code{xin}.
#' @param qin A matrix or list holding the quantile functions of the response. If \code{qin} is a matrix, each row holds the quantile function corresponding to a predictor value in \code{xin}.
#' @param xout A numeric vector holding the output grid for the predictor. Default is \code{xin}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bwReg}{A scalar used as the bandwidth or a character indicating the selection method: \code{"CV"} (default) or \code{"GCV"}.}
#' \item{kernelReg}{A character holding the type of kernel functions for local Fréchet regression for densities; \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' \item{qSup}{A numeric vector holding the grid on [0,1] quantile functions take value on. Default is an equidistant grid.}
#' \item{nqSup}{A scalar giving the length of \code{qSup}. Default is 201.}
#' \item{lower}{A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.}
#' \item{upper}{A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.}
#' \item{bwRange}{Bandwidth selection range if \code{bwReg} equals \code{"CV"} or \code{"GCV"}.}
#' \item{bwDen}{The bandwidth value for the smoothed mean function; positive numeric - default: determine automatically based on the data-driven bandwidth selector proposed by Sheather and Jones (1991).}
#' \item{ndSup}{The number of support points the kernel density estimation; numeric - default: 101.}
#' \item{dSup}{User defined output grid for the support of kernel density estimation, it overrides \code{nRegGrid}; numeric - default: \code{NULL}}
#' \item{delta}{The size of the bin to be used; numeric - default: \code{diff(range(y))/1000}. It only works when the raw sample is available.}
#' \item{kernelDen}{A character holding the type of kernel functions for density estimation; \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' \item{infSupport}{logical if we expect the distribution to have infinite support or not; logical - default: \code{TRUE}}
#' \item{denLowerThreshold}{\code{FALSE} or a positive value giving the lower threshold of the densities; default: \code{0.001 * mean(qin[,ncol(qin)] - qin[,1])}.}
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
#' \item{optns}{A list of control options used.}
#'
#' @examples
#' xin = seq(0,1,0.05)
#' yin = lapply(xin, function(x) {
#'   rnorm(100, rnorm(1,x + x^2,0.005), 0.05)
#' })
#' qSup = seq(0,1,0.02)
#' xout = seq(0,1,0.2)
#' res1 <- LocDenReg(xin=xin, yin=yin, xout=xout, optns = list(bw = 0.12, qSup = qSup))
#' plot(res1)
#'\donttest{
#' xout <- xin
#' hin = lapply(yin, function(y) hist(y, breaks = 50))
#' res2 <- LocDenReg(xin=xin, hin=hin, xout=xout, optns = list(qSup = qSup))
#' plot(res2)
#'}
#' @references
#' \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' @export

LocDenReg <- function(xin=NULL, yin=NULL, hin=NULL, qin=NULL, xout=NULL, optns=list()) {
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

  if (!is.null(optns$bwReg)) {
    if (length(optns$bwReg) > 1) {
      warning ("optns$bwReg has length larger than 1 and only the first element is used.")
      optns$bwReg <- optns$bwReg[1]
    }
  }

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

  if (is.null(optnsReg$bw)) {
    optnsReg$bw <- "CV"
  }
  if (!is.numeric(optnsReg$bw)) {
    if (optnsReg$bw == "CV") {
      optnsReg$bw <- bwCV(xin=xin, qin=qin, xout=xout, optns=optnsReg)
    } else {
      if (optnsReg$bw != "GCV")
        warning("optns$bwReg was mis-specified and is reset to be chosen by GCV.")
      optnsReg$bw <- bwGCV(xin=xin, qin=qin, xout=xout, optns=optnsReg)
    }
    #if (optnsReg$bw < diff(range(xin))*0.1)
    #  optnsReg$bw <- diff(range(xin))*0.1
  } else {
    if (optnsReg$bw < max(diff(sort(xin))) & !is.null(optnsReg$ker)) {
      if(optnsReg$ker %in% c("rect","quar","epan")) {
        warning("optns$bwReg was set too small and is reset to be chosen by GCV.")
        optnsReg$bw <- bwGCV(xin=xin, qin=qin, xout=xout, optns=optnsReg)
      }
    }
  }
  optns$bwReg <- optnsReg$bw

  qout <- LocWassReg(xin = xin, qin = qin, xout = xout, optns = optnsReg)

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
  class(res) <- "denReg"
  return(res)
}

# set up bandwidth range
SetBwRange <- function(xin, xout, kernel_type) {
  xinSt <- sort(xin)
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) - xinSt[length(xin)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) * ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin))/3 / (ifelse(kernel_type == "gauss", 3, 1) * ifelse(kernel_type == "gausvar", 2.5, 1))
  if (bw.max < bw.min) {
    if (bw.min > bw.max*3/2) {
      #warning("Data is too sparse.")
      bw.max <- bw.min*1.01
    } else bw.max <- bw.max*3/2
  }
  return(list(min=bw.min, max = bw.max))
}

# bandwidth selection via cross validation
bwCV <- function(xin, qin, xout, optns) {
  qin <- qin[order(xin),]
  xin <- sort(xin)
  compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)

  # k-fold
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    folds <- numeric(length(xin))
    n <- sum(compareRange)
    numFolds <- ifelse(n > 30, 10, sum(compareRange))

    tmp <- c(sapply(1:ceiling(n/numFolds), function(i)
      sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
    tmp <- tmp[1:n]
    repIdx <- which(diff(tmp) == 0)
    for (i in which(diff(tmp) == 0)) {
      s <- tmp[i]
      tmp[i] <- tmp[i-1]
      tmp[i-1] <- s
    }
    #tmp <- cut(1:n,breaks = seq(0,n,length.out = numFolds+1), labels=FALSE)
    #tmp <- tmp[sample(seq_len(n), n)]

    folds[compareRange] <- tmp

    qfit <- lapply(seq_len(numFolds), function(foldidx) {
      testidx <- which(folds == foldidx)
      res <- LocWassReg(xin = xin[-testidx], qin = qin[-testidx,], xout = xin[testidx],
                        optns = optns1)
      res # each row is a qt function
    })
    qfit <- do.call(rbind, qfit)
    mean(apply((qfit - qin[which(compareRange)[order(tmp)],])^2, 1, pracma::trapz, x = optns1$qSup))
  }
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
  if (!is.null(optns$bwRange)) {
    if (min(optns$bwRange) < bwRange$min) {
      message("Minimum bandwidth is too small and has been reset.")
    } else bwRange$min <- min(optns$bwRange)
    if (max(optns$bwRange) >  bwRange$min) {
      bwRange$max <- max(optns$bwRange)
    } else {
      message("Maximum bandwidth is too small and has been reset.")
    }
  }
  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  res$minimum
}
if(0){
  bwCV <- function(xin, qin, optns) {
    compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)
    # leave-one-out cv
    objFctn <- function(bw) {
      optns1 <- optns
      optns1$bw <- bw
      qfit <- sapply(seq_along(xin)[compareRange], function(i) {
        res <- LocWassReg(xin = xin[-i], qin = qin[-i,], xout = xin[i], optns = optns1)
        as.vector(res)
      })
      mean(apply((qfit - t(qin[compareRange,]))^2, 2, pracma::trapz, x = optns1$qSup))
    }
    bwRange <- SetBwRange(xin = xin, ker = optns$ker)
    res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
    res$minimum
  }
}

# bandwidth selection via GCV implemented for the linear combination of quantile functions before optimization
# GCV = n^{-1}\sum_i \|g(x_i) - q(x_i)\|^2_{L^2} / (1 - n^{-1}\sum_i s_i(x_i))^2
# S = [s_j(x_i)]
# g(x) = \sum_i s_i(x) q(x_i)
bwGCV <- function(xin, qin, xout, optns) {
  compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    S <- LocWt(xin = xin[compareRange], qin = qin[compareRange,], xout = xin[compareRange], optns = optns1)
    g <- S %*% qin[compareRange,]
    mean(apply((g - qin[compareRange,])^2, 1, pracma::trapz, x = optns1$qSup)) /
      (1 - mean(diag(S)))^2
  }
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
  if (!is.null(optns$bwRange)) {
    if (min(optns$bwRange) < bwRange$min) {
      message("Minimum bandwidth is too small and has been reset.")
    } else bwRange$min <- min(optns$bwRange)
    if (max(optns$bwRange) >  bwRange$min) {
      bwRange$max <- max(optns$bwRange)
    } else {
      message("Maximum bandwidth is too small and has been reset.")
    }
  }
  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  res$minimum
}

# Weight used in local Wasserstein regression
LocWt = function(xin, qin, xout, optns = list()){

  if (is.null(optns$bw))
    stop ("optns$bw has no default values and must be input by users.")
  bw <- optns$bw
  if(!is.numeric(optns$bw) | (length(optns$bw)>1))
    stop("optns$bw should be a numerical vector of length 1.")
  if(!is.vector(xin) | !is.numeric(xin))
    stop("xin should be a numerical vector.")
  if(!is.matrix(qin) | !is.numeric(qin))
    stop("qin should be a numerical matrix.")
  if(length(xin)!=nrow(qin))
    stop("The length of xin should be the same as the number of rows of qin.")

  if (is.null(optns$ker)) {
    optns$ker <- 'gauss'
  }
  ker <- kerFctn(optns$ker)

  n = length(xout)
  #m = ncol(qin)

  S = sapply(1:n, function(j){
    mu0 = mean(ker((xout[j] - xin) / bw))
    mu1 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
    mu2 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
    s = ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
      (mu0 * mu2 - mu1^2)
    #gx = colMeans(qin * s)
    return(s)
  })
  S = t(S)/n
  return(S)
}
