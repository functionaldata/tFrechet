#' @title Create density functions from raw data, histogram objects or frequency tables with bins
#'
#' @description Create kernel density estimate along the support of the raw data using the HADES method.
#'
#' @param y A vector of raw readings.
#' @param histogram A \code{histogram} object in R. Use this option when histogram object is only available, but not the raw data \code{y}. The default is \code{NULL}.
#' @param freq A frequency vector. Use this option when frequency table is only available, but not the raw sample or the histogram object. The corresponding \code{bin} should be provided together. The default is \code{NULL}.
#' @param bin A bin vector having its length with \code{length(freq)+1}. Use this option when frequency table is only available, but not the raw sample or the histogram object. The corresponding \code{freq} should be provided together.The default is \code{NULL}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{userBwMu}{The bandwidth value for the smoothed mean function; positive numeric - default: determine automatically based on the data-driven bandwidth selector proposed by Sheather and Jones (1991)}
#' \item{nRegGrid}{The number of support points the KDE; numeric - default: 101.}
#' \item{delta}{The size of the bin to be used; numeric - default: \code{diff(range(y))/1000}. It only works when the raw sample is available.}
#' \item{kernel}{smoothing kernel choice, \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' \item{infSupport}{logical if we expect the distribution to have infinite support or not; logical - default: \code{FALSE}.}
#' \item{outputGrid}{User defined output grid for the support of the KDE, it overrides \code{nRegGrid}; numeric - default: \code{NULL}.}
#' }
#'
#' @return A list containing the following fields:
#' \item{bw}{The bandwidth used for smoothing.}
#' \item{x}{A vector of length \code{nRegGrid} with the values of the KDE's support points.}
#' \item{y}{A vector of length \code{nRegGrid} with the values of the KDE at the support points.}
#'
#' @examples
#'
#' ### compact support case
#'
#' # input: raw sample
#' set.seed(100)
#' n <- 100
#' x0 <-seq(0,1,length.out=51)
#' Y <- rbeta(n,3,2)
#' f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
#'
#' # input: histogram
#' histY <- hist(Y)
#' f2 <- CreateDensity(histogram=histY,optns = list(outputGrid=x0))
#'
#' # input: frequency table with unequally spaced (random) bins
#' binY <- c(0,sort(runif(9)),1)
#' freqY <- c()
#' for (i in 1:(length(binY)-1)) {
#'   freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
#' }
#' f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
#'
#' # plot
#' plot(f1$x,f1$y,type='l',col=2,lty=2,lwd=2,
#'      xlim=c(0,1),ylim=c(0,2),xlab='domain',ylab='density')
#' points(f2$x,f2$y,type='l',col=3,lty=3,lwd=2)
#' points(f3$x,f3$y,type='l',col=4,lty=4,lwd=2)
#' points(x0,dbeta(x0,3,2),type='l',lwd=2)
#' legend('topleft',
#'        c('true','raw sample','histogram','frequency table (unequal bin)'),
#'        col=1:4,lty=1:4,lwd=3,bty='n')
#'
#' ### infinite support case
#'
#' # input: raw sample
#' set.seed(100)
#' n <- 200
#' x0 <-seq(-3,3,length.out=101)
#' Y <- rnorm(n)
#' f1 <- CreateDensity(y=Y,optns = list(outputGrid=x0))
#'
#' # input: histogram
#' histY <- hist(Y)
#' f2 <- CreateDensity(histogram=histY,optns = list(outputGrid=x0))
#'
#' # input: frequency table with unequally spaced (random) bins
#' binY <- c(-3,sort(runif(9,-3,3)),3)
#' freqY <- c()
#' for (i in 1:(length(binY)-1)) {
#'   freqY[i] <- length(which(Y>binY[i] & Y<=binY[i+1]))
#' }
#' f3 <- CreateDensity(freq=freqY, bin=binY,optns = list(outputGrid=x0))
#'
#' # plot
#' plot(f1$x,f1$y,type='l',col=2,lty=2,lwd=2,
#'      xlim=c(-3,3),ylim=c(0,0.5),xlab='domain',ylab='density')
#' points(f2$x,f2$y,type='l',col=3,lty=3,lwd=2)
#' points(f3$x,f3$y,type='l',col=4,lty=4,lwd=2)
#' points(x0,dnorm(x0),type='l',lwd=2)
#' legend('topright',
#'        c('true','raw sample','histogram','frequency table (unequal bin)'),
#'        col=1:4,lty=1:4,lwd=3,bty='n')
#'
#' @references
#' \itemize{
#' \item \cite{H.-G. Müller, J.L. Wang and W.B. Capra (1997). "From lifetables to hazard rates: The transformation approach." Biometrika 84, 881--892.}
#' \item \cite{S.J. Sheather and M.C. Jones (1991). "A reliable data-based bandwidth selection method for kernel density estimation." JRSS-B 53, 683--690.}
#' \item \cite{H.-G. Müller, U. Stadtmüller, and T. Schmitt. (1987) "Bandwidth choice and confidence intervals for derivatives of noisy data." Biometrika 74, 743--749.}
#' }
#' @export

CreateDensity <- function(y=NULL, histogram=NULL, freq=NULL, bin=NULL, optns = list()){

  if (is.null(y)==TRUE) {

    if (is.null(histogram)==FALSE) {
      if (is.null(histogram$breaks)) {
        mids <- histogram$mids
        bin <- mids[1]-(mids[2]-mids[1])/2
        for (i in 2:length(mids)) {
          bin[i] <- (mids[i-1]+mids[i])/2
        }
        bin[length(mids)+1] <- mids[length(mids)]+(mids[length(mids)]-mids[length(mids)-1])/2
      } else bin <- histogram$breaks

      freq <- histogram$counts
    }

    if (is.null(freq)==TRUE) {

      if ((length(freq)+1)!=length(bin)) {
        stop('length(bin) should equal to length(freq)+1.')
      }

      # mids <- c()
      # for (i in 1:length(freq)) {
      #   mids[i] <- (bin[i]+bin[i+1])/2
      # }
    }

    y <- c()
    for (i in 1:length(freq)) {
      if (freq[i]!=0) {
        yTmp <- seq(bin[i],bin[i+1],length.out=(freq[i]+2))
        y <- c(y,yTmp[-c(1,length(yTmp))])
      }
    }
  }

  if(is.null(optns$kernel)){
    kernel = 'gauss'
  } else {
    kernel =  optns$kernel
  }

  if(is.null(optns$delta)){
    delta = diff(range(y))/1000
    #delta = max(c( diff(range(y))/1000, min(diff(sort(unique(y)))) ))
  } else {
    delta = optns$delta
  }

  if(is.null(optns$nRegGrid)){
    nRegGrid = 101
  } else {
    nRegGrid = optns$nRegGrid
  }

  if(is.null(optns$outputGrid)){
    outputGrid = NULL
  } else {
    outputGrid = optns$outputGrid
  }

  if(is.null(optns$infSupport)){
    infSupport = FALSE
  } else {
    infSupport = optns$infSupport
  }

  N = length(y)
  #histgrid = seq( min(y)-delta*0.5, max(y)+delta*0.5, by = delta );
  histgrid = c(seq( min(y), max(y), by = delta) - delta*0.5, max(y)+delta*0.5)
  if(  (max(y)+delta*0.5) > histgrid[length(histgrid)] ){
    histgrid[length(histgrid)] =  max(y)+delta*0.5;
  }
  M = length(histgrid)
  histObj =  hist(y, breaks = histgrid, plot = FALSE);
  yin = histObj$density; xin = histObj$mids
  #yin = histObj$counts[1:M-1] / N / delta;
  #xin = seq( min(y), max(y), by = delta);

  if( is.null(optns$userBwMu)){

    densTmp <- density(y,kernel='epanechnikov')
    bw <- 2*densTmp$bw
    #bw <- 2*densTmp$bw

    # bwCV = fdapace:::CVLwls1D(y = yin, t = xin, kernel = kernel, npoly = 1, nder = 0, dataType = 'Dense', kFolds = 5)
    # bw <- bwCV
  } else {
    bw = optns$userBwMu
  }

  densObj <- list()
  densObj$bw <- bw
  densObj$x <- outputGrid
  if( infSupport ){
    if(is.null(outputGrid)){
      densObj$x <- seq(min(y)-delta, max(y)+delta, length.out = nRegGrid)
    }
    qpadding = 10
    mu = fdapace::Lwls1D(bw = bw, kernel_type = kernel, win = rep(1,M+(-1+2*qpadding)),
                         xin = c( min(y) - 11:2 *delta, xin, max(y) + 2:11 * delta), yin = c(rep(0,qpadding),yin,rep(0,qpadding)), xout = densObj$x)
  } else {
    if(is.null(outputGrid)){
      densObj$x <- outputGrid <- seq(min(y), max(y), length.out = nRegGrid)
    }
    if (min(xin) > outputGrid[1]) {
      qpad1 <- seq(outputGrid[1],min(xin),by=delta)
    } else qpad1 <- numeric()
    if (max(xin) < outputGrid[length(outputGrid)]) {
      qpad2 <- seq(outputGrid[length(outputGrid)],max(xin),by=-delta)
      qpad2 <- qpad2[length(qpad2):1]
    } else qpad2 <- numeric()
    xinNew <- c(qpad1,xin,qpad2)
    yinNew <- c(rep(0,length(qpad1)),yin,rep(0,length(qpad2)))

    # note that Lwls1D requires xin to be at least of length 6
    mu = fdapace::Lwls1D(bw = bw, kernel_type = kernel, win = rep(1,length(xinNew)), xin = xinNew, yin = yinNew, xout = densObj$x)
  }
  mu[mu<0] = 0
  #densObj$y = mu / fdapace:::trapzRcpp(densObj$x, mu)
  densObj$y = mu / pracma::trapz(densObj$x, mu)

  return(densObj)
}






