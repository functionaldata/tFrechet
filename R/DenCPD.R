#' @title Fréchet Change Point Detection for Densities
#' @description Fréchet change point detection for densities with respect to
#'   \eqn{L^2}-Wasserstein distance.
#' @param yin A matrix or data frame or list holding the sample of measurements
#'   for the observed distributions. If \code{yin} is a matrix or data frame,
#'   each row holds the measurements for one distribution.
#' @param hin A list holding the histograms for the observed distributions.
#' @param din A matrix or data frame or list holding the density functions.
#'   If \code{din} is a matrix or data frame, each row of \code{din} holds
#'   the density function for one distribution.
#' @param qin A matrix or data frame or list holding the quantile functions.
#'   If \code{qin} is a matrix or data frame, each row of \code{qin} holds
#'   the quantile function for one distribution.
#' Note that the input can be only one of the four \code{yin}, \code{hin},
#' \code{din}, and \code{qin}. If more than one of them are specified,
#' \code{yin} overwrites \code{hin}, \code{hin} overwrites \code{din},
#' and \code{din} overwrites \code{qin}.
#' @param supin A matrix or data frame or list holding the support grids of
#'   the density functions in \code{din} or the quantile functions in \code{qin}.
#'   If \code{supin} is a matrix or data frame, each row of \code{supin} holds
#'   the support grid of the corresponding density function or quantile function.
#'   Ignored if the input is \code{yin} or \code{hin}.
#'   It can also be a vector if all density functions in \code{din} or
#'   all quantile functions in \code{qin} have the same support grid.
#' @param optns A list of control parameters specified by
#'   \code{list(name = value)}. See `Details`.
#' @details Available control options are
#' \describe{
#' \item{cutOff}{A scalar between 0 and 1 indicating the interval,
#'   i.e., [cutOff, 1 - cutOff], in which candidate change points lie.}
#' \item{Q}{A scalar representing the number of Monte Carlo simulations to run
#'   while approximating the critical value (stardized Brownian bridge).
#'   Default is 1000.}
#' \item{boot}{Logical, also compute bootstrap \eqn{p}-value if \code{TRUE}.
#'   Default is \code{FALSE}.}
#' \item{R}{The number of bootstrap replicates. Only used when \code{boot}
#'   is \code{TRUE}. Default is 1000.}
#' \item{nqSup}{A scalar giving the number of the support points for
#'   quantile functions based on which the \eqn{L^2} Wasserstein distance
#'   (i.e., the \eqn{L^2} distance between the quantile functions) is computed.
#'   Default is 201.}
#' \item{qSup}{A numeric vector holding the support grid on [0, 1] based on
#'   which the \eqn{L^2} Wasserstein distance (i.e., the \eqn{L^2} distance
#'   between the quantile functions) is computed. It overrides \code{nqSup}.}
#' \item{bwDen}{The bandwidth value used in \code{CreateDensity()} for
#'   density estimation; positive numeric - default: determine automatically
#'   based on the data-driven bandwidth selector proposed by
#'   Sheather and Jones (1991).}
#' \item{ndSup}{A scalar giving the number of support points the kernel density
#'   estimation used in \code{CreateDensity()}; numeric - default: 101.}
#' \item{dSup}{User defined output grid for the support of
#'   kernel density estimation used in \code{CreateDensity()},
#'   it overrides \code{ndSup}.}
#' \item{delta}{A scalar giving the size of the bin to be used used in
#'   \code{CreateDensity()}; numeric - default: \code{diff(range(y))/1000}.
#'   It only works when the raw sample is available.}
#' \item{kernelDen}{A character holding the type of kernel functions used in
#'   \code{CreateDensity()} for density estimation; \code{"rect"},
#'   \code{"gauss"}, \code{"epan"}, \code{"gausvar"},
#'   \code{"quar"} - default: \code{"gauss"}.}
#' \item{infSupport}{logical if we expect the distribution to have
#'   infinite support or not, used in \code{CreateDensity()} for
#'   density estimation; logical - default: \code{FALSE}}
#' \item{denLowerThreshold}{\code{FALSE} or a positive value giving
#'   the lower threshold of the densities used in \code{CreateDensity()};
#'   default: \code{0.001 * mean(qin[,ncol(qin)] - qin[,1])}.}
#' }
#' @return A \code{DenCPD} object --- a list containing the following fields:
#' \item{tau}{a scalar holding the estimated change point.}
#' \item{pvalAsy}{a scalar holding the asymptotic \eqn{p}-value.}
#' \item{pvalBoot}{a scalar holding the bootstrap \eqn{p}-value.
#'   Returned if \code{optns$boot} is TRUE.}
#' \item{optns}{the control options used.}
#' @examples
#' \donttest{
#' set.seed(1)
#' n1 <- 100
#' n2 <- 200
#' delta <- 0.75
#' qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
#' mu1 <- rnorm(n1, mean = delta, sd = 0.5)
#' mu2 <- rnorm(n2, mean = 0, sd = 0.5)
#' Y1 <- lapply(1:n1, function(i) {
#'   qnorm(qSup, mu1[i], sd = 1)
#' })
#' Y2 <- lapply(1:n2, function(i) {
#'   qnorm(qSup, mu2[i], sd = 1)
#' })
#' Ly <- c(Y1, Y2)
#' Lx <- qSup
#' res <- DenCPD(qin = Ly, supin = Lx, optns = list(boot = TRUE))
#' res$tau # returns the estimated change point
#' res$pvalAsy # returns asymptotic pvalue
#' res$pvalBoot # returns bootstrap pvalue
#' } 
#' @references
#' \itemize{
#' \item \cite{Dubey, P. and Müller, H.G., 2020. Fréchet change-point detection. The Annals of Statistics, 48(6), pp.3312-3335.}
#' }
#' @importFrom e1071 rbridge
#' @importFrom fdadensity dens2quantile
#' @importFrom boot boot
#' @export

DenCPD <- function(yin = NULL, hin = NULL, din = NULL, qin = NULL,
                   supin = NULL, optns = list()) {
  if (is.null(yin) & is.null(qin) & is.null(hin) & is.null(din)) {
    stop("one of the four arguments, yin, hin, din, and qin, should be inputted by users")
  }
  if (!is.null(yin)) {
    if (!(is.null(hin) & is.null(din) & is.null(qin))) {
      warning("hin, din, and qin are redundant when yin is available")
    }
    tin <- 1 # type of input
    ls <- yin
  } else if (!is.null(hin)) {
    if (!(is.null(din) & is.null(qin))) {
      warning("din and qin are redundant when hin is available")
    }
    tin <- 2 # type of input
    ls <- hin
  } else if (!is.null(din)) {
    if (!is.null(qin)) {
      warning("qin is redundant when din is available")
    }
    tin <- 3 # type of input
    ls <- din
    if (!is.list(din)) {
      if (is.matrix(din) | is.data.frame(din)) {
        din <- lapply(1:nrow(din), function(i) din[i, ])
      } else {
        stop("din must be a matrix or data frame or list")
      }
    }
  } else {
    tin <- 4 # type of input
    ls <- qin
  }
  if (tin == 2) {
    if (!is.list(ls)) {
      stop("hin must be a list")
    }
    for (histogram in ls) {
      if (!is.list(histogram)) {
        stop("each element of hin must be a list")
      }
      if (is.null(histogram$breaks) & is.null(histogram$mids)) {
        stop("each element of hin must be a list with at least one of the components breaks or mids")
      }
      if (is.null(histogram$counts)) {
        stop("each element of hin must be a list with component counts")
      }
    }
  } else {
    if (!is.list(ls)) {
      if (is.matrix(ls) | is.data.frame(ls)) {
        ls <- lapply(1:nrow(ls), function(i) ls[i, ])
      } else {
        if (tin == 1) {
          stop("yin must be a matrix or data frame or list")
        } else if (tin == 3) {
          stop("din must be a matrix or data frame or list")
        } else {
          stop("qin must be a matrix or data frame or list")
        }
      }
    }
  }
  n <- length(ls)
  if (is.null(optns$cutOff)) {
    optns$cutOff <- 0.1
  }
  if (is.null(optns$Q)) {
    optns$Q <- 1000
  }
  if (is.null(optns$boot)) {
    optns$boot <- FALSE
  }
  if (optns$boot) {
    if (is.null(optns$R)) {
      optns$R <- 1000
    }
  }
  if (tin > 2) { # if din or qin
    if (!is.null(supin)) {
      if (!is.list(supin)) {
        if (is.matrix(supin) | is.data.frame(supin)) {
          supin <- lapply(1:nrow(supin), function(i) supin[i, ])
        } else if (is.vector(supin)) {
          supin <- rep(list(supin), n)
        } else {
          stop("supin must be a vector or matrix or data frame or list")
        }
      }
      if (length(supin) != n) {
        if (tin == 3) {
          stop("the number of support grids in supin is not equal to the number of observed distributions in din")
        } else {
          stop("the number of support grids in supin is not equal to the number of observed distributions in qin")
        }
      }
      if (sum(sapply(supin, length) - sapply(ls, length))) {
        if (tin == 3) {
          stop("the number of support points must be equal to the number of observations for each density function in din")
        } else {
          stop("the number of support points must be equal to the number of observations for each quantile function in qin")
        }
      }
    } else {
      if (tin == 3) {
        stop("requires the input of supin for din")
      } else {
        stop("requires the input of supin for qin")
      }
    }
  }
  if (!is.null(optns$qSup)) {
    if (min(optns$qSup) != 0 | max(optns$qSup) - 1 != 0) {
      stop("optns$qSup must have minimum 0 and maximum 1")
    }
    if (sum(duplicated(optns$qSup)) > 0) {
      optns$qSup <- unique(optns$qSup)
      warning("optns$qSup has duplicated elements which has been removed")
    }
    if (is.unsorted(optns$qSup)) {
      optns$qSup <- sort(optns$qSup)
      warning("optns$qSup has been reordered to be increasing")
    }
  } else {
    if (tin < 4) { # if not qin
      if (is.null(optns$nqSup)) {
        optns$nqSup <- 201
      }
      optns$qSup <- seq(0, 1, length.out = optns$nqSup)
    } else { # if qin
      diffSupp <- TRUE
      if (length(unique(sapply(supin, length))) == 1) {
        if (!sum(diff(matrix(unlist(supin), nrow = n, byrow = TRUE)))) {
          diffSupp <- FALSE
        }
      }
      if (diffSupp) {
        if (is.null(optns$nqSup)) {
          optns$nqSup <- 201
        }
        optns$qSup <- seq(0, 1, length.out = optns$nqSup)
        ls <- lapply(1:n, function(i) approx(x = supin[[i]], y = ls[[i]], xout = optns$qSup)$y)
      } else {
        optns$qSup <- supin[[1]]
        optns$nqSup <- length(optns$qSup)
      }
    }
  }
  qSup <- optns$qSup
  if (tin == 3) {
    ls <- lapply(1:n, function(i) fdadensity::dens2quantile(ls[[i]], dSup = supin[[i]], qSup = qSup))
  }
  if (tin < 3) {
    optnsDen <- optns
    if (!is.null(optnsDen$kernelDen)) {
      names(optnsDen)[which(names(optnsDen) == "kernelDen")] <- "kernel"
    }
    if (!is.null(optnsDen$bwDen)) {
      names(optnsDen)[which(names(optnsDen) == "bwDen")] <- "userBwMu"
    }
    if (!is.null(optnsDen$ndSup)) {
      names(optnsDen)[which(names(optnsDen) == "ndSup")] <- "nRegGrid"
    }
    if (!is.null(optnsDen$dSup)) {
      names(optnsDen)[which(names(optnsDen) == "dSup")] <- "outputGrid"
    }
    if (tin == 1) {
      den <- lapply(ls, CreateDensity, optns = optnsDen)
    } else {
      den <- lapply(ls, function(histogram) {
        CreateDensity(histogram = histogram, optns = optnsDen)
      })
    }
    ls <- lapply(den, function(deni) fdadensity::dens2quantile(deni$y, dSup = deni$x, qSup = qSup))
  }

  nc <- ceiling(n * optns$cutOff)
  sbb <- sapply(1:optns$Q, function(i) { # standardized Brownian bridge
    bu <- e1071::rbridge(frequency = n)[nc:(n - nc)]
    u <- (nc:(n - nc)) / n
    max(bu^2 / (u * (1 - u)))
  })
  tTau <- DenCPDStatistic(ls, 1:n, optns$cutOff, qSup, n)
  maxnTn <- tTau[1] # the test statistic
  tau <- tTau[2] + nc - 1 # estimated change point
  pvalAsy <- length(which(sbb > maxnTn)) / optns$Q
  if (optns$boot) {
    bootSize <- n # check bootstrap scheme in the paper
    bootRes <- boot::boot(
      data = ls, statistic = DenCPDStatistic,
      R = optns$R, cutOff = optns$cutOff, qSup = qSup, bootSize = bootSize)
    pvalBoot <- length(which(bootRes$t[, 1] > bootRes$t0[1])) / optns$R
    res <- list(tau = tau, pvalAsy = pvalAsy, pvalBoot = pvalBoot, optns = optns)
  } else {
    res <- list(tau = tau, pvalAsy = pvalAsy, optns = optns)
  }
  class(res) <- "DenCPD"
  res
}
