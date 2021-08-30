#' @title Fréchet ANOVA for Densities
#' @description Fréchet analysis of variance for densities with respect to 
#'   \eqn{L^2}-Wasserstein distance.
#' @param Ly a list of \eqn{n} vectors containing the observed values for each 
#'   density/quantile function. See `Details`.
#' @param Lx a list of \eqn{n} vectors containing the support grid for each 
#'   density/quantile function corresponding to \code{Ly} or a vector 
#'   if all density/quantile functions have the same support grid.
#' @param group a vector containing the group memberships of the corresponding 
#'   densities/quantiles in \code{Ly}.
#' @param optns a list of control parameters specified by 
#'   \code{list(name = value)}. See `Details`.
#' @details Note that the type of functions representing the distributions 
#'   in \code{Ly} should be the same --- either all are density functions, 
#'   or all are quantile functions.
#'
#' If all support grids in \code{Lx} are the same, \code{qSup} will be set as 
#'   \code{Lx[[1]]}, no matter the value of \code{nqSup}.
#'
#' Available control options are:
#' \describe{
#' \item{fctn_type}{function type in \code{Ly} representing the distributions: 
#'   \code{"density"} (default), \code{"quantile"}.}
#' \item{boot}{logical, also compute bootstrap \eqn{p}-value if \code{TRUE}. 
#'   Default is \code{FALSE}.}
#' \item{R}{the number of bootstrap replicates. Only used when \code{boot} 
#'   is \code{TRUE}. Default is 1000.}
#' \item{nqSup}{a scalar giving the length of the support grid of 
#'   quantile functions based on which the \eqn{L^2} Wasserstein distance 
#'   (i.e., the \eqn{L^2} distance between the quantile functions) is computed. 
#'   Default is 201.}
#' }
#' @return A \code{DenANOVA} object --- a list containing the following fields:
#' \item{pvalAsy}{a scalar holding the asymptotic \eqn{p}-value.}
#' \item{pvalBoot}{a scalar holding the bootstrap \eqn{p}-value.}
#' \item{qSup}{a numeric vector holding the support grid used.}
#' \item{optns}{the control options used.}
#' @examples
#' set.seed(1)
#' n1 <- 100
#' n2 <- 100
#' delta <- 1
#' qSup <- seq(0.01, 0.99, (0.99 - 0.01) / 50)
#' mu1 <- rnorm(n1, mean = 0, sd = 0.5)
#' mu2 <- rnorm(n2, mean = delta, sd = 0.5)Ly
#' Y1 <- lapply(1:n1, function(i) {
#'   qnorm(qSup, mu1[i], sd = 1)
#' })
#' Y2 <- lapply(1:n2, function(i) {
#'   qnorm(qSup, mu2[i], sd = 1)
#' })
#' Ly <- c(Y1, Y2)
#' Lx <- qSup
#' group <- c(rep(1, n1), rep(2, n2))
#' res <- DenANOVA(Ly, Lx, group, optns = list(fctn_type = "quantile", boot = TRUE))
#' res$pvalAsy # returns asymptotic pvalue
#' res$pvalBoot # returns bootstrap pvalue
#' @references
#' \itemize{
#' \item \cite{Dubey, P. and Müller, H.G., 2019. Fréchet analysis of variance for random objects. Biometrika, 106(4), pp.803-821.}
#' }
#' @export

DenANOVA <- function(Ly = NULL, Lx = NULL, group = NULL, optns = list()) {
  if (is.null(Ly) | is.null(Lx) | is.null(group)) {
    stop("requires the input of Ly, Lx, and group")
  }
  if (!is.list(Ly)) {
    stop("Ly should be a list")
  }
  if (!is.list(Lx)) {
    Lx <- rep(list(Lx), length(Ly))
  } else{
    if (length(Ly) != length(Lx)) {
      stop("Ly and Lx should have the same length")
    }
  }
  if (length(Ly) != length(group)) {
    stop("Ly and group should have the same length")
  }
  if (sum(sapply(Ly, length) - sapply(Lx, length))) {
    stop("each vector in Ly and its corresponding vector in Lx should have the same length")
  }
  if (is.null(optns$fctn_type)) {
    fctn_type <- "density"
  } else {
    fctn_type <- optns$fctn_type
  }
  if (length(fctn_type) > 1) {
    fctn_type <- fctn_type[1]
    warning("fctn_type has length greater than 1---only the first element is used")
  }
  if (!fctn_type %in% c("density", "quantile")) {
    stop("unrecognized value of fctn_type")
  }
  if (is.null(optns$nqSup)) {
    nqSup <- 201
  } else {
    nqSup <- optns$nqSup
  }
  if (is.null(optns$boot)) {
    boot <- FALSE
  } else {
    boot <- optns$boot
  }
  if (is.null(optns$R)) {
    R <- 1000
  } else {
    R <- optns$R
  }
  n <- length(Ly)
  # tol <- 1e-5
  if (fctn_type == "density") {
    # if (any(sapply(1:n, function(i) {
    #   abs(pracma::trapz(Lx[[i]], Ly[[i]]) - 1) > tol
    # }))) {
    #   stop("each element of Ly should be a density function (integrates to $1$ with tolerance of ", tol, ")")
    # }
    # if (any(sapply(1:n, function(i) {
    #   any(Ly[[i]] < 0)
    # }))) {
    #   stop("each vector in Ly should be all non-negative")
    # }
    # fdadensity::dens2quantile() will detect if each Ly[[i]] is indeed a density function
    qSup <- seq(0, 1, length.out = nqSup)
    Ly <- lapply(1:n, function(i) fdadensity::dens2quantile(Ly[[i]], dSup = Lx[[i]], qSup = qSup))
  } else {
    if (any(sapply(Lx, function(Lxi) {
      any(Lxi < 0 | Lxi > 1)
    }))) {
      stop("each vector in Lx should lie in [0, 1]")
    }
    for (i in 1:n) {
      if (is.unsorted(Lx[[i]])) {
        Ly[[i]] <- Ly[[i]][order(Lx[[i]])]
        Lx[[i]] <- sort(Lx[[i]])
      }
    }
    if (any(sapply(1:n, function(i) is.unsorted(Ly[[i]])))) {
      stop("quantile functions given in Ly are not monotonic")
    }
    diffSupp <- TRUE
    if (length(unique(sapply(Lx, length))) == 1) {
      if (sum(diff(matrix(unlist(Lx), nrow = n, byrow = TRUE))) == 0) diffSupp <- FALSE
    }
    if (diffSupp) {
      qSup <- seq(0, 1, length.out = nqSup)
      Ly <- lapply(1:n, function(i) approx(x = Lx[[i]], y = Ly[[i]], xout = qSup)$y)
    } else {
      qSup <- Lx[[1]]
    }
  }
  sizes <- as.vector(table(group))
  k <- length(sizes) # number of groups
  data <- Ly[order(group)]
  # data <- unlist(split(Ly, group), recursive = FALSE)
  if (boot) {
    bootRes <- boot::boot(data = data, statistic = DenANOVAStatistic, 
                          R = R, sizes = sizes, qSup = qSup)
    pvalBoot <- length(which(bootRes$t > bootRes$t0)) / R
    pvalAsy <- 1 - pchisq(bootRes$t0, df = k - 1)
    res <- list(pvalAsy = pvalAsy, pvalBoot = pvalBoot, qSup = qSup, optns = optns)
  } else {
    t0 <- DenANOVAStatistic(data, 1:length(group), sizes, qSup)
    pvalAsy <- 1 - pchisq(t0, df = k - 1)
    res <- list(pvalAsy = pvalAsy, qSup = qSup, optns = optns)
  }
  class(res) <- "DenANOVA"
  res
}
