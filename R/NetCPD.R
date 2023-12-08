#' @title Fréchet Change Point Detection for Networks
#' @description Fréchet change point detection for graph Laplacian matrices, 
#'   covariance matrices, or correlation matrices 
#'   with respect to the Frobenius distance.
#' @param Ly A list (length n) of m by m matrices or a m by m by n array where
#'   \code{Ly[, , i]} contains an m by m matrix, which can be either graph 
#'   Laplacian matrices or covariance matrices or correlation matrices.
#' @param optns A list of control parameters specified by 
#'   \code{list(name = value)}. See `Details`.
#' @details Available control options are:
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
#' }
#' @return A \code{NetCPD} object --- a list containing the following fields:
#' \item{tau}{a scalar holding the estimated change point.}
#' \item{pvalAsy}{A scalar holding the asymptotic \eqn{p}-value.}
#' \item{pvalBoot}{A scalar holding the bootstrap \eqn{p}-value. 
#'   Returned if \code{optns$boot} is TRUE.}
#' \item{optns}{The control options used.}
#' @examples
#' \donttest{
#' set.seed(1)
#' n1 <- 100
#' n2 <- 100
#' gamma1 <- 2
#' gamma2 <- 3
#' Y1 <- lapply(1:n1, function(i) {
#'   igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma1, 
#'                                              directed = FALSE), 
#'                            sparse = FALSE)
#' })
#' Y2 <- lapply(1:n2, function(i) {
#'   igraph::laplacian_matrix(igraph::sample_pa(n = 10, power = gamma2, 
#'                                              directed = FALSE), 
#'                            sparse = FALSE)
#' })
#' Ly <- c(Y1, Y2)
#' res <- NetCPD(Ly, optns = list(boot = TRUE))
#' res$tau # returns the estimated change point
#' res$pvalAsy # returns asymptotic pvalue
#' res$pvalBoot # returns bootstrap pvalue
#' }
#' @references
#' \itemize{
#' \item \cite{Dubey, P. and Müller, H.G., 2020. Fréchet change-point detection. The Annals of Statistics, 48(6), pp.3312-3335.}
#' }
#' 
#' @importFrom e1071 rbridge
#' @importFrom boot boot
#' @export

NetCPD <- function(Ly = NULL, optns = list()) {
  if (is.null(Ly)) {
    stop("requires the input of Ly")
  }
  if (!is.list(Ly)) {
    if (is.array(Ly)) {
      Ly <- lapply(seq(dim(Ly)[3]), function(i) Ly[, , i])
    } else {
      stop("Ly must be a list or an array")
    }
  }
  if (is.null(optns$cutOff)) {
    optns$cutOff <- 0.1
  }
  if (is.null(optns$Q)) {
    optns$Q <- 1000
  }
  if (is.null(optns$boot)) {
    optns$boot <- FALSE
  }
  if (is.null(optns$R)) {
    optns$R <- 1000
  }
  n <- length(Ly)
  nc <- ceiling(n * optns$cutOff)
  sbb <- sapply(1:optns$Q, function(i) { # standardized Brownian bridge
    bu <- e1071::rbridge(frequency = n)[nc:(n - nc)]
    u <- (nc:(n - nc)) / n
    max(bu^2 / (u * (1 - u)))
  })
  tTau <- NetCPDStatistic(Ly, 1:n, optns$cutOff, n)
  maxnTn <- tTau[1] # the test statistic
  tau <- tTau[2] + nc - 1 # estimated change point
  pvalAsy <- length(which(sbb > maxnTn)) / optns$Q
  if (optns$boot) {
    bootSize <- n # check bootstrap scheme in the paper
    bootRes <- boot::boot(
      data = Ly, statistic = NetCPDStatistic,
      R = optns$R, cutOff = optns$cutOff, bootSize = bootSize)
    pvalBoot <- length(which(bootRes$t[, 1] > bootRes$t0[1])) / optns$R
    res <- list(tau = tau, pvalAsy = pvalAsy, pvalBoot = pvalBoot, optns = optns)
  } else {
    res <- list(tau = tau, pvalAsy = pvalAsy, optns = optns)
  }
  class(res) <- "NetCPD"
  res
}
