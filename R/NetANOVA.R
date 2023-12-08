#' @title Fréchet ANOVA for Networks
#' @description Fréchet analysis of variance for graph Laplacian matrices, 
#'   covariance matrices, or correlation matrices 
#'   with respect to the Frobenius distance.
#' @param Ly A list (length n) of m by m matrices or a m by m by n array where
#'   \code{Ly[, , i]} contains an m by m matrix, which can be either graph 
#'   Laplacian matrices or covariance matrices or correlation matrices.
#' @param group A vector containing the group memberships of the corresponding 
#'   matrices in \code{Ly}.
#' @param optns A list of control parameters specified by 
#'   \code{list(name = value)}. See `Details`.
#' @details Available control options are:
#' \describe{
#' \item{boot}{Logical, also compute bootstrap \eqn{p}-value if \code{TRUE}. 
#'   Default is \code{FALSE}.}
#' \item{R}{The number of bootstrap replicates. Only used when \code{boot} 
#'   is \code{TRUE}. Default is 1000.}
#' }
#' @return A \code{NetANOVA} object --- a list containing the following fields:
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
#' group <- c(rep(1, n1), rep(2, n2))
#' res <- NetANOVA(Ly, group, optns = list(boot = TRUE))
#' res$pvalAsy # returns asymptotic pvalue
#' res$pvalBoot # returns bootstrap pvalue
#' }
#' @references
#' \itemize{
#' \item \cite{Dubey, P. and Müller, H.G., 2019. Fréchet analysis of variance for random objects. Biometrika, 106(4), pp.803-821.}
#' }
#' @importFrom boot boot
#' @importFrom stats pchisq
#' @export

NetANOVA <- function(Ly = NULL, group = NULL, optns = list()) {
  if (is.null(Ly) | is.null(group)) {
    stop("requires the input of both Ly and group")
  }
  if (!is.list(Ly)) {
    if (is.array(Ly)) {
      Ly <- lapply(seq(dim(Ly)[3]), function(i) Ly[, , i])
    } else {
      stop("Ly must be a list or an array")
    }
  }
  if (length(Ly) != length(group)) {
    stop("Ly and group should have the same length")
  }
  if (length(unique(sapply(Ly, length))) > 1) {
    stop("each matrix in Ly should be of the same dimension")
  }
  if (any(sapply(Ly, function(Lyi) nrow(Lyi) != ncol(Lyi)))) {
    stop("each matrix in Ly should be a square matrix")
  }
  if (is.null(optns$boot)) {
    optns$boot <- FALSE
  }
  if (is.null(optns$R)) {
    optns$R <- 1000
  }
  n <- length(Ly)
  sizes <- as.vector(table(group))
  k <- length(sizes) # number of groups
  data <- Ly[order(group)]
  # data <- unlist(split(Ly, group), recursive = FALSE)
  if (optns$boot) {
    bootRes <- boot::boot(data = data, statistic = NetANOVAStatistic, 
                          R = optns$R, sizes = sizes)
    pvalBoot <- length(which(bootRes$t > bootRes$t0)) / optns$R
    pvalAsy <- 1 - pchisq(bootRes$t0, df = k - 1)
    res <- list(pvalAsy = pvalAsy, pvalBoot = pvalBoot, optns = optns)
  } else {
    t0 <- NetANOVAStatistic(data, 1:length(group), sizes)
    pvalAsy <- 1 - pchisq(t0, df = k - 1)
    res <- list(pvalAsy = pvalAsy, optns = optns)
  }
  class(res) <- "NetANOVA"
  res
}
