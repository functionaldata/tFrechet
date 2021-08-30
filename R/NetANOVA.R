#' @title Fréchet ANOVA for Networks
#' @description Fréchet analysis of variance for graph Laplacian matrices, 
#'   covariance matrices, or correlation matrices 
#'   with respect to the Frobenius distance.
#' @param Ly a list (length \eqn{n}) of \eqn{m} by \eqn{m} matrices, 
#'   which can be either graph Laplacian matrices or covariance matrices 
#'   or correlation matrices. 
#' @param group a vector containing the group memberships of the corresponding 
#'   matrices in \code{Ly}.
#' @param optns a list of control parameters specified by 
#'   \code{list(name = value)}. See `Details`.
#' @details Available control options are:
#' \describe{
#' \item{boot}{logical, also compute bootstrap \eqn{p}-value if \code{TRUE}. 
#'   Default is \code{FALSE}.}
#' \item{R}{the number of bootstrap replicates. Only used when \code{boot} 
#'   is \code{TRUE}. Default is 1000.}
#' }
#' @return A \code{NetANOVA} object --- a list containing the following fields:
#' \item{pvalAsy}{a scalar holding the asymptotic \eqn{p}-value.}
#' \item{pvalBoot}{a scalar holding the bootstrap \eqn{p}-value.}
#' \item{optns}{the control options used.}
#' @examples
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
#' @references
#' \itemize{
#' \item \cite{Dubey, P. and Müller, H.G., 2019. Fréchet analysis of variance for random objects. Biometrika, 106(4), pp.803-821.}
#' }
#' @export

NetANOVA <- function(Ly = NULL, group = NULL, optns = list()) {
  if (is.null(Ly) | is.null(group)) {
    stop("requires the input of both Ly and group")
  }
  if (!is.list(Ly)) {
    stop("Ly should be a list")
  }
  if (length(Ly) != length(group)) {
    stop("Ly and group should have the same length")
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
  sizes <- as.vector(table(group))
  k <- length(sizes) # number of groups
  data <- Ly[order(group)]
  # data <- unlist(split(Ly, group), recursive = FALSE)
  if (boot) {
    bootRes <- boot::boot(data = data, statistic = NetANOVAStatistic, 
                          R = R, sizes = sizes)
    pvalBoot <- length(which(bootRes$t > bootRes$t0)) / R
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
