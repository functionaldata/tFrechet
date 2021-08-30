#' @title Fréchet Variance for Networks
#' @description Obtain Fréchet variance for graph Laplacian matrices, 
#'   covariance matrices, or correlation matrices 
#'   with respect to the Frobenius distance.
#' @param Ly a list (length \eqn{n}) of \eqn{m} by \eqn{m} matrices, 
#'   which can be either graph Laplacian matrices or covariance matrices 
#'   or correlation matrices. 
#' @return A list containing the following fields:
#' \item{NetFVar}{a scalar holding the Fréchet variance.}
#' \item{NetFMean}{a matrix holding the Fréchet mean.}
#' @examples
#' set.seed(1)
#' n <- 100
#' U <- pracma::randortho(10)
#' Ly <- lapply(1:n, function(i) {
#'   U %*% diag(rexp(10, (1:10)/2)) %*% t(U)
#' })
#' res <- NetFVar(Ly)
#' res$NetFVar
#' @export

NetFVar <- function(Ly = NULL) {
  if (is.null(Ly)) {
    stop("requires the input of Ly")
  }
  if (!is.list(Ly)) {
    stop("Ly should be a list")
  }
  if (length(unique(sapply(Ly, length))) > 1) {
    stop("each matrix in Ly should be of the same dimension")
  }
  if (any(sapply(Ly, function(Lyi) nrow(Lyi) != col(Lyi)))) {
    stop("each matrix in Ly should be a square matrix")
  }
  n <- length(Ly)
  mup <- rowMeans(matrix(unlist(Ly), ncol = n))
  Vp <- mean(sapply(Ly, function(Lyi) sum((Lyi - mup)^2)))
  res <- list(NetFVar = Vp, NetFMean = mup)
  res
}
