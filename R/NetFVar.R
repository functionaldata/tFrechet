#' @title Fréchet Variance for Networks
#' @description Obtain Fréchet variance for graph Laplacian matrices, 
#'   covariance matrices, or correlation matrices 
#'   with respect to the Frobenius distance.
#' @param Ly A list (length n) of m by m matrices or a m by m by n array where
#'   \code{Ly[, , i]} contains an m by m matrix, which can be either graph 
#'   Laplacian matrices or covariance matrices or correlation matrices.
#' @return A list containing the following fields:
#' \item{NetFVar}{A scalar holding the Fréchet variance.}
#' \item{NetFMean}{A matrix holding the Fréchet mean.}
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
    if (is.array(Ly)) {
      Ly <- lapply(seq(dim(Ly)[3]), function(i) Ly[, , i])
    } else {
      stop("Ly must be a list or an array")
    }
  }
  if (length(unique(sapply(Ly, length))) > 1) {
    stop("each matrix in Ly should be of the same dimension")
  }
  if (any(sapply(Ly, function(Lyi) nrow(Lyi) != ncol(Lyi)))) {
    stop("each matrix in Ly should be a square matrix")
  }
  n <- length(Ly)
  mup <- rowMeans(matrix(unlist(Ly), ncol = n))
  Vp <- mean(sapply(Ly, function(Lyi) sum((Lyi - mup)^2)))
  res <- list(NetFVar = Vp, NetFMean = mup)
  res
}
