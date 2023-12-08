#' @description Helper function computing bootstrap replications of
#'   the Fréchet CPD test statistics for densities/quantiles.
#' @importFrom pracma trapz
#' @noRd

DenCPDStatistic <- function(data, indices, cutOff, qSup, bootSize = length(indices)) {
  LyBoot <- data[indices[1:bootSize]] # booted sample. m see bootstrap scheme in the paper
  n <- length(LyBoot) # number of observations
  scope <- ceiling(cutOff * n):(n - ceiling(cutOff * n))
  nTn <- sapply(scope, function(i) {
    n1 <- i
    n2 <- n - i
    lambda <- n1 / n
    LyBoot1 <- LyBoot[1:n1]
    LyBoot2 <- LyBoot[(n1 + 1):n]
    mup <- rowMeans(matrix(unlist(LyBoot), ncol = n)) # overall Frechet mean
    mu1 <- rowMeans(matrix(unlist(LyBoot1), ncol = n1)) # first part Frechet mean
    mu2 <- rowMeans(matrix(unlist(LyBoot2), ncol = n2)) # second part Frechet mean
    V1 <- mean(sapply(LyBoot1, function(LyBoot1i) {
      pracma::trapz(qSup, (LyBoot1i - mu1)^2)
    }))
    V2 <- mean(sapply(LyBoot2, function(LyBoot2i) {
      pracma::trapz(qSup, (LyBoot2i - mu2)^2)
    }))
    V1C <- mean(sapply(LyBoot1, function(LyBoot1i) {
      pracma::trapz(qSup, (LyBoot1i - mu2)^2)
    }))
    V2C <- mean(sapply(LyBoot2, function(LyBoot2i) {
      pracma::trapz(qSup, (LyBoot2i - mu1)^2)
    }))
    Di <- sapply(LyBoot, function(LyBooti) {
      pracma::trapz(qSup, (LyBooti - mup)^2)
    })
    sigma2 <- mean(Di^2) - mean(Di)^2
    nTni <- n * lambda * (1 - lambda) * ((V1 - V2)^2 + (V1C - V1 + V2C - V2)^2) / sigma2 # formula 2.5 nTn
    nTni
  })
  c(maxnTn = max(nTn), tau = which.max(nTn)) # test statistic and location of change point
}
