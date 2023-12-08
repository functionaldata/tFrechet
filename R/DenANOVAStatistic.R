#' @description Helper function computing bootstrap replications of
#'   the Fréchet ANOVA test statistics for densities/quantiles.
#' @importFrom pracma trapz
#' @noRd

DenANOVAStatistic <- function(data, indices, sizes, qSup) {
  k <- length(sizes) # number of groups
  n <- sum(sizes) # total number of observations
  LyBoot <- data[indices] # booted sample
  lambda <- sizes / n # ratio for each group
  groupData <- split(LyBoot, rep(1:k, sizes))
  mup <- rowMeans(matrix(unlist(LyBoot), ncol = n))
  Vp <- mean(sapply(LyBoot, function(LyBooti) {
    pracma::trapz(qSup, (LyBooti - mup)^2)
  }))
  V <- rep(0, k)
  sigma2 <- rep(0, k)
  for (i in 1:k) {
    mui <- rowMeans(matrix(unlist(groupData[[i]]), ncol = sizes[i]))
    Di <- sapply(groupData[[i]], function(Lyi) {
      pracma::trapz(qSup, (Lyi - mui)^2)
    })
    V[i] <- mean(Di)
    sigma2[i] <- mean(Di^2) - (mean(Di))^2
  }
  Fn <- Vp - sum(lambda * V)
  Un <- 0
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      Un <- Un + lambda[i] * lambda[j] * (V[i] - V[j])^2 / (sigma2[i] * sigma2[j])
    }
  }
  Tn <- n * Un / sum(lambda / sigma2) + n * Fn^2 / sum(lambda^2 * sigma2)
  Tn
}
