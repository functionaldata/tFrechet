#' frechet: Statistical Analysis for Random Objects and Non-Euclidean Data
#' @description Provides implementation of statistical methods for random objects 
#' lying in various metric spaces, which are not necessarily linear spaces. 
#' The core of this package is Fréchet regression for random objects with 
#' Euclidean predictors, which allows one to perform regression analysis 
#' for non-Euclidean responses under some mild conditions. 
#' Examples include distributions in \eqn{L^2}-Wasserstein space, 
#' covariance matrices endowed with power metric (with Frobenius metric as a special case), Cholesky and log-Cholesky metrics.  
#' References: Petersen, A., & Müller, H.-G. (2019) <doi:10.1214/17-AOS1624>.
#' @docType package
#' @name frechet
#' @importFrom boot boot
#' @importFrom e1071 rbridge
#' @importFrom fdadensity dens2quantile
#' @importFrom grDevices colorRampPalette dev.new palette
#' @importFrom graphics abline axis barplot boxplot grid hist layout legend lines matlines matplot par plot points polygon rect text
#' @importFrom igraph laplacian_matrix sample_pa
#' @importFrom methods as
#' @importFrom Matrix Matrix
#' @importFrom pracma trapz
#' @importFrom stats aggregate approx approxfun binomial cov cor density dist dnorm dunif fitted glm kmeans lm median na.omit optim optimize poly predict quantile rnorm runif spline var sd weighted.mean
#' @importFrom utils head tail
#' @importFrom stats pchisq
NULL

utils::globalVariables(c("y"))
