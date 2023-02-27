#' Iterative Smooth Backfitting Algorithm
#'
#' Smooth backfitting procedure for nonparametric functional additive models for density responses
#'
#' @param Ly A list of \emph{n} vectors holding the observed random samples independently drawn from each density response
#' @param X An \emph{n}-by-\emph{d} design matrix whose row vectors consist of regressors for each component function
#' @param x An \emph{M}-by-\emph{d} matrix whose row vectors consist of evaluation points for each component function
#' @param h0 A scalar bandwidth for kernel smoothing marginal mean function
#' @param h A \emph{d}-dimensional vector bandwidths for kernel smoothing each component function
#' @param densSup A \emph{2}-dimensional vector that specify the lower and upper limits of the common support for density responses. Defaults to the range of the observed random samples.
#'
#' @details \code{DensityFAM} fits component functions of functional additive models for density responses associated with multiple regressors (Han et al., 2020). 
#' The current version only supports the log-quantile density (LQD) transformation approach (Petersen and Mueller, 2016). 
#'
#' @seealso SBFitting() in the \emph{fdapace} pacakge
#'
#' @return A list holding the following fields:
#' \item{lqdGrid}{A grid where the LQD transformation of density responses is evaluated. Defaults to an equally space grid on [0,1].}
#' \item{densGrid}{A grid where density responses are evaluated. Defaults to to an equally space grid on the range of the observed random samples.}
#' \item{lqdSbfMean}{The estimates of the marginal mean function evaluated on \code{lqdGrid}}
#' \item{LlqdSbfComp}{A list of \emph{d} matrices holding the estimates of each component functions evaluated on \code{lqdGrid}}
#' \item{densSbfMean}{The density inversion of \code{lqdSbfMean} evaluated on \code{densGrid}}
#' \item{LdensSbfComp}{A list of \emph{d} matrices holding the density inversion of \code{lqdSbfMean} together with each component function evaluated on \code{lqdGrid}}
#' @examples
#' library(MASS)
#' 
#' ### generating random samples
#' g1 <- function (u, x1) sin(2*pi*u) * (2*x1 - 1)
#' g2 <- function (u, x2) sin(2*pi*u) * sin(2*pi*x2)
#' 
#' g <- function (u, x) g1(u, x[1]) + g2(u, x[2])
#' 
#' GenLqdNoise <- function (u, e) e[1]*sin(pi*u) + e[2]*sin(2*pi*u) 
#' 
#' GenQdResp <- function (u, x, e) exp(g(u, x) + GenLqdNoise(u, e))
#' 
#' GenCondQResp <- function (u, x, e) {
#'   
#'   tmp1 <- integrate(GenQdResp, lower = 0, upper = u, x = x, e = e)$value
#'   tmp2 <- integrate(GenQdResp, lower = 0, upper = 1, x = x, e = e)$value
#'   
#'   return (tmp1 / tmp2)
#' }
#' 
#' set.seed(999)
#' n <- 100
#' N <- 200
#' 
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
#' X <- pnorm(mvrnorm(n, rep(0, 2), Sigma))
#' 
#' Ly <- list()
#' for (i in 1:n) {
#'   U_i <- runif(N)
#'   E_i <- c(rnorm(1, 0, 0.1), rnorm(1, 0, 0.05))
#'   Ly[[i]] <- sapply(1:N, function (l) GenCondQResp(U_i[l], X[i,], E_i))
#' }
#' 
#' M <- 51
#' x0 <- seq(0, 1, length.out = M)
#' x <- cbind(x0, x0)
#' 
#' h0 <- 0.05
#' h <- c(0.075, 0.075)
#' densSup <- c(0, 1)
#' 
#' 
#' ### component function estimation
#' # smooth backfitting
#' estDensityFAM <- DensityFAM(Ly = Ly, X = X, x = x, h0 = h0, h = h, densSup = densSup)
#' 
#' # lqd evaluation grid
#' lqdGrid <- estDensityFAM$lqdGrid
#' lqdGridLen <- length(lqdGrid)
#' 
#' # density evaluation grid
#' densGrid <- estDensityFAM$densGrid
#' densGridLen <- length(densGrid)
#' 
#' # true LQD component functions
#' g1Eval <- g2Eval <- matrix(nrow = lqdGridLen, ncol = M)
#' for (l in 1:lqdGridLen) {
#'   for (m in 1:M) {
#'     g1Eval[l,m] <- g1(lqdGrid[l], x[m,1])
#'     g2Eval[l,m] <- g2(lqdGrid[l], x[m,2])
#'   }
#' }
#' 
#' # LQD component function estimates
#' g0Sbf <- estDensityFAM$lqdSbfMean
#' gjSbf <- estDensityFAM$LlqdSbfComp
#' 
#' # true density component functions
#' dens1Eval <- dens2Eval <- matrix(nrow = densGridLen, ncol = M)
#' for (m in 1:M) {
#'   dens1Eval[,m] <- fdadensity::lqd2dens(lqd = g1Eval[,m],
#'                                         lqSup = lqdGrid,
#'                                         dSup = densGrid)
#'   dens2Eval[,m] <- fdadensity::lqd2dens(lqd = g2Eval[,m],
#'                                         lqSup = lqdGrid,
#'                                         dSup = densGrid)
#' }
#' 
#' # density component function estimates
#' dens0Sbf <- estDensityFAM$densSbfMean
#' densjSbf <- estDensityFAM$LdensSbfComp
#' 
#' # graphical illustration of LQD component function estimates
#' par(mfrow = c(2,2))
#' par(mar=rep(0.5, 4)+0.1)
#' persp(lqdGrid, x[,1], g1Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n LQD component (g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(lqdGrid, x[,2], g2Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n LQD component  (g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(lqdGrid, x[,1], gjSbf[[1]],
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n SBF estimate (g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(lqdGrid, x[,2], gjSbf[[2]],
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x2', zlab = '\n SBF estimate (g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' # graphical illustration of density component function estimates
#' par(mfrow = c(2,2))
#' par(mar=rep(0.5, 4)+0.1)
#' persp(densGrid, x[,1], dens1Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n Component density \n (Inversion of g0 + g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(densGrid, x[,2], dens2Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n Component density \n (Inversion of g0 + g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(densGrid, x[,1], densjSbf[[1]],
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n SBF estimate \n (Inversion of g0 + g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(densGrid, x[,2], densjSbf[[2]],
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n SBF estimate \n (Inversion of g0 + g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' 
#' ### fitted density responses
#' # smooth backfitting
#' fitDensityFAM <- DensityFAM(Ly = Ly, X = X, h0 = h0, h = h, densSup = densSup)
#' 
#' # lqd evaluation grid
#' lqdGrid <- fitDensityFAM$lqdGrid
#' lqdGridLen <- length(lqdGrid)
#' 
#' # density evaluation grid
#' densGrid <- fitDensityFAM$densGrid
#' densGridLen <- length(densGrid)
#' 
#' # fitted LQD component functions
#' g0SBFit <- fitDensityFAM$lqdSbfMean
#' gjSBFit <- fitDensityFAM$LlqdSbfComp
#' 
#' # fitted density responses
#' densSBFit <- lapply(1:n, 
#'                     function (i) {
#'                       gSBFit <- g0SBFit + gjSBFit[[1]][,i] + gjSBFit[[2]][,i]
#'                       densSBFit_i <- fdadensity::lqd2dens(lqd = gSBFit, 
#'                                                           lqSup = lqdGrid,
#'                                                           dSup = densGrid)
#'                       return (densSBFit_i)
#'                     }
#' )
#' 
#' # graphical illustration of fitted density responses
#' set.seed(1)
#' ind <- sample(1:n, 12)
#' par(mfrow = c(3, 4))
#' par(mar=c(4, 4, 4, 1)+0.1)
#' for (i in ind) {
#'   hist_i <- hist(Ly[[i]], plot = FALSE)
#'   hist(Ly[[i]], probability = TRUE, 
#'        ylim = range(c(hist_i$density, densSBFit[[i]])),
#'        xlab = 'Y',
#'        main = paste(i, '-th random sample \n with X = (', round(X[i,1],2), ', ', round(X[i,2],2), ')', sep = ''))
#'   lines(densGrid, densSBFit[[i]], col = 2, lwd = 2)
#' }
#' 
#' @references
#' \cite{Han, K., Mueller, H.-G., and Park, B. U. (2020), "Additive functional regression for densities as responses", Journal of the American Statistical Association , 115 (530), pp.997-1010.}
#'
#' \cite{Peterson, A. and Mueller, H.-G. (2019), "Frechet regression for random objects with Eucliean predictors", The Annals of Statistics, 47(2), pp.691-719.}
#' 
#' \cite{Peterson, A. and Mueller, H.-G. (2016), "Functional Data Analysis for Density Functions by Transformation to a Hilbert space", The Annals of Statistics, 44(1), pp.183-218}
#'
#' @export

DensityFAM <- function (Ly, X, x = NULL, h0 = NULL, h = NULL, densSup = NULL) {
  
  if (is.list(Ly) == FALSE) {
    
    return (message('The response input should be a list of random samples.'))
    
  }
  
  if (is.null(ncol(X))) {
    
    return (message('The design matrix must be multi-dimensional.'))
    
  }
  
  if (length(Ly) != nrow(X)) {
    
    return (message('The number of subjects in Ly and X are different.'))
    
  }
  
  if (is.null(x) == TRUE) {
    
    message('The evaluation grid will be replaced by the design matrix.')
    x <- X

  } else {
    
    if (sum(apply(x, 2, diff) < 0) > 0) {
      
      message('The evaluation grid must be of increasing order. Sorted by increasing for each column.')
      x <- apply(x, 2, sort)
      
    }
    if (is.null(ncol(x))) {
      
      return (message('The evaluation grid must be multi-dimensional corresponding to the support of each additive component.'))
      
    } else if (ncol(X) != ncol(x)) {
      
      return (message('The lengths of columns between X and x are different.'))
      
    }
    
  }
  
  if (is.null(h0)) {
    
    h0 <- 0.05
    
  } else if (h0 <= 0 || h0 >= 1) {
    
    return (message('The bandwidth for smoothing transformed densities should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(h)) {
    
    h <- rep(0.25 * n^(-1/5), ncol(X)) * apply(apply(x, 2, range), 2, diff)
    
  } else if (length(h) < 2) {
    
    return (message('The bandwidth must be multi-dimensional.'))
    
  } else if (min(h) <= 0 || max(h) >= 1) {
    
    return (message('Bandwidths for additive component functions should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(densSup)) {
   
    message('The density response is assumed to have the common support. The observed min and max will be used for the lower/upper limits of the support.') 
    densSup <- range(unlist(Ly))
    
  }
  
  # common arguments
  M <- nrow(x)
  n <- nrow(X)
  d <- ncol(X)
  
  # minimum bandwidth
  hMin <- apply(apply(apply(X, 2, sort), 2, diff), 2, max)
  for (j in 1:d) {
    if (hMin[j]/2 > h[j]) {
      h[j] <- hMin[j]/2
    }
  }
  
  # common evaluation grids for density and LQD responses
  densGridLen <- 101
  lqdGridLen <- 201
  densGrid <- seq(0, 1, length.out = densGridLen)
  lqdGrid <- seq(0, 1, length.out = lqdGridLen)
  
  # normalize random samples to be supported on [0,1]
  normalizeLy <- lapply(1:n, 
                          function (i) {
                            normalizeY_i <- (Ly[[i]] - densSup[1]) / diff(densSup)
                          }
                        )
  
  # density response reconstruction
  message('Estimating density resopnses...')
  Ldens <- lapply(1:n,
                function (i) {
                  f_i <- fdadensity::CreateDensity(y = normalizeLy[[i]], 
                                                   optns = list(outputGrid = densGrid))$y  
                  return (f_i)
                }
               )
  
  densMat <- matrix(unlist(Ldens), nrow = n, ncol = densGridLen, byrow = TRUE)
  
  # LQD transformation
  message('Transforming density resopnses...')
  Llqd <- lapply(1:n, 
                  function (i) {
                    
                    f_i <- Ldens[[i]]
                    if (min(f_i) < 1e-8) {
                      
                      f_i <- fdadensity::RegulariseByAlpha(densGrid, f_i)
                      
                    }
                    
                    lqd_i <- fdadensity::dens2lqd(dens = f_i, 
                                                  dSup = densGrid,
                                                  lqSup = lqdGrid)
                    return (lqd_i)
                  }
                 )
  
  lqdMat <- matrix(unlist(Llqd), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smoothing LQD responses
  LlqdSmooth <- lapply(1:n, 
                        function (i) {
                          lqdSmooth_i <- fdapace::Lwls1D(bw = h0,
                                                         kernel_type = 'epan',
                                                         xin = lqdGrid,
                                                         yin = Llqd[[i]],
                                                         xout = lqdGrid)
                          return (lqdSmooth_i)
                        }
                       )
  
  lqdSmoothMat <- matrix(unlist(LlqdSmooth), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smooth backfitting
  message('Smooth backfitting...')
  g0Sbf <- c()
  gjSbf <- lapply(1:d,
                    function (j) {
                      return (matrix(NA, nrow = lqdGridLen, ncol = M))
                    }
                  ) 
  gjSbfAddMean <- lapply(1:d,
                         function (j) {
                           return (matrix(NA, nrow = lqdGridLen, ncol = M))
                         }
  ) 
  
  
  for (l in 1:lqdGridLen) {
    
    sbfSurf <- fdapace::SBFitting(lqdSmoothMat[,l], x, X, h = h)
    
    g0Sbf[l] <- sbfSurf$mY
    for (j in 1:d) {
      
      gjSbf[[j]][l,] <- sbfSurf$SBFit[,j]
      gjSbfAddMean[[j]][l,] <- sbfSurf$mY + sbfSurf$SBFit[,j]
      
    }
  }
  
  # LQD inversion to density
  message('Inverting to density resopnses...')
  dens0Sbf <- fdadensity::lqd2dens(lqd = g0Sbf,
                                   lqSup = lqdGrid,
                                   dSup = densGrid)
  densjSbf <- lapply(1:d, 
                      function (j) {
                        return (matrix(NA, nrow = densGridLen, ncol = M))  
                      }
                     )
  
  for (m in 1:M) {
    for (j in 1:d) {
      
      densSbf_j <- fdadensity::lqd2dens(lqd = gjSbfAddMean[[j]][, m],
                                        lqSup = lqdGrid,
                                        dSup = densGrid)
      
      densjSbf[[j]][, m] <- fdapace::Lwls1D(bw = h0,
                                            kernel_type = 'epan',
                                            xin = densGrid,
                                            yin = densSbf_j,
                                            xout = densGrid)
    }
  }
  
  return (list(Ly = Ly,
               X = X,
               x = x,
               h0 = h0,
               h = h,
               densSup = densSup,
               lqdGrid = lqdGrid,
               densGrid = (densGrid*diff(densSup) + densSup[1]),
               lqdSbfMean = g0Sbf,
               LlqdSbfComp = gjSbf,
               densSbfMean = dens0Sbf,
               LdensSbfComp = densjSbf
               )
          )
}



