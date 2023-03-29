#' @title Additive functional regression for densities as responses
#'
#' @description Smooth backfitting procedure for estimating nonparametric functional additive models for density responses
#'
#' @param Ly A list of \emph{n} vectors holding random samples independently drawn from each density response
#' @param X An \emph{n}-by-\emph{d} design matrix whose row vectors consist of regressors for component functions
#' @param x A grid matrix whose row vectors consist of evaluation points for component functions. Defaults to the observed design matrix.
#' @param hu A scalar bandwidth for kernel smoothing marginal mean function of the functional additive model
#' @param huA \emph{d}-dimensional vector bandwidth for kernel smoothing each component function
#' @param dSup A \emph{2}-dimensional vector that specify the lower and upper limits of the common support for density responses. Defaults to the range of the observed random samples.
#'
#' @details \code{AddDensReg} fits additive functional regression models for density responses, where the conditional mean function of a transformed density response is given by the summation of nonparametric univariate functions associated with \emph{d} covariates, respectively.
#' Instead of having functional responses as infinite-dimensional random objects, \code{AddDensReg} has inputs with random samples independently drawn from each random density, and density responses are reconstructed by kernel density estimation.
#' The current version only supports the log-quantile density (LQD) transformation proposed by Petersen and Mueller (2016).
#'
#' @return A list holding the following fields:
#' \item{lqdGrid}{A grid where the LQD transformed density responses are evaluated. Defaults to an equally space grid over \code{dSup}.}
#' \item{lqdSbfMean}{The marginal mean function estimates evaluated on \code{lqdGrid}}
#' \item{LlqdSbfComp}{A list of \emph{d} matrices holding the estimates of each component function evaluated on \code{lqdGrid}}
#' \item{densGrid}{A grid where density responses are evaluated. Defaults to an equally space grid over the range of the observed random samples.}
#' \item{densSbfMean}{The density inversion of \code{lqdSbfMean} evaluated on \code{densGrid}}
#' \item{LdensSbfComp}{A list of \emph{d} matrices holding the density inversion of each component function together with \code{lqdSbfMean} evaluated on \code{lqdGrid}}
#'
#' @references
#' \cite{Han, K., Mueller, H.-G., and Park, B. U. (2020), "Additive functional regression for densities as responses", Journal of the American Statistical Association , 115 (530), pp.997-1010.}
#'
#' \cite{Peterson, A. and Mueller, H.-G. (2016), "Functional data analysis for density functions by transformation to a Hilbert space", The Annals of Statistics, 44(1), pp.183-218}
#'
#' @seealso \code{SBFitting} in the \code{fdapace} package
#'
#' @examples
#' library(MASS)
#' 
#' # additive component functions
#' g1 <- function (u, x1) sin(2*pi*u) * (2*x1 - 1)
#' g2 <- function (u, x2) sin(2*pi*u) * sin(2*pi*x2)
#' 
#' g <- function (u, x) g1(u, x[1]) + g2(u, x[2])
#' 
#' # generating random samples from conditional quantile functions
#' GenLqdNoise <- function (u, e) e[1]*sin(pi*u) + e[2]*sin(2*pi*u) 
#' GenQdensResp <- function (u, x, e) exp(g(u, x) + GenLqdNoise(u, e))
#' 
#' GenQuantileResp <- function (u, x, e) {
#'   
#'   tmp1 <- integrate(GenQdensResp, lower = 0, upper = u, x = x, e = e)$value
#'   tmp2 <- integrate(GenQdensResp, lower = 0, upper = 1, x = x, e = e)$value
#'   
#'   return (tmp1 / tmp2)
#' }
#' 
#' set.seed(999)
#' n <- 150
#' N <- 250
#' 
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
#' X <- pnorm(mvrnorm(n, rep(0, 2), Sigma))
#' 
#' Ly <- list()
#' for (i in 1:n) {
#'   U_i <- runif(N)
#'   E_i <- c(rnorm(1, 0, 0.1), rnorm(1, 0, 0.05))
#'   Ly[[i]] <- sapply(1:N, function (l) GenQuantileResp(U_i[l], X[i,], E_i))
#' }
#' 
#' M <- 51
#' x1 <- x2 <- seq(0, 1, length.out = M)
#' x <- cbind(x1, x2)
#' 
#' hu <- 0.05
#' hx <- c(0.075, 0.075)
#' dSup <- c(0, 1)
#' 
#' # estimating the functional additive model
#' estAddDensReg <- AddDensReg(Ly = Ly, X = X, x = x, hu = hu, hx = hx, dSup = dSup)
#' 
#' # true LQD component functions
#' g1Eval <- g2Eval <- matrix(nrow = length(estAddDensReg$lqdGrid), ncol = M)
#' for (l in seq(estAddDensReg$lqdGrid)) {
#'   for (m in seq(M)) {
#'     g1Eval[l,m] <- g1(estAddDensReg$lqdGrid[l], x1[m])
#'     g2Eval[l,m] <- g2(estAddDensReg$lqdGrid[l], x2[m])
#'   }
#' }
#' 
#' # LQD component function estimates
#' g0Sbf <- estAddDensReg$lqdSbfMean
#' gjSbf <- estAddDensReg$LlqdSbfComp
#' 
#' # true density component functions
#' dens1Eval <- dens2Eval <- matrix(nrow = length(estAddDensReg$densGrid), ncol = M)
#' for (m in seq(M)) {
#'   dens1Eval[,m] <- fdadensity::lqd2dens(lqd = g1Eval[,m],
#'                                         lqdSup = estAddDensReg$lqdGrid,
#'                                         dSup = estAddDensReg$densGrid)
#'   dens2Eval[,m] <- fdadensity::lqd2dens(lqd = g2Eval[,m],
#'                                         lqdSup = estAddDensReg$lqdGrid,
#'                                         dSup = estAddDensReg$densGrid)
#' }
#' 
#' # density component function estimates
#' dens0Sbf <- estAddDensReg$densSbfMean
#' densjSbf <- estAddDensReg$LdensSbfComp
#' 
#' # graphical illustration of LQD component function estimates
#' par(mfrow = c(2,2))
#' par(mar=rep(0.5, 4)+0.1)
#' persp(estAddDensReg$lqdGrid, x1, g1Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n LQD component (g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$lqdGrid, x2, g2Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n LQD component  (g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$lqdGrid, x1, gjSbf[[1]],
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x1', zlab = '\n SBF estimate (g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$lqdGrid, x2, gjSbf[[2]],
#'       theta = 35, phi = 35,
#'       xlab = '\n u', ylab = '\n x2', zlab = '\n SBF estimate (g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' # graphical illustration of density component function estimates
#' par(mfrow = c(2,2))
#' par(mar=rep(0.5, 4)+0.1)
#' persp(estAddDensReg$densGrid, x1, dens1Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n Component density \n (Inversion of g0 + g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$densGrid, x2, dens2Eval,
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n Component density \n (Inversion of g0 + g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$densGrid, x1, densjSbf[[1]],
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n SBF estimate \n (Inversion of g0 + g1)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' persp(estAddDensReg$densGrid, x2, densjSbf[[2]],
#'       theta = 35, phi = 35,
#'       xlab = '\n y', ylab = '\n x1', zlab = '\n\n SBF estimate \n (Inversion of g0 + g2)',
#'       border = NA, shade = 0.5,
#'       ticktype = 'detailed')
#' 
#' # fitted density responses
#' fitAddDensReg <- AddDensReg(Ly = Ly, X = X, hu = hu, hx = hx, dSup = dSup)
#' 
#' # fitted LQD component functions
#' g0SBFit <- fitAddDensReg$lqdSbfMean
#' gjSBFit <- fitAddDensReg$LlqdSbfComp
#' 
#' # fitted density responses
#' densSBFit <- lapply(1:n, 
#'                     function (i) {
#'                       gSBFit <- g0SBFit + gjSBFit[[1]][,i] + gjSBFit[[2]][,i]
#'                       densSBFit_i <- fdadensity::lqd2dens(lqd = gSBFit, 
#'                                                           lqdSup = fitAddDensReg$lqdGrid,
#'                                                           dSup = fitAddDensReg$densGrid)
#'                       return (densSBFit_i)
#'                     }
#' )
#' 
#' # graphical illustration of fitted density responses
#' set.seed(999)
#' ind <- sample(1:n, 12)
#' par(mfrow = c(3, 4))
#' par(mar=c(4, 4, 4, 1)+0.1)
#' for (i in ind) {
#'   hist_i <- hist(Ly[[i]], plot = FALSE)
#'   hist(Ly[[i]], probability = TRUE, 
#'        ylim = range(c(hist_i$density, densSBFit[[i]])),
#'        xlab = 'Y',
#'        main = paste(i, '-th random sample \n with X = (', round(X[i,1],2), ', ', round(X[i,2],2), ')', sep = ''))
#'   lines(fitAddDensReg$densGrid, densSBFit[[i]], col = 2, lwd = 2)
#' }
#' 
#' @export

AddDensReg <- function (Ly, X, x = NULL, hu = NULL, hx = NULL, dSup = NULL) {
  
  if (is.list(Ly) == FALSE) {
    
    return (message('The response input should be a list of random samples.'))
    
  }
  
  if (is.null(ncol(X))) {
    
    return (message('The design matrix must be multi-dimensional corresponding to additive component.'))
    
  }
  
  if (length(Ly) != nrow(X)) {
    
    return (message('The sample sizes of responses and regressors are different.'))
    
  }
  
  if (is.null(x) == TRUE) {
    
    message('The evaluation grid will be replaced by the observed design matrix.')
    x <- X

  } else {
    
    # if (sum(apply(x, 2, diff) < 0) > 0) {
    #   
    #   message('The evaluation grid must be in increasing order. Sorted by increasing for each column.')
    #   x <- apply(x, 2, sort)
    #   
    # }
    
    if (is.null(ncol(x))) {
      
      return (message('The evaluation grid must be multi-dimensional corresponding to additive component.'))
      
    } else if (ncol(X) != ncol(x)) {
      
      return (message('The lengths of columns between X and x are different.'))
      
    }
    
  }
  
  if (is.null(hu)) {
    
    hu <- 0.05
    
  } else if (hu <= 0 || hu >= 1) {
    
    return (message('The bandwidth for smoothing transformed densities should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(hx)) {
    
    hx <- rep(0.25 * n^(-1/5), ncol(X)) * apply(apply(x, 2, range), 2, diff)
    
  } else if (length(hx) < 2) {
    
    return (message('The bandwidth must be multi-dimensional.'))
    
  } else if (min(hx) <= 0 || max(hx) >= 1) {
    
    return (message('Bandwidths for additive component functions should be chosen by a positive number less than 1.'))
    
  }
  
  if (is.null(dSup)) {
   
    message('The density response is assumed to have the common support. The observed min and max will be used for the lower/upper limits of the support.') 
    dSup <- range(unlist(Ly))
    
  }
  
  # common arguments
  M <- nrow(x)
  n <- nrow(X)
  d <- ncol(X)
  
  # minimum bandwidth
  hxMin <- apply(apply(apply(X, 2, sort), 2, diff), 2, max)
  for (j in 1:d) {
    if (hx[j] < hxMin[j]) {
      hx[j] <- hxMin[j]
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
                            normalizeY_i <- (Ly[[i]] - dSup[1]) / diff(dSup)
                          }
                        )
  
  # density response reconstruction
  message('Estimating density resopnses...(1/4)')
  Ldens <- lapply(1:n,
                function (i) {
                  f_i <- frechet::CreateDensity(y = normalizeLy[[i]], 
                                                   optns = list(outputGrid = densGrid))$y  
                  return (f_i)
                }
               )
  
  densMat <- matrix(unlist(Ldens), nrow = n, ncol = densGridLen, byrow = TRUE)
  
  # LQD transformation
  message('Transforming density resopnses...(2/4)')
  Llqd <- lapply(1:n, 
                  function (i) {
                    
                    f_i <- Ldens[[i]]
                    if (min(f_i) < 1e-8) {
                      
                      f_i <- fdadensity::RegulariseByAlpha(densGrid, f_i)
                      
                    }
                    
                    lqd_i <- fdadensity::dens2lqd(dens = f_i, 
                                                  dSup = densGrid,
                                                  lqdSup = lqdGrid)
                    return (lqd_i)
                  }
                 )
  
  lqdMat <- matrix(unlist(Llqd), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smoothing LQD responses
  LlqdSmooth <- lapply(1:n, 
                        function (i) {
                          lqdSmooth_i <- fdapace::Lwls1D(bw = hu,
                                                         kernel_type = 'epan',
                                                         xin = lqdGrid,
                                                         yin = Llqd[[i]],
                                                         xout = lqdGrid)
                          return (lqdSmooth_i)
                        }
                       )
  
  lqdSmoothMat <- matrix(unlist(LlqdSmooth), nrow = n, ncol = lqdGridLen, byrow = TRUE)
  
  # smooth backfitting
  message('Smooth backfitting...(3/4)')
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
    
    sbfSurf <- fdapace::SBFitting(lqdSmoothMat[,l], x, X, h = hx)
    
    g0Sbf[l] <- sbfSurf$mY
    for (j in 1:d) {
      
      gjSbf[[j]][l,] <- sbfSurf$SBFit[,j]
      gjSbfAddMean[[j]][l,] <- sbfSurf$mY + sbfSurf$SBFit[,j]
      
    }
  }
  
  # LQD inversion to density
  message('Inverting to density resopnses...(4/4)')
  dens0Sbf <- fdadensity::lqd2dens(lqd = g0Sbf,
                                   lqdSup = lqdGrid,
                                   dSup = densGrid)
  densjSbf <- lapply(1:d, 
                      function (j) {
                        return (matrix(NA, nrow = densGridLen, ncol = M))  
                      }
                     )
  
  for (m in 1:M) {
    for (j in 1:d) {
      
      densSbf_j <- fdadensity::lqd2dens(lqd = gjSbfAddMean[[j]][, m],
                                        lqdSup = lqdGrid,
                                        dSup = densGrid)
      
      densjSbf[[j]][, m] <- fdapace::Lwls1D(bw = hu,
                                            kernel_type = 'epan',
                                            xin = densGrid,
                                            yin = densSbf_j,
                                            xout = densGrid)
    }
  }
  
  return (list(Ly = Ly,
               X = X,
               x = x,
               hu = hu,
               hx = hx,
               dSup = dSup,
               lqdGrid = lqdGrid,
               densGrid = (dSup[1] + densGrid*diff(dSup)),
               lqdSbfMean = g0Sbf,
               LlqdSbfComp = gjSbf,
               densSbfMean = dens0Sbf,
               LdensSbfComp = densjSbf
               )
          )
}
