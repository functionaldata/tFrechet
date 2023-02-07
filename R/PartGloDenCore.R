#' @title Partially global concurrent object regresion (CORE)
#' @description Concurrent object regression (CORE) model for time-varying density responses and time-varying real covariates, by modeling
#' the global dependence of the response on the predictor and the local dependence on the time direction
#' through weighted Frechet means.
#' @param xin A list holding the time-varying real predictors. 
#' Each element of \code{xin} is a matrix where the rows hold the predictor values at each time point in the corresponding element of \code{tin} and
#' the number of columns is same as the predictor dimension.
#' @param tin A list holding the time points at which the response and predictors are observed.
#' Each element of \code{tin} is a vector holding the time points for each subject.
#' @param qin A list holding the quantile functions of the time-varying response. 
#' Each element of \code{qin} is a matrix where the rows holds the quantile function at  each time point in the corresponding element of \code{tin} and 
#' the number of columns is same as the length of the support for the quantile functions-
#' the support of the quantile functions should be the same (i.e., \code{optns$qSup}).
#' jth row of the ith element of \code{qin} qin holds the quantile response corresponding to jth row of the ith element of \code{xin}.
#' If the quantile functions are evaluated on different grids (i.e., \code{optns$diff_qSup} == TRUE}),
#' \code{qin} should be a list of list holding the time-varying quantile functions corresponding to \code{xin} at the time point \cdoe{tin},
#' each element of the sub-list consisting of two components \code{x} and \code{y} holding the support grid and the corresponding values of the quantile functions, respectively.
#' @param yin A list holding the sample of observations of the time-varying response.
#' Each element of \code{yin} is a matrix where the rows holds the response observations at  each time point in the corresponding element of \code{tin}.
#' @param hin A list holding the histograms of the  time-varying response.
#' #' Each element of \code{hin} is a list of histograms of the response observations at each time point in the corresponding element of \code{tin}. the response observations at  each time point in the corresponding element of \code{tin}.
#' Note that only one of the three \code{yin}, \code{hin}, and \code{qin} needs to be input.
#' If more than one of them are specified, \code{yin} overwrites \code{hin}, and \code{hin} overwrites \code{qin}.
#' @param xout is either NULL or is a vector with the same dimension as any element of the list \code{xin}
#' if \code{xout} is NULL, the model is fitted at each point of \code{xin}.
#' @param tout is either NULL or is a vector with the same dimension as any element of the list \code{tin}
#' if \code{xout} is NULL, the model is fitted at each point of \code{tin}.
#' @details Available control options are
#' \describe{
#' \item{bw_t}{A scalar used as the bandwidth for the local regression fit in the time direction or \code{"CV"} (default), i.e., a data-adaptive selection done by cross-validation.}
#' \item{kernelReg}{A character holding the type of kernel functions for local regression in the time direction; \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' \item{qSup}{A numeric vector holding the grid on [0,1] quantile functions take value on, if the time-varying responses have same support. Default is an equidistant grid.}
#' \item{nqSup}{A scalar giving the length of \code{qSup}. Default is 101.}
#' \item{lower}{A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.}
#' \item{upper}{A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.}
#' \item{bwRange}{A vector of length 2 containing the bandwidth selection range for the time direction \code{tin} for the case when \code{bwReg} equals \code{"CV"}. Default is \code{NULL} and is automatically chosen by a data-adaptive method.}
#' \item{bwDen}{The bandwidth value used in \code{CreateDensity()} for density estimation; positive numeric - default: determine automatically based on the data-driven bandwidth selector proposed by Sheather and Jones (1991).}
#' \item{ndSup}{The number of support points the kernel density estimation uses in \code{CreateDensity()}; numeric - default: 101.}
#' \item{dSup}{User defined output grid for the support of kernel density estimation used in \code{CreateDensity()}, it overrides \code{nRegGrid}; numeric - default: \code{NULL}}
#' \item{delta}{The size of the bin to be used used in \code{CreateDensity()}; numeric - default: \code{diff(range(y))/1000}. It only works when the raw sample is available.}
#' \item{kernelDen}{A character holding the type of kernel functions used in \code{CreateDensity()} for density estimation; \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} - default: \code{"gauss"}.}
#' \item{infSupport}{logical if we expect the distribution to have infinite support or not, used in \code{CreateDensity()} for density estimation; logical - default: \code{FALSE}}
#' \item{denLowerThreshold}{\code{FALSE} or a positive value giving the lower threshold of the densities used in \code{CreateDensity()}; default: \code{0.001 * mean(sapply(1:length(qin), function(ind){ ni = nrow(qin[[ind]]); mean(sapply(1:ni, function(l){  qin[[ind]][l,length(qin[[ind]][l,])] - qin[[ind]][l,1]}))}))  qin[,1])}.}
#' \item{diff_qSup} {logical if we expect the distributions to have different support. logical - default: \code{FALSE}}
#' }
#' @return A list containing the following components:
#' \item{xout}{Input \code{xout}. If \code{xout} was NULL return \code{xin}}
#' \item{xout}{Input \code{xout}. If \code{tout} was NULL return \code{tin}}
#' \item{dout}{A matrix or list holding the output densities corresponding to \code{xout} and \code{tout}.
#' \code{dout} list of list holding the time-varying density functions corresponding to \code{xout} at the time point \cdoe{tout},
#' each element of the sub-list consisting of two components \code{x} and \code{y}, giving the domain grid and density function values, respectively.}
#' \item{dSup}{A numeric vector giving the domain grid of each sub-list of the list \code{dout}}
#' \item{qout}{A list holding the quantile functions of the output densities corresponding to a value in \code{xout},
#' Each element of \code{qout} is a  matrix with rows corresponding to a value in \code{tout}.}
#' \item{qSup}{A numeric vector giving the domain grid of \code{qout}.}
#' \item{xin}{Input \code{xin}.}
#' \item{din}{Densities corresponding to the input \code{yin}, \code{hin} or \code{qin} for each suject and each time point of observation}
#' \item{qin}{Quantile functions corresponding to the input \code{yin}, \code{hin} or \code{qin}.}
#' \item{optns}{A list of control options used.}
#' @examples
#' \donttest{
#'  n =  50; m = 100; K = 10; dim = 2
#'  e_val_fn = function(x,k){
#'    phi1 = function(x) -cos(pi*x/10)/sqrt(5) First eigenfunction from PACE paper. Domain is [0,10]
#'    phi_k = function(x, k) sin((2*k - 1)*pi*x/10)/sqrt(5)
#'    mu = function(x) x+sin(x); mean function
#'    lamb1 = seq(8,1,length.out = 7)
#'    lamb = .7^(k-1)
#'    return(list(phi1 = phi1, phi_k = phi_k, mu = mu, lamb1 = lamb1, lamb = lamb))
#'  }
#'  generateData_K = function(n, K){
#'    xi = matrix(0, nrow = n, ncol = K)
#'    X_tilde = list() 
#'    T = list()
#'    lamb1 = e_val_fn(1,1)$lamb1
#'    lambda = c(lamb1, sapply((length(lamb1)+1) :K, function(j) e_val_fn(x,j)$lamb))
#'    Lambda = diag(lambda)
#'    for(i in 1:n){
#'      N = sample(2:3,1, replace = FALSE)
#'      T[[i]] = runif(N,0,1)
#'      X_tilde[[i]]  = replicate(dim,{
#'        L = rnorm(K,0,1)
#'        xi[i,] = sqrt(Lambda)%*%L  rnorm(dim(Lambda)[1]) Generate true scores
#'        phi = matrix(0, nrow = N, ncol = K)
#'        phi[,1] = sapply(T[[i]], function(x) e_val_fn(x, K)$phi1(x))
#'        for( j in 2:K){
#'          phi[,j] = sapply(T[[i]], function(x) e_val_fn(x, j)$phi_k(x, j))
#'        }
#'        mu = sapply(T[[i]], function(x) e_val_fn(x, K)$mu(x))
#'        X1 = (mu + phi %*%(xi[i,]))
#'      }, simplify = TRUE)
#'    }
#'    return(list(X_tilde = X_tilde , T = T))
#'  }
#'   calculates true quantile model at given x and t values
#'  true_reg = function(x,t){
#'    zeta_xt = .1 + .2*x + .5*t^2 
#'    nu_xt = .6 + .2*sin(10*pi*t)
#'    reg =  zeta_xt+  nu_xt * qnorm(seq(0.01,.99,length.out = m))
#'    reg = sort(reg,decreasing = FALSE)
#'    return(reg)
#'  }
#'   function to generate quantiles
#'  gen_quantile = function(x,t){
#'    mu_given = rnorm(n = 1, mean = .1 + .2%*%x + .5*t^2 , .1)
#'    sig_given = rgamma(1,shape = ((.6 + .2*sin(10*pi*t))^2/.1), scale = (.1/(.6 + .2*sin(10*pi*t))))
#'    Q_Y = mu_given + sig_given*qnorm(seq(0.01,.99,length.out = m))
#'    dens_Y = qf2pdf(Q_Y)
#'    return(Q_Y)
#'  }
#'  data = generateData_K(n,K)
#'  gen_pred = data$X_tilde
#'  xin = t(sapply(1:n, function(i){
#'    gen_pred[[i]]
#'  }))
#'  tin = data$T
#'  qin = lapply(1:n, function(i){
#'    ni = nrow(gen_pred[[i]])
#'    t(sapply(1:ni, function(j){
#'      gen_quantile(gen_pred[[i]][j,], tin[[i]][j])
#'    }))
#'  })
#'  xout = rowMeans(as.matrix(sapply(1:n, function(ind) colMeans(xin[[ind]]))))
#'  tout = mean(sapply(1:n, function(ind) mean(tin[[ind]])))
#'  qSup <- seq(0,1,length.out = m)
#'  optns = list(bw_t = .1, kernelReg = "gauss", lower = NULL,upper = NULL,
#'               qSup =qSup,nqSup =  100, bwRange = NULL, diff_qSup = FALSE)
#'  res = partial_glob_den_CORE(xin, tin, qin = qin, xout = xout, tout = tout, optns = optns)
#'  plot(res$qout)
#'  plot(res$dout$x, res$dout$y)
#'  nobs = 100
#'  gen_yin = function(x,t){
#'   mu_given = rnorm(n = 1, mean = .1 + .2%*%x + .5*t^2 , .1)
#'   sig_given = rgamma(1,shape = ((.6 + .2*sin(10*pi*t))^2/.1), scale = (.1/(.6 + .2*sin(10*pi*t))))
#'   Y = rnorm(nobs,mu_given ,sig_given)
#'   dens_Y = qf2pdf(Q_Y)
#'   return(Y)
#' }
#' yin = lapply(1:n, function(i){
#'   ni = nrow(xin[[i]])
#'   t(sapply(1:ni, function(j){
#'     gen_yin(xin[[i]][j,], tin[[i]][j])
#'   }))
#' })
#' hin = lapply(1:n, function(i){
#'   ni = nrow(yin[[i]])
#'   lapply(1:ni, function(j){
#'     hist(yin[[i]][j,], breaks = 50)
#'   })
#' })
#' optns = list(bw_t = .1, kernelReg = "gauss", lower = NULL,upper = NULL, qSup =qSup)
#' res = partial_glob_den_CORE(xin, tin, yin = yin, optns = optns)
#' plot(res$qout[[1]][2,])
#' plot(res$dout[[1]][[2]]$x, res$dout[[1]][[2]]$y)
#' res = partial_glob_den_CORE(xin, tin, hin = hin, tout = .3, optns = optns)
#' plot(res$qout[[1]][2,])
#' plot(res$dout[[1]][[2]]$x, res$dout[[1]][[2]]$y)
#' }
#' @references
#' \cite{Bhattacharee, S., & Müller, H.-G. (2022). "Concurrent object regression." Electronic Journal of Statistics16(2): 4031-4089.}
#' @export



partial_glob_den_CORE <- function(xin = NULL, tin = NULL, 
                                  qin = NULL, yin = NULL, hin = NULL,
                                  xout = NULL, tout = NULL, optns = list()){
  
  if (is.null(xin)){
    stop ("xin no default and must be input by users.")
  }
  if (is.null(tin)){
    stop ("tin no default and must be input by users.")
  }
  
  if(!is.list(xin)){
    stop('xin must be a list - xin is a list of length n (for n subjects) each element of the list is a matrix of size ni * p, 
         where ni is the number of observations for the ith subject and p is the dimension of each predictor')
  }
  if(!is.list(tin)){
    stop('tin must be a list - tin is a list of length n (for n subjects) each element of the list is a vector of lenght ni, 
         where ni is the number of observations for the ith subject.')
  }
  
  if(!is.null(qin)){
    if(!is.list(qin)){
      stop('qin must be a list - qin is a list of length n (for n subjects) each element of the list is a matrix of size ni * m, 
         where ni is the number of observations for the ith subject and m is the length of the grid over which the quantile function is observed')
    } else{
      if(is.null(optns$diff_qSup)){
        optns$diff_qSup = FALSE
      } else{
        if(optns$diff_qSup == TRUE){
          qin = lapply(1:length(qin), function(ind){
            ni = length(qin[[ind]])
            t(sapply(1:ni, function(l){
              if (!is.list(qin1[[ind]][[l]])) {
                stop ("If qin are entered with diff. support, then for each subject and 
                    at each time point where the quantiles are observed, a list wiht two components
                    x and y containing the support and values of the quantile function.")
              } else if (is.null(qin1[[ind]][[l]]$x) | is.null(qin1[[ind]][[l]]$y)) {
                stop ("If qin are entered with diff. support, then for each subject and 
                    at each time point where the quantiles are observed, a list wiht two components
                    x and y containing the support and values of the quantile function.")
              }
              qin_il <- approx(x = qin1[[ind]][[l]]$x, y = qin1[[ind]][[l]]$y, xout = qSup, rule = 2)$y
            }))
          })
        }
      }
    }
  }
  
  
  if(length(unique(
    sapply(1:length(xin), function(ind) ncol(xin[[ind]]) )
  )) !=1){
    stop('xin for all subjects must have the same dimension')
  }
  
  if(is.list(qin)){
    if(length(qin)!=length(xin)){
      stop('qin (as list) and xin must have the same length (same number of subjects)')
    } 
    if(sum(sapply(1:n, function(ind) nrow(xin[[ind]]) != nrow(qin[[ind]]) ))!=0){
      stop('qin and xin must be observed at same number of times for each subject')
    }
  }
  
  if (is.null(xout)){
    warning("xout not provided by the user- the model will be fitted at each point of xin")
  }
  
  if (is.null(tout)){
    warning("tout not provided by the user- the model will be fitted at each point of tin")
  }
  
  if(!is.null(xout)){
    if(!is.vector(xout)){
      stop('xout if enterted by the user must be a vector')
    }
    if(ncol(xin[[1]])!= length(xout)){
      stop(' xin and xout must have the same dimension')
    }
  }
  
  if(!is.null(tout)){
    if(!is.vector(tout)){
      stop('tout if enterted by the user must be a vector')
    }
  }
  
  if(any(sapply(1:n, function(ind) nrow(xin[[ind]]))<2)){
    stop('The data is too sparse')
  }
  
  #bwRange : vector of length 2
  if(!is.null(optns$bwRange)){
    if(!is.vector(optns$bwRange)){
      stop('bwRange must be a vector')
    }
    if(length(optns$bwRange)!=2){
      stop('bwRange must have the lower and upper bound for the bandwidth range')
    }
    if(sum(optns$bwReg<=0)>0){
      stop('bwReg must contain positive bandwidths')
    }
  }
  
  ####
  if (is.null(yin) & is.null(qin) & is.null(hin)){
    stop ("One of the three arguments, yin, hin and qin, should be input by users.")
  }
  
  if (!is.null(optns$qSup)) {
    if (min(optns$qSup) != 0 | max(optns$qSup) - 1 != 0)
      stop ("optns$qSup must have minimum 0 and maximum 1.")
    if (sum(duplicated(optns$qSup)) > 0) {
      optns$qSup <- unique(optns$qSup)
      warning ("optns$qSup has duplicated elements which has been removed.")
    }
    if (is.unsorted(optns$qSup)) {
      optns$qSup <- sort(optns$qSup)
      warning ("optns$qSup has been reordered to be increasing.")
    }
  } else {
    if (!(is.null(yin) & is.null(hin))) {
      if(is.null(optns$nqSup)) {
        optns$nqSup <- 101
      }
      optns$qSup <- seq(0,1,length.out = optns$nqSup)
    } else {
      if (is.list(qin)) {
        optns$qSup <- seq(0,1,length.out = ncol(qin[[1]]))
        warning ("optns$qSup is missing and is set by default as an equidistant grid on [0,1] with length equal to the number of columns in list qin.")
      } else {
        if(is.null(optns$nqSup)) {
          optns$nqSup <- 101
        }
        optns$qSup <- seq(0,1,length.out = optns$nqSup)
      }
    }
  }
  qSup <- optns$qSup
  
  optnsRegIdx <- match(c("bw_t","kernelReg","lower","upper","qSup","nqSup","bwRange"), names(optns))
  optnsRegIdx <- optnsRegIdx[!is.na(optnsRegIdx)]
  optnsReg <- optns[optnsRegIdx]
  if (is.null(optnsReg$kernelReg)){
    optnsReg$kernelReg <- "gauss"
  }
  names(optnsReg)[which(names(optnsReg) == "kernelReg")] <- "kernelReg"
  if (!is.null(optnsReg$bw_t)){
    names(optnsReg)[which(names(optnsReg) == "bw_t")] <- "bw_t"
  }
  
  optnsDen <- optns[-optnsRegIdx]
  if (!is.null(optnsDen$kernelDen))
    names(optnsDen)[which(names(optnsDen) == "kernelDen")] <- "kernel"
  if (!is.null(optnsDen$bwDen))
    names(optnsDen)[which(names(optnsDen) == "bwDen")] <- "userBwMu"
  
  #####
  if (!(is.null(yin) & is.null(hin))) {
    #require(fdadensity)
  
    if (!is.null(yin)) {
      if (!is.null(hin) | !is.null(qin))
        warning ("hin and qin are redundant when yin is available.")
      
      ##
      if (!is.list(yin))
        stop ("yin must be a list.")
      den <- lapply(1:length(yin), function(ind){
        ni = nrow(yin[[ind]])
        lapply(1:ni, function(l) CreateDensity(yin[[ind]][l,], optns = optnsDen) )
      })
    } else if (!is.null(hin)) {
      if (!is.null(qin))
        warning ("qin is redundant when hin is available.")
      if (!is.list(hin) | length(hin) != length(xin))
        stop ("hin must be a list of the same length as xin.")
      den = lapply(1:length(hin), function(ind){
        ni = length(hin[[ind]])
        lapply(1:ni, function(l){
          if (!is.list(hin[[ind]][[l]]))
            stop ("Each element of hin must be a list.")
          if (is.null(hin[[ind]][[l]]$breaks) & is.null(hin[[ind]][[l]]$mids))
            stop ("Each element of hin must be a list with at least one of the components breaks or mids.")
          if (is.null(hin[[ind]][[l]]$counts))
            stop ("Each element of hin must be a list with component counts.")
          CreateDensity(histogram = hin[[ind]][[l]], optns = optnsDen)
          })
        })
    }
    qin <- lapply(1:length(den), function(ind){
      ni = length(den[[ind]])
      t(sapply(1:ni, function(l){
        fdadensity::dens2quantile(dens = den[[ind]][[l]]$y, dSup = den[[ind]][[l]]$x, qSup = qSup)
      }))
    })
  }
  den <- lapply(1:length(qin), function(ind){
    ni = nrow(qin[[ind]])
    lapply(1:ni, function(l){
      frechet:::qf2pdf(qf = sort(qin[[ind]][l,]),prob = qSup)
    })
  })
  mean_qin_aux =  mean(sapply(1:length(qin), function(ind){
    ni = nrow(qin[[ind]])
    mean(sapply(1:ni, function(l){
      qin[[ind]][l,length(qin[[ind]][l,])] - qin[[ind]][l,1]
    }))
  }))
  if (is.null(optns$denLowerThreshold)) {
    optns$denLowerThreshold <- 0.001 *  mean_qin_aux
     
  } else if (optns$denLowerThreshold) {
    if(!is.numeric(optns$denLowerThreshold) | optns$denLowerThreshold < 0)
      optns$denLowerThreshold <- 0.001 *  mean_qin_aux
  }
  
  if (optns$denLowerThreshold) {
    if (sum(sapply(1:length(den), function(ind){
        ni = length(den[[ind]])
        sum(sapply(1:ni, function(l){
          sum(den[[ind]][[l]]$y < optns$denLowerThreshold/diff(range(den[[ind]][[l]]$x)))
        }))
      }))>0) {
    # density thresholding from below
    lapply(1:length(den), function(ind){
      ni = length(den[[ind]])
      lapply(1:ni, function(l){
        d = den[[ind]][[l]]
          lower <- optns$denLowerThreshold/diff(range(d$x))
          if (sum(d$y < lower) > 0) {
            d$y[d$y < lower] <- lower
            d$y <- d$y / pracma::trapz(d$x,d$y)
          }
          list(x=d$x, y=d$y)
          })
      })
      
    qin <- lapply(1:length(den), function(ind){
      ni = length(den[[ind]])
      t(sapply(1:ni, function(l){
        fdadensity::dens2quantile(dens = den[[ind]][[l]]$y, dSup = den[[ind]][[l]]$x, qSup = qSup)
      }))
    })
    }
  }
  ######
  
  #bwReg A vector of length p used as the bandwidth for the Fréchet regression or \code{"CV"} (default), i.e., a data-adaptive selection done by cross-validation.}
  if (!("bw_t"%in%names(optnsReg))) {
    optnsReg$bw_t <- "CV"
  }
  if (!is.numeric(optnsReg$bw_t)) {
    if (optnsReg$bw_t == "CV") {
      optnsReg$bw_t <- bwCV_pgm(xin, tin, qin, xout, tout, optns = optnsReg)
    } else {
      warning("optns$bw_t was mis-specified and is reset to be chosen by CV.")
      optnsReg$bw_t <- bwCV_pgm(xin, tin, qin, xout, tout, optns = optnsReg)
    }
  } else {
    if (optnsReg$bw_t < max(diff(sort(unlist(tin)))) & !is.null(optnsReg$ker)) {
      if(optnsReg$ker %in% c("rect","quar","epan")) {
          warning("optns$bw_t was set too small and is reset to be chosen by CV.")
          optnsReg$bw_t <- bwCV_pgm(xin, tin, qin, xout, tout, optns = optnsReg)
      }
    }
  }
  optns$bw_t <- optnsReg$bw_t
  qout <- partial_glob_CORE(xin, tin, qin, xout, tout, optns = optnsReg)
  
  if (!is.null(optnsDen$ndSup))
    names(optnsDen)[which(names(optnsDen) == "ndSup")] <- "nRegGrid"
  if (!is.null(optnsDen$dSup))
    names(optnsDen)[which(names(optnsDen) == "dSup")] <- "outputGrid"
  
###
  output_dout = function(qout){
    dout <- frechet:::qf2pdf(qout,prob = qSup, optns = optnsDen)
    dout[c("x","y")]
  }
  

  if (is.null(optnsDen$outputGrid)) {
    if(!is.null(tout) & !is.null(xout)){
      dout = output_dout(qout)
    } else {
      dout = lapply(1:length(qout), function(ind){
        ni = nrow(qout[[ind]])
        lapply(1:ni, function(l){
          output_dout(qout[[ind]][l,])
        })
      })
    }
    res <- list(xout = xout, dout = dout, qout = qout, qSup = qSup, xin=xin, tin = tin, din=den, qin=qin, optns=optns)
  } else {
    dSup <- optnsDen$outputGrid
    if(!is.null(tout) & !is.null(xout)){
      dout = qf2pdf(qout, prob = qSup, optns = optnsDen)$y
    } else {
      dout = lapply(1:length(qout), function(ind){
        ni = nrow(qout[[ind]])
        lapply(1:ni, function(l){
          qf2pdf(qout[[ind]][l,], prob = qSup, optns = optnsDen)$y
        })
      })
    }
    res <- list(xout = xout, dout = dout, dSup = dSup, qout = qout, qSup = qSup, xin=xin, tin = tin, din=den, qin=qin, optns=optns)
  }
  if(is.null(res$xout)){
    res$xout = xin
  }
  if(is.null(res$tout)){
    res$tout = tin
  }
  class(res) <- "denReg"
  return(res)
}
  