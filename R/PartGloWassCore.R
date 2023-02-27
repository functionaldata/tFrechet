PartGloWassCore = function(xin, tin, qin, xout, tout, optns = list()){ #ker = ker.gauss, bw_t= NULL, lower=NULL, upper=NULL){
  
  if(!is.list(xin)){
    stop('xin must be a list')
  }
  if(!is.list(qin)){
    stop('qin must be a list')
  }
  if(!is.list(tin)){
    stop('tin must be a list')
  }
  
  if(length(unique(
    sapply(1:length(xin), function(ind) ncol(xin[[ind]]) )
  )) !=1){
    stop('xin for all subjects must have the same dimension')
  }
  
  if(length(qin)!=length(xin)){
    stop('qin (as list) and xin must have the same length (same number of subjects)')
  } 
  if(sum(sapply(1:length(xin), function(ind) nrow(xin[[ind]]) != nrow(qin[[ind]]) ))!=0){
    stop('qin and xin must be observed at same number of times for each subject')
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
      stop('tout if enterted by the user must be a number')
      }
  }
  
  if(is.null(optns$bw_t)){
    stop ("optns$bw_t has no default values and must be input by user.")
  }
  if(!is.numeric(optns$bw_t)){
    stop("optns$bw_t should be a number")
  }
  
  if(is.null(optns$kernelReg)){
    optns$kernelReg <- 'gauss'
  }
  ker <- frechet:::kerFctn(optns$kernelReg)
  
  
  n = length(xin)
  m = ncol(qin[[1]])
  p = ncol(xin[[1]])


##get the weights for the time-varying density regression    
  getLFRweights=function(x0, t0){
    Kt <- lapply(1:n, function(ind) ker((tin[[ind]] - t0)/ optns$bw_t))
    t <- lapply(1:n, function(ind) (tin[[ind]] - t0))
    
    mu00 = mean(sapply(1:n, function(ind) mean(Kt[[ind]])))
    mu01 = mean(sapply(1:n, function(ind) mean(Kt[[ind]]*t[[ind]])))
    mu02 = mean(sapply(1:n, function(ind) mean(Kt[[ind]]*t[[ind]]^2)))
    
    mean_x = sapply(1:p, function(j){
      Ly = lapply(1:n, function(ind) xin[[ind]][,j])
      Lt = lapply(1:n, function(ind)  sort(tin[[ind]]))
      bw_mu = suppressWarnings(fdapace:::CVLwls1D(Ly, Lt, kernel= "gauss", 
                                 npoly=1, nder=0, dataType= "sparse", kFolds = 5, 
                                 useBW1SE = FALSE))
      
      xin_vec = as.vector(unlist(Ly))[order(as.vector(unlist(Lt)))]
      tin_vec = sort(as.vector(unlist(Lt)))
      win = rep(1, length(tin_vec))
      mean_x_dim = fdapace::Lwls1D(bw_mu, kernel_type = "gauss", npoly = 1, nder = 0,
                                   xin = tin_vec, yin= xin_vec, xout = t0, win = win)
      return(mean_x_dim)
    })
    
    Sigma_20  = Reduce('+', lapply(1:n, function(ind){
      ni = nrow(xin[[ind]])
      Reduce('+', lapply(1:ni, function(l){
        Kt[[ind]][l] * t(xin[[ind]] - mean_x) %*% (xin[[ind]] - mean_x)
      }))/ni
    }))/n
    
    sigma0 = mu02*mu00 - mu01^2
    
    s_il = sapply(1:n, function(i){
      ni = nrow(xin[[i]])
      sapply(1:ni, function(l){
        Kt[[i]][l] * (
          ((xin[[i]][l,] - mean_x) %*% solve(Sigma_20) %*% (x0 - mean_x))
          + ((1/sigma0^2)*(mu02 - t[[i]][l]*mu01))
        )
      })
    })
    s = sum(unlist(s_il))
    s_il = sapply(1:n, function(i){
      ni = nrow(xin[[i]])
      sapply(1:ni, function(l) s_il[[i]][l]/s )
    })
    return(s_il)
  }
  A = cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  if(!is.null(optns$upper) & !is.null(optns$lower)){
    b0 = c(optns$lower, rep(0,m-1), -optns$upper)
  }else if(!is.null(optns$upper)){
    A = A[,-1]
    b0 = c(rep(0,m-1), -optns$upper)
  }else if(!is.null(optns$lower)){
    A = A[,-ncol(A)]
    b0 = c(lower,rep(0,m-1))
  }else{
    A = A[,-c(1,ncol(A))]
    b0 = rep(0,m-1)
  }
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")
  
  ##For any points xout (in dimesnion p), and tout (a number), compute the output quantile function
  compute_res = function(xout,tout,qin){
    n = length(qin)
    ss=getLFRweights(xout,tout)
    gx = sapply(1:n, function(i){
      ni = nrow(xin[[i]])
      # rowMeans(sapply(1:ni, function(l){
      #   qin[[i]][l,]* ss[[i]][l]
      # }))
      colMeans(qin[[i]] * ss[[i]]) *n*ni
    })
    gx = rowMeans(gx)
    res <- do.call(osqp::solve_osqp,
                   list(P=Pmat, q= -gx, A=Amat, l=b0, 
                        pars = osqp::osqpSettings(verbose = FALSE)))
    
    qout = sort(res$x)
    return(qout)
  }
  ##if either of tout or xout is NULL, compute the output quantile function at the inputted xin and/or tin  
  if(!is.null(tout) & !is.null(xout)){
    qout = compute_res(xout, tout, qin) 
  }
  
  if(is.null(xout) & !is.null(tout)){
    xout = xin
    qout = sapply(1:length(xout), function(ind){
      ni = nrow(xout[[ind]])
      qq = t(sapply(1:ni, function(l) compute_res(xout[[ind]][l,], tout, qin) ))
    })
  }
  if(is.null(tout) & !is.null(xout)){
    tout = tin
    qout = lapply(1:length(tout), function(ind){
      ni = length(tout[[ind]])
      qq = t(sapply(1:ni, function(l) compute_res(xout, tout[[ind]][l], qin) ))
    })
  }
  if(is.null(tout) & is.null(xout)){
    tout = tin; xout = xin
    qout = lapply(1:length(xout), function(ind){
      ni = nrow(xout[[ind]])
      qq = t(sapply(1:ni, function(l) compute_res(xout[[ind]][l,], tout[[ind]][l], qin) ))
    })
  }
  return(qout)
}

##For selecting the bandwidth in the t-direction by CV
SetBwRange <- function(xin, xout,kernel_type) {
  xinSt <- unique(sort(unlist(xin)))
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) -
                  xinSt[length(xin)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) *
                                                 ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin))/3
  if (bw.max < bw.min) {
    if (bw.min > bw.max*3/2) {
      #warning("Data is too sparse.")
      bw.max <- bw.min*1.01
    } else bw.max <- bw.max*3/2
  }
  return(list(min=bw.min, max = bw.max))
}

bwCV_pgm <- function(xin, tin, qin, xout, tout, optns = optns){
  tin_vec = sapply(1:n, function(ind) mean(tin[[ind]])) #rowMeans(tin)
  compareRange_t <- (tin_vec  > (min(tin_vec ) + diff(range(tin_vec))/5)) & 
    (tin_vec  < (max(tin_vec ) - diff(range(tin_vec ))/5))
  # k-fold
  objFctn <- function(bw_t) {
    optns1 <- optns
    optns1$bw_t <- bw_t
    folds_t <- numeric(length(tin))
    nn_t <- sum(compareRange_t)
    numFolds_t <- ifelse(nn_t > 30, 10, nn_t)
    
    tmp_t <- c(sapply(1:ceiling(nn_t/numFolds_t), function(i)
      sample(x = seq_len(numFolds_t), size = numFolds_t, replace = FALSE)))
    tmp_t <- tmp_t[1:nn_t]
    repIdx_t <- which(diff(tmp_t) == 0)
    for (i in which(diff(tmp_t) == 0)) {
      s <- tmp_t[i]
      tmp_t[i] <- tmp_t[i-1]
      tmp_t[i-1] <- s
    }
    
    folds_t[compareRange_t] <- tmp_t
    qfit <- lapply(seq_len(numFolds_t), function(foldidx_t){
      testidx_t <- which(folds_t == foldidx_t)
      res = lapply(testidx_t, function(k){
        ni = nrow(xin[[k]])
        rr = t(sapply(1:ni, function(j){
          PartGloWassCore(xin = xin[-k], tin = tin[-k], qin = qin[-k],
                            xout = xin[[k]][j,],
                            tout = tin[[k]][j],
                            optns = optns1)
                           # ker = ker.gauss, bw_t, lower=NULL, upper=NULL)
        }))
        return(rr)
      })
    })
    qobs =  lapply(seq_len(numFolds_t), function(foldidx_t){
      testidx_t <- which(folds_t == foldidx_t)
      return(qin[testidx_t])
    })
    
    
    err = mean(sapply(seq_len(numFolds_t), function(ff){
      testidt = length(qfit[[ff]])
      mean(sapply(1:testidt, function(k){
        mean(apply((qfit[[ff]][[k]] - qobs[[ff]][[k]])^2,1,pracma::trapz, x = optns1$qSup))
      }))
    }))
    return(err)
  }
  
  if(is.null(tout)){
    tout = tin
  }
  if(is.list(tout)){
    tout = unlist(tout)
  }
  aux=SetBwRange(xin = tin, xout = tout, kernel_type = optns$kernelReg)
  bwRange <- c(aux$min,aux$max)
  if(!is.null(optns$bwRange)){
    if (min(optns$bwRange) < min(bwRange)) {
      message("Minimum bandwidth is too small and has been reset.")
    }else{
      bwRange[1] <- min(optns$bwRange)
    }
    if (max(optns$bwRange) >  min(bwRange)) {
      bwRange[2] <- max(optns$bwRange)
    }else {
      message("Maximum bandwidth is too small and has been reset.")
    }
  }
  res <- optimize(f = objFctn, interval = bwRange)$minimum
  res
}
