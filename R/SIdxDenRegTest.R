library(parallel)

SIdxDenTest <- function(est, b0, xin, qin, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxDenReg
  ## b0: true value or H0
  ## xin, qin: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for DenBoot_est; should be given, no cv just for now
  ## verbose: print the iteration counts?
  
  DenBoot_est <- function(est, xin, qin, reps, bw, M){
    n = nrow(xin)
    p = length(est)
    samp_ind = mclapply(1:reps, function(i) {
      #set.seed(i)
      sample(1:n,n,replace=T)
    })
    
    xin_resamp = mclapply(1:reps, function(i){
      xin[samp_ind[[i]],]
    })
    
    qin_resamp = mclapply(1:reps, function(i){
      qin[samp_ind[[i]],]
    })  
    
    b_est <- mclapply(
      1:reps, 
      FUN = function(i){
        fit = SIdxDenReg(xin_resamp[[i]], qin_resamp[[i]], bw, M, verbose = verbose)
        fit$est
      })
    
    est_signed <- mclapply(b_est, function(x) {x[2:p] * sign(sum(est[2:p] * x[2:p]))})
    return(est_signed)
  }
  tt <- DenBoot_est(est = est, xin = xin, qin = qin, 
                    reps = nboot, bw = bw, M = M)
  cova_boot  = var(do.call(rbind, tt)) # Covariance matrix bootstrap estimate
  p = length(est)
  test_stat = c(M * t(est[2:p] - b0[2:p]) %*% solve(cova_boot, est[2:p] - b0[2:p]))
  test_stat_boot = sapply(tt, function(x){
    M * t(x - est[2:p]) %*% solve(cova_boot, x - est[2:p])
  })
  p_val = mean(test_stat_boot >= test_stat)
  res = list(cov_boot = cova_boot, test_stat = test_stat, 
       pval_bootstrap = p_val, pval_chisq = 1 - pchisq(test_stat, df = p-1), df = p-1)
  return(res)
}

#### Test

set.seed(10101)
b <- c(4, 1.3, -2.5, 1.7)
b0 <- normalize(b)
dat_den <- DenGen_data_setting(100, b0, function(x) x)
fit_den <- SIdxDenReg(xin = dat_den$xin, qin = dat_den$qin, iter = 500)

SIdxDenTest(fit_den$est, b0, xin = dat_den$xin, qin = dat_den$qin,
            nboot = 50, bw = 0.8, M = 4)
