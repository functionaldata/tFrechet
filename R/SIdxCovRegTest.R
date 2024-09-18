library(parallel)

SIdxCovTest <- function(est, b0, xin, Min, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxCovReg
  ## b0: true value or H0
  ## xin, Min: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for CovBoot_est; should be given, no cv just for now
  ## verbose: print the iteration counts?
  CovBoot_est <- function(est, xin, Min, reps, bw, M){
    n = nrow(xin)
    p = length(est)
    samp_ind = mclapply(1:reps, function(i) {
      #set.seed(i)
      sample(1:n,n,replace=T)
    })
    
    xin_resamp = mclapply(1:reps, function(i){
      xin[samp_ind[[i]],]
    })
    
    Min_resamp = mclapply(1:reps, function(i){
      Min[,,samp_ind[[i]]]
    })  
    
    b_est <- mclapply(
      1:reps, 
      FUN = function(i){
        fit = SIdxCovReg(xin_resamp[[i]], Min_resamp[[i]], bw, M, verbose = verbose)
        fit$est
      })
    
    est_signed <- mclapply(b_est, function(x) {x[2:p] * sign(sum(est[2:p] * x[2:p]))})
    return(est_signed)
  }
  tt <- CovBoot_est(est = est, xin = xin, Min = Min, 
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
set.seed(1000)
b <- c(3, -1.3, -3, 1.7)
b0 <- normalize(b)
b0 #0.6313342 -0.2735781 -0.6313342  0.3577560

dat_cov <- CovGen_data_setting(100, b0, function(x) x)
fit_cov <- SIdxCovReg(xin = dat_cov$xin, Min = dat_cov$Min, bw= 0.25, M = 4, iter = 1000, verbose = F)
fit_cov$est
SIdxCovTest(fit_cov$est, b0, xin = dat_cov$xin, Min = dat_cov$Min,
            nboot = 500, bw = fit_cov$bw, M = fit_cov$M)

