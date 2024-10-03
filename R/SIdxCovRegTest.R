library(parallel)

SIdxCovTest <- function(est, b0, xin, Min, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxCovReg
  ## b0: true value or H0
  ## xin, Min: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for CovBoot_est; should be given, no cv just for now
  ## verbose: print the iteration counts
  
  if (abs(sum(b0^2) -1) > 1e-04 ){
    
    warning("The norm of b0 is not 1, we normalize it")
    b0 = b0 / sqrt(sum(b0^2))
    
  }
  
  CovBoot_est <- function(est, xin, Min, reps, bw, M, verbose = verbose){
    n = nrow(xin)
    p = length(est)
    samp_ind = mclapply(1:reps, function(i) {
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
  
  boot_res = CovBoot_est(est = est, xin = xin, Min = Min, reps = nboot, bw = bw, M = M, verbose = verbose)
  
  p = length(est)
  
  cov_mat = matrix(0, nrow = p-1, ncol = p-1)
  for(i in 1:length(boot_res)){
    
    cov_mat = cov_mat + (boot_res[[i]] - est[2:p]) %*% t(boot_res[[i]] - est[2:p])
    
  }
  cova_boot = cov_mat / length(boot_res)
  
  test_stat = c(t(est[2:p] - b0[2:p]) %*% solve(cova_boot, est[2:p] - b0[2:p]))
  
  p_val = 1 - pchisq(test_stat, df = p-1)
  
  res = list(cov_boot = (M * cova_boot), test_stat = test_stat, pval_chisq = p_val, df = p-1)
  return(res)
}

#### Test
b <- c(3, -1.3, -3, 1.7)
b0 <- normalize(b)
b0 #0.6313342 -0.2735781 -0.6313342  0.3577560

set.seed(99)
dat <- CovGen_data_setting(500, b0, function(x) x)
res_cov <- SIdxCovReg(dat$xin, dat$Min, iter = 500)

SIdxCovTest(res_cov$est, b0, xin = dat$xin, Min = dat$Min,
            nboot = 50, bw = res_cov$bw, M = res_cov$M)

