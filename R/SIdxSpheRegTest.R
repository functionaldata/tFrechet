library(parallel)

SIdxSpheTest <- function(est, b0, xin, yin, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxDenReg
  ## b0: true value or H0
  ## xin, yin: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for DenBoot_est; should be given, no cv just for now,
  #### Can use output from SIdxDenReg
  ## verbose: print the iteration counts?
  
  SpheBoot_est <- function(est, xin, yin, reps, bw, M){
    n = nrow(xin)
    p = length(est)
    samp_ind = mclapply(1:reps, function(i) {
      sample(1:n,n,replace=T)
    })
    
    xin_resamp = mclapply(1:reps, function(i){
      xin[samp_ind[[i]],]
    })
    
    yin_resamp = mclapply(1:reps, function(i){
      yin[samp_ind[[i]],]
    })  
    
    b_est <- mclapply(
      1:reps, 
      FUN = function(i){
        fit = SIdxSpheReg(xin_resamp[[i]], yin_resamp[[i]], bw, M, verbose = verbose)
        fit$est
      })
    
    est_signed <- mclapply(b_est, function(x) {x[2:p] * sign(sum(est[2:p] * x[2:p]))})
    return(est_signed)
  }
  tt <- SpheBoot_est(est = est, xin = xin, yin = yin, 
                    reps = nboot, bw = bw, M = M)
  cova_boot  = var(do.call(rbind, tt)) # Covariance matrix bootstrap estimate
  p = length(est)
  test_stat = c(M * t(est[2:p] - b0[2:p]) %*% solve(cova_boot, est[2:p] - b0[2:p]))
  test_stat_boot = sapply(tt, function(x){
    M * t(x - est[2:p]) %*% solve(cova_boot, x - est[2:p])
  })
  p_val = mean(test_stat_boot >= test_stat)
  res = list(cov_boot = cova_boot, test_stat = test_stat, 
             p_boot = p_val, p_asymp = pchisq(test_stat, df = p-1), df = p-1)
  return(res)
}

#### Test

set.seed(10101)
b <- c(4, 1.3, -2.5, 1.7)
b0 <- normalize(b) #0.7722353  0.2509765 -0.4826471  0.3282000
dat_sphe <- SpheGenerate_data(100, 0, b0, function(x) x)
fit_sphe <- SIdxSpheReg(xin = dat_sphe$xin, yin = dat_sphe$yin, iter = 500)

SIdxSpheTest(fit_sphe$est, b0, xin = dat_sphe$xin, yin = dat_sphe$yin,
            nboot = 50, bw = 0.2, M = 4)
