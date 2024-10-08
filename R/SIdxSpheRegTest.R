library(parallel)

SIdxSpheTest <- function(est, b0, xin, yin, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxDenReg
  ## b0: true value or H0
  ## xin, yin: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for DenBoot_est; should be given, no cv just for now,
  #### Can use output from SIdxDenReg
  ## verbose: print the iteration counts
  
  if (length(est) != length(b0)){
    
    stop("The length of est and b0 should be the same.")
    
  }
  
  if (abs(sum(b0^2) -1) > 1e-04 ){
    
    warning("The norm of b0 is not 1, we normalize it")
    b0 = b0 / sqrt(sum(b0^2))
    
  }
  
  SpheBoot_est <- function(est, xin, yin, reps, bw, M, verbose = verbose){
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
  
  
  boot_res = SpheBoot_est(est = est, xin = xin,  yin = yin, reps = nboot, bw =bw, M = M, verbose = verbose)
  
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
set.seed(100)
b <- c(3, -1.3, -3, 1.7)
b0 <- normalize(b)
b0 #0.6313342 -0.2735781 -0.6313342  0.3577560

dat <- SpheGenerate_data(100, 0, b0, function(x) x)
res_sphe <- SIdxSpheReg(xin = dat$xin, yin = dat$yin)

test_res = SIdxSpheTest(res_sphe$est, b0, xin = dat$xin, yin = dat$yin,
                       nboot = 50, bw = res_sphe$bw, M = res_sphe$M)
test_res
