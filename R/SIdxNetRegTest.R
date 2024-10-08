library(parallel)

SIdxNetTest <- function(est, b0, xin, Min, nboot = 500, bw, M, verbose = F){
  ## est: estimate from SIdxCovReg
  ## b0: true value or H0
  ## xin, Min: data input
  ## nboot: number of bootstrap replicate
  ## bw, M: parameters for CovBoot_est; should be given, no cv just for now
  ## verbose: print the iteration counts
  
  if (length(est) != length(b0)){
    
    stop("The length of est and b0 should be the same.")
    
  }
  
  if (abs(sum(b0^2) -1) > 1e-04 ){
    
    warning("The norm of b0 is not 1, we normalize it")
    b0 = b0 / sqrt(sum(b0^2))
    
  }
  
  NetBoot_est <- function(est, xin, Min, reps, bw, M, verbose = verbose){
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
        fit = SIdxNetReg(xin_resamp[[i]], Min_resamp[[i]], bw, M, verbose = verbose, iter = 100)
        fit$est
      })
    
    est_signed <- mclapply(b_est, function(x) {x[2:p] * sign(sum(est[2:p] * x[2:p]))})
    return(est_signed)
  }
  
  boot_res = NetBoot_est(est = est, xin = xin, Min = Min, reps = nboot, bw = bw, M = M, verbose = verbose)
  
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

set.seed(1)
n = 100
m=3

# Set Parameters
X = matrix(c(runif(n, -1, 1), 
             runif(n, -1, 1),
             runif(n, 1, 2),
             
             rgamma(n, 3, 1),
             rgamma(n, 4, 1),
             rgamma(n, 5, 1),
             
             rbinom(n, 1, 0.2),
             rbinom(n, 1, 0.3),
             rbinom(n, 1, 0.5)), n) 

y = lapply(1:n, function(i){
  a = 2 * sin(pi*X[i,1])^2*X[i,7] + cos(pi*X[i,2])^2*(1-X[i,7])
  b = X[i,4]*X[i,8]+X[i,5]*(1-X[i,8])
  
  Vec = -rbeta(m*(m-1)/2, shape1 = a, shape2 = b)
  temp = matrix(0, nrow = m, ncol = m)
  temp[lower.tri(temp)] = Vec
  temp <- temp + t(temp)
  diag(temp) = -colSums(temp)
  return(temp)
})

y_mat = array(0, c(m, m, n))
for(j in 1:n){
  y_mat[,,j] = y[[j]]
}

res_net = SIdxNetReg(xin = X, Min = y_mat)
b0 = c(1,0,0,0,0,0,0,0,0)
test_res = SIdxNetTest(res_net$est, b0, xin = X, Min = y_mat,
            nboot = 20, bw = res_net$bw, M = res_net$M)
test_res
