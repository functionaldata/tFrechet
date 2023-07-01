#rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#library(devtools)
#library(testthat)
#devtools::load_all("../../R/")
#devtools::document("../../R")
#devtools::load_all("./tFrechet/")
#devtools::document(pkg = "./tFrechet/")

#devtools::test_active_file("./tFrechet/tests/testthat/test_PartGloDenCore.R")


n =  30; m = 50; K = 2; dim = 2
e_val_fn = function(x,k){
  phi1 = function(x) -cos(pi*x/10)/sqrt(5) ##First eigenfunction from PACE paper. Domain is [0,10]
  phi_k = function(x, k) sin((2*k - 1)*pi*x/10)/sqrt(5)
  mu = function(x) x+sin(x); ##mean function
  lamb1 = seq(4,1,length.out = max((K-3),1))
  lamb = .7^(k-1)
  return(list(phi1 = phi1, phi_k = phi_k, mu = mu, lamb1 = lamb1, lamb = lamb))
}
generateData_K = function(n, K, dim){
  xi = matrix(0, nrow = n, ncol = K)
  X_tilde = list()
  T = list()
  lamb1 = e_val_fn(1,1)$lamb1
  lambda = c(lamb1, sapply((length(lamb1)+1) :K, function(j) e_val_fn(x,j)$lamb))
  Lambda = diag(lambda)
  for(i in 1:n){
    N = sample(2:20,1, replace = FALSE)
    T[[i]] = runif(N,0,1)
    X_tilde[[i]]  = replicate(dim,{
      L = rnorm(K,0,1)
      xi[i,] = sqrt(Lambda)%*%L  ##rnorm(dim(Lambda)[1]) Generate true scores
      phi = matrix(0, nrow = N, ncol = K)
      phi[,1] = sapply(T[[i]], function(x) e_val_fn(x, K)$phi1(x))
      for( j in 2:K){
        phi[,j] = sapply(T[[i]], function(x) e_val_fn(x, j)$phi_k(x, j))
      }
      mu = sapply(T[[i]], function(x) e_val_fn(x, K)$mu(x))
      X1 = (mu + phi %*%(xi[i,]))
    }, simplify = TRUE)
  }
  return(list(X_tilde = X_tilde , T = T))
}
##calculates true quantile model at given x and t values
true_reg = function(x,t){
  zeta_xt = mean(x) + .5*t^2
  nu_xt = .6 + .1*sin(10*pi*t)
  reg =  zeta_xt +  nu_xt * qnorm(seq(0.01,.99,length.out = m))
  reg = sort(reg,decreasing = FALSE)
  return(reg)
}
##function to generate quantiles
gen_quantile = function(x,t){
  mu_given = rnorm(n = 1, mean = mean(x) + .5*t^2 , .1)
  sig_given = rgamma(1,shape = ((.6 + .1*sin(10*pi*t))^2/.1), scale = (.1/(.6 + .1*sin(10*pi*t))))
  Q_Y = mu_given + sig_given*qnorm(seq(0.01,.99,length.out = m))
  dens_Y = frechet:::qf2pdf(Q_Y)
  return(Q_Y)
}
data = generateData_K(n,K, dim)
gen_pred = data$X_tilde
xin = t(sapply(1:n, function(i){
  gen_pred[[i]]
}))
tin = data$T
qSup <- seq(1e-6,1-1e-6,length.out = (m-2))
qin = lapply(1:n, function(i){
  ni = nrow(gen_pred[[i]])
  t(sapply(1:ni, function(j){
    gen_quantile(gen_pred[[i]][j,], tin[[i]][j])
  }))
})
nobs = 100
gen_yin = function(x,t){
  mu_given = rnorm(n = 1, mean = x + .5*t^2 , .1)
  sig_given = rgamma(1,shape = ((.6 + .1*sin(10*pi*t))^2/.1), scale = (.1/(.6 + .1*sin(10*pi*t))))
  Y = rnorm(nobs,mu_given ,sig_given)
  return(Y)
}
yin = lapply(1:n, function(i){
  ni = nrow(xin[[i]])
  t(sapply(1:ni, function(j){
    gen_yin(xin[[i]][j,], tin[[i]][j])
  }))
})
hin = lapply(1:n, function(i){
  ni = nrow(yin[[i]])
  lapply(1:ni, function(j){
    hist(yin[[i]][j,], breaks = 50)
  })
})
qtrue = lapply(1:n, function(i){
  ni = length(tin[[i]])
  t(sapply(1:ni, function(j){
    true_reg(xin[[i]][j,], tin[[i]][j])
  }))
})
mean_qin = mean(sapply(1:n, function(i){
  ni = nrow(qin[[i]])
  mean(sapply(1:ni, function(j) mean(qin[[i]][j,])))
}))
#########################################
test_that("Works with fully observed distributions", {
  xout = tout = NULL
  res = suppressWarnings(PartGloDenCore(xin = xin, tin = tin, qin = qin, xout = xout, tout = tout,
                       optns = list(qSup = c(0,qSup,1))))
  qout_compare = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
     res$qout[[i]][j,] #[j,-c(1,length(res$qSup))]
    }))
  })
  expect_true((mean(sapply(1:n, function(i){
    ni = nrow(qtrue[[i]])
    mean(sqrt(sapply(1:ni, function(j){
      (qtrue[[i]][j,] - qout_compare[[i]][j,])^2 * diff(qSup)[2] 
    })))
  })))/mean_qin <.1)
})


#########################################
test_that("Works with discrete noisy measurements", {
  xout <- NULL; tout = NULL
  res = suppressWarnings(PartGloDenCore(xin = xin, tin = tin, yin = yin, xout = xout, tout = tout,
                                        optns = list(qSup = c(0,qSup,1))))
  qout_compare = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
      res$qout[[i]][j,] #[j,-c(1,length(res$qSup))]
    }))
  })
  expect_true((mean(sapply(1:n, function(i){
    ni = nrow(qtrue[[i]])
    mean(sqrt(sapply(1:ni, function(j){
      (qtrue[[i]][j,] - qout_compare[[i]][j,])^2 * diff(qSup)[2]
      #pracma::trapz(x = c(0,qSup,1), y = (qtrue[[i]][j,] - qout_compare[[i]][j,])^2)
    })))
  })))/mean_qin <.1)
})
#########################################
test_that("Works with specifying outputGrid", {
  xout <- NULL; tout = NULL
  dSup <- seq(-1.5,1.5,0.1)
  res = suppressWarnings(PartGloDenCore(xin = xin, tin = tin, qin = qin, xout = xout, tout = tout,
                                        optns = list(qSup = c(0,qSup,1), outputGrid = dSup)))
  expect_true("dSup" %in% names(res))
  expect_equal(sum(abs(res$dSup - dSup)), 0)
})
#########################################
test_that("Generates warnings when more than one of the three, yin, hin, and qin, are specified and priority order is yin, hin, qin", {
  xout <- NULL; tout = NULL
  expect_warning(res_yh <- PartGloDenCore(xin = xin, tin = tin, yin = yin, hin = hin,
                                          xout = xout, tout = tout, optns = list(qSup = c(0,qSup,1))))
  expect_warning(res_yq <- PartGloDenCore(xin = xin, tin = tin, yin = yin, qin = qin,
                                          xout = xout, tout = tout, optns = list(qSup = c(0,qSup,1))))
  expect_warning(res_hq <- PartGloDenCore(xin = xin, tin = tin, qin = qin, hin = hin,
                                          xout = xout, tout = tout, optns = list(qSup = c(0,qSup,1))))
  
  res_y <- suppressWarnings(PartGloDenCore(xin = xin, tin = tin, yin = yin,
                                           xout = xout, tout = tout, optns = list(qSup = c(0,qSup,1))))
  res_h <- suppressWarnings(PartGloDenCore(xin = xin, tin = tin, hin = hin,
                                           xout = xout, tout = tout, optns = list(qSup = c(0,qSup,1))))
  expect_equal(
    Reduce('+', lapply(1:n, function(ind){
      ni = nrow(qin[[ind]])
      sum(sapply(1:ni, function(j){
        abs(res_y$qin[[ind]][j,] - res_yq$qin[[ind]][j,])
      }))
    })), 0)
  
  expect_equal(
    Reduce('+', lapply(1:n, function(ind){
      ni = nrow(qin[[ind]])
      sum(sapply(1:ni, function(j){
        abs(res_h$qin[[ind]][j,] - res_hq$qin[[ind]][j,])
      }))
    })), 0)
})
#########################################
test_that("Works with fully observed distributions when only providing bandwidth range; bandwidth selected by CV", {
  xout <- NULL; tout = NULL
  res <- suppressWarnings(PartGloDenCore(xin = xin, tin = tin, qin = qin,
                                         xout = xout, tout = tout,
                                         optns = list(qSup = c(0,qSup,1), bwRange = c(0.2,0.7))))
  
  qout_compare = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
      res$qout[[i]][j,] #[j,-c(1,length(res$qSup))]
    }))
  })
  expect_true((mean(sapply(1:n, function(i){
    ni = nrow(qtrue[[i]])
    mean(sqrt(sapply(1:ni, function(j){
      (qtrue[[i]][j,] - qout_compare[[i]][j,])^2 * diff(qSup)[2] 
      #pracma::trapz(x = c(0,qSup,1), y = (qtrue[[i]][j,] - qout_compare[[i]][j,])^2)
    })))
  })))/mean_qin <.1)
})
###################################
test_that("bw found by CV is in range when providing fully observed distributions and bandwidth range", {
  xout <- NULL; tout = NULL
  bwRange = c(0.2,0.7)
  res <- suppressWarnings(PartGloDenCore(xin = xin, tin = tin, qin = qin,
                                         xout = xout, tout = tout,
                                         optns = list(qSup = c(0,qSup,1), bwRange = bwRange)))
  
  expect_true((res$optns$bw_t >= bwRange[1])&&(res$optns$bw_t<=bwRange[2]))
})
#########################################
test_that("Works with fully observed distributions for higher dimension of xin (p > 2)", {
  dim = 4
  data = generateData_K(n,K, dim)
  xin = data$X_tilde
  tin = data$T
  qSup <- seq(1e-6,1-1e-6,length.out = (m-2))
  qin = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
      gen_quantile(xin[[i]][j,], tin[[i]][j])
    }))
  })
  
  qtrue = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
      true_reg(xin[[i]][j,], tin[[i]][j])
    }))
  })
  mean_qin = mean(sapply(1:n, function(i){
    ni = nrow(qin[[i]])
    mean(sapply(1:ni, function(j) mean(qin[[i]][j,])))
  }))
  xout = tout = NULL
  res = suppressWarnings(PartGloDenCore(xin = xin, tin = tin, qin = qin, xout = xout, tout = tout,
                                        optns = list(qSup = c(0,qSup,1))))
  qout_compare = lapply(1:n, function(i){
    ni = length(tin[[i]])
    t(sapply(1:ni, function(j){
      res$qout[[i]][j,] 
    }))
  })
  expect_true((mean(sapply(1:n, function(i){
    ni = nrow(qtrue[[i]])
    mean(sqrt(sapply(1:ni, function(j){
      (qtrue[[i]][j,] - qout_compare[[i]][j,])^2 * diff(qSup)[2]
    })))
  })))/mean_qin <.1)
})

