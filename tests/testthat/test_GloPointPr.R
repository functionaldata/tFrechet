library(testthat)


test_that("Works with p=1 and no options", {
  set.seed(1)
  n=100
  alpha_n=sqrt(n)
  alpha1=2.0
  beta1=1.0
  gridQ=seq(0,1,length.out=500+2)[2:(500+1)]
  X=runif(n,0,1)#p=1
  tau=matrix(0,nrow=n,ncol=1)
  for(i in 1:n){
    tau[i]=alpha1+beta1*X[i]+truncnorm::rtruncnorm(1, a=-0.3, b=0.3, mean = 0, sd = 1.0)
  }
  Ni_n=matrix(0,nrow=n,ncol=1)
  u0=0.4
  u1=0.5
  u2=0.05
  u3=-0.01
  tin=list()
  for(i in 1:n){
    Ni_n[i]=rpois(1,alpha_n*tau[i])
    mu_x=u0+u1*X[i]+truncnorm::rtruncnorm(1,a=-0.1,b=0.1,mean=0,sd=1)
    sd_x=u2+u3*X[i]+truncnorm::rtruncnorm(1,a=-0.02,b=0.02,mean=0,sd=0.5)
    if(Ni_n[i]==0){
      tin[[i]]=c()
    }else{
      tin[[i]]=truncnorm::rtruncnorm(Ni_n[i],a=0,b=1,mean=mu_x,sd=sd_x) #Sample from truncated normal
    }
  }
  res=GloPointPrReg(xin=matrix(X,ncol=1),tin=tin,T0=1,xout=matrix(seq(0,1,length.out=10),ncol=1))
  
  expect_true(all(dim(res$intensityReg)==c(length(res$dSup),10)))
})

test_that("Works with p=2 and no options", {
  set.seed(1)
  n=100
  alpha_n=sqrt(n)
  alpha1=2.0
  beta1=1.0
  gridQ=seq(0,1,length.out=500+2)[2:(500+1)]
  X=cbind(runif(n,0,1),runif(n,0,1))#p=2
  tau=matrix(0,nrow=n,ncol=1)
  for(i in 1:n){
    tau[i]=alpha1+beta1*X[i,1]+0.5*beta1*X[i,2]+truncnorm::rtruncnorm(1, a=-0.3, b=0.3, mean = 0, sd = 1.0)
  }
  Ni_n=matrix(0,nrow=n,ncol=1)
  u0=0.4
  u1=0.5
  u2=0.05
  u3=-0.01
  tin=list()
  for(i in 1:n){
    Ni_n[i]=rpois(1,alpha_n*tau[i])
    mu_x=u0+u1*X[i,1]+0.1*X[i,2]+truncnorm::rtruncnorm(1,a=-0.1,b=0.1,mean=0,sd=1)
    sd_x=u2+u3*X[i,1]+0.01*X[i,2]+truncnorm::rtruncnorm(1,a=-0.02,b=0.02,mean=0,sd=0.5)
    if(Ni_n[i]==0){
      tin[[i]]=c()
    }else{
      tin[[i]]=truncnorm::rtruncnorm(Ni_n[i],a=0,b=1,mean=mu_x,sd=sd_x) #Sample from truncated normal
    }
  }
  res=GloPointPrReg(xin=matrix(X,ncol=2),tin=tin,T0=1,xout=matrix(cbind(seq(0,1,length.out=10),0),ncol=2))
  
  expect_true(all(dim(res$intensityReg)==c(length(res$dSup),10)))
})



test_that("Works with specifying dSup and p=2", {
  set.seed(1)
  n=100
  alpha_n=sqrt(n)
  alpha1=2.0
  beta1=1.0
  gridQ=seq(0,1,length.out=500+2)[2:(500+1)]
  X=cbind(runif(n,0,1),runif(n,0,1))#p=2
  tau=matrix(0,nrow=n,ncol=1)
  for(i in 1:n){
    tau[i]=alpha1+beta1*X[i,1]+0.5*beta1*X[i,2]+truncnorm::rtruncnorm(1, a=-0.3, b=0.3, mean = 0, sd = 1.0)
  }
  Ni_n=matrix(0,nrow=n,ncol=1)
  u0=0.4
  u1=0.5
  u2=0.05
  u3=-0.01
  tin=list()
  for(i in 1:n){
    Ni_n[i]=rpois(1,alpha_n*tau[i])
    mu_x=u0+u1*X[i,1]+0.1*X[i,2]+truncnorm::rtruncnorm(1,a=-0.1,b=0.1,mean=0,sd=1)
    sd_x=u2+u3*X[i,1]+0.01*X[i,2]+truncnorm::rtruncnorm(1,a=-0.02,b=0.02,mean=0,sd=0.5)
    if(Ni_n[i]==0){
      tin[[i]]=c()
    }else{
      tin[[i]]=truncnorm::rtruncnorm(Ni_n[i],a=0,b=1,mean=mu_x,sd=sd_x) #Sample from truncated normal
    }
  }
  res=GloPointPrReg(xin=matrix(X,ncol=2),tin=tin,T0=1,xout=matrix(cbind(seq(0,1,length.out=10),0),ncol=2),optns=list(dSup=seq(0,1,length.out=100)))
  
  expect_equal((nrow(res$intensityReg)-100)^2+(length(res$dSup)-100)^2, 0)
})
