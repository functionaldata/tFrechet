# L2 norm
l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}