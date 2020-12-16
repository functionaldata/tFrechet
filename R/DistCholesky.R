#' @title the Log-Cholesky and Cholesky distance between symmetric positive definite 
#' @noRd
#' @description the Log-Cholesky and Cholesky distance between two matrices
#' @param A an p by p matrix
#' @param B an p by p matrix
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{metric}{Metric type choice, "log_cholesky", "cholesky" - default: \code{log_cholesky} for log Cholesky metric}
#' }
#' @return A list containing the following fields:
#' \item{dist}{the distance between matrices \code{A} and \code{B}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples
#' p <- 3
#' a <- matrix(rnorm(p*p),p,p)
#' b <- matrix(rnorm(p*p),p,p)
#' A <- a %*% t(a)
#' B <- b %*% t(b)
#' res <- DistCholesky(A,B)

DistCholesky <- function(A=NULL,B=NULL,optns = list()){
  
  if(!is.matrix(A) | !is.matrix(B) ){
    stop('A and B must be of matrix class')
  }
  if(nrow(A)!=ncol(A) | nrow(B)!=ncol(B)){
    stop('Both A and B must be square matrices')
  }
  if(sum(dim(A)!=dim(B))>0){
    stop('Both A and B must have the same dimension')
  }

  
  if(is.null(optns$metric)){
    metric = 'log_cholesky'
  } else {
    metric =  optns$metric
  }
  
  p <- dim(A)[1]
  M.a <- chol(A)
  D.a <- diag(M.a)
  L.a <- M.a - diag(D.a)
  
  M.b <- chol(B)
  D.b <- diag(M.b)
  L.b <- M.b - diag(D.b)  
  
  # if( (D.a[p] <= 0) | (D.b[p] <= 0) | sum(A!=t(A))>0 |sum(B!=t(B))>0 ){
  #   stop('Both A and B must be symmetric positive definite matrices')
  # }
  # 
  if(metric == 'log_cholesky'){
    dist <- sqrt( sum((L.a-L.b)^2)+sum((log(D.a)-log(D.b))^2) ) 
  }
  if(metric == 'cholesky'){
    dist <- sqrt(sum((M.a-M.b)^2)) 
  }
  return(list(dist=dist,optns=list(metric=metric)))
}
