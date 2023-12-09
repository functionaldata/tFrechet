#' @title Generalized Fréchet integrals of 1D distribution 
#' @description Calculating generalized Fréchet integrals of 1D distribution (equipped with Wasserstein distance) 
#' @param phi An eigenfunction along which we want to project the distribution
#' @param t_out Support of \code{phi}
#' @param Q A \code{length(t_out)} X \code{length(Qout)} matrix whose jth row corresponds to the quantile function on grid \code{Qout} for the jth time point.   
#' @param Qout Support of the quantile valued process
#' @return A list of the following:
#' \item{f}{Quantile function corresponding to the frechet integral of \code{Q} along \code{phi}}
#' @examples 
#' #simulation as in the paper Dubey, P., & Müller, H. G. (2020).
#' # Functional models for time‐varying random objects. 
#' # JRSSB, 82(2), 275-327.
#' \donttest{
#' n <- 100
#' N <- 50
#' t_out <- seq(0,1,length.out = N)
#' 
#' phi1 <- function(t){
#'   (t^2-0.5)/0.3416
#' }
#' phi2 <- function(t){
#'   sqrt(3)*t
#' }
#' phi3 <- function(t){
#'   (t^3 - 0.3571*t^2 - 0.6*t + 0.1786)/0.0895
#' }
#' 
#' Z <- cbind(rnorm(n)*sqrt(12), rnorm(n), runif(n)*sqrt(72), runif(n)*sqrt(9))
#' mu_vec <- 1 + Z[,1] %*% t(phi1(t_out)) + Z[,2] %*% t(phi3(t_out))
#' sigma_vec <- 3 + Z[,3] %*% t(phi2(t_out)) + Z[,4] %*% t(phi3(t_out))
#' 
#' # grids of quantile function
#' Nq <- 40
#' eps <- 0.00001
#' Qout <- seq(0+eps,1-eps,length.out=Nq)
#' 
#' # I: four dimension array of n x n matrix of squared distances 
#' # between the time point u of the ith process and 
#' # process and the time point v of the jth object process, 
#' # e.g.: I[i,j,u,v] <- d_w^2(X_i(u) X_j(v)).
#' I <- array(0, dim = c(n,n,N,N))
#' for(i in 1:n){
#'   for(j in 1:n){
#'     for(u in 1:N){
#'       for(v in 1:N){
#'         #wasserstein distance between distribution X_i(u) and X_j(v) 
#'         I[i,j,u,v] <- (mu_vec[i,u] - mu_vec[j,v])^2 + (sigma_vec[i,u] - sigma_vec[j,v])^2
#'       }
#'     }
#'   }
#' }
#' 
#' # check ObjCov work 
#' Cov_result <- ObjCov(t_out, I, 3)
#' #Cov_result$lambda #12 6 1.75
#' 
#' # calculate Q 
#' i <- 6 # for the ith subject
#' Q <- t(sapply(1:N, function(t){
#'   qnorm(Qout, mean = mu_vec[i,t], sd = sigma_vec[i,t])
#' }))
#' 
#' score_result <- WassFIntegral(Cov_result$phi[,1], t_out, Q, Qout)
#' score_result$f
#' }
#' @references 
#' \cite{Dubey, P., & Müller, H. G. (2020). Functional models for time‐varying random objects. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(2), 275-327.}
#' @export
#' @import fdapace
#' @import quadprog
#' @import pracma

WassFIntegral <- function(phi, t_out, Q, Qout){

  m <- length(Qout)
  if(m<2){
    stop("length of Q out too small, Please increase ")
  }
  if(length(phi)!=length(t_out)){
    stop("length of t_out and phi inconsistent")
  }
  if(dim(Q)[2]!=length(Qout)){
    stop("column length of Qout and length Qout inconsistent")
  }
  if(dim(Q)[1]!=length(t_out)){
    stop("row length of Qout and length t_out inconsistent")
  }
  g_mini <- rep(0, m)
  phi_out <- phi / fdapace::trapzRcpp(t_out, phi)
  #g_mini = \int Q(S(t))(u)phi(t)dt
  for (i in 1:m){
    g_mini[i] <- fdapace::trapzRcpp(t_out, phi_out * Q[,i])
  }
  D <- pracma::Toeplitz(a=c(-1,rep(0,m-2)),
                       b=c(-1,1, rep(0,m-2)))
  f <- quadprog::solve.QP(diag(m),g_mini,t(D),rep(0,m-1))
  return(list(f=f$solution))
}




