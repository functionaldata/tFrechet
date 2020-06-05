#' @title Converting a quantile function to a histogram
#' @noRd
#' @param qf a numerical vector holding the values of a quantile function at a probability grid.
#' @param prob a numerical vector holding the probability grid at which the quantile function takes values. By default, \code{prob} is an equidistant sequence with the same length as \code{qf}.
#' @param breaks a numerical vector holding the breakpoints between histogram bins. By default, \code{breaks} is an equidistant sequence from the minimum to the maximum values of \code{qf}.

qf2hist <- function(qf=NULL, prob=NULL, breaks=NULL){#, tol=1e-2){
  if(is.null(qf))
    stop("qf is missing.")
  if(!is.vector(qf))
    stop("qf should be a vector.")
  if (!is.numeric(qf))
    stop("qf should be a numerical vector.")
  if (is.unsorted(qf))
    stop("qf should be an increasingly sorted numerical vector.")
  if(is.null(prob))
    prob = seq(0,1,length.out=length(qf))
  if(length(prob)!=length(qf))
    stop("The length of prob should be the same as qf.")
  if(!is.vector(prob) | !is.numeric(prob) | is.unsorted(prob))
    stop("prob should be an increasingly sorted numerical vector.")

  #if(min(prob)>tol | max(prob) < 1-tol)
  #  stop("prob should be a vector with minimum 0 and maximum 1.")
  if(is.null(breaks))
    breaks = seq(min(qf), max(qf), length.out = 1e3)
  if(!is.vector(breaks) | !is.numeric(breaks) | is.unsorted(breaks) )
    stop("breaks should be an increasingly sorted numerical vector.")
  if(min(breaks) > min(qf) | max(breaks) < max(qf))
    stop("The range of breaks should cover that of qf.")

  cdf = approx(x=qf, y=prob, xout=breaks, ties = mean)$y
  if (min(qf) > min(breaks)) cdf[breaks < min(qf)] <- 0
  if (max(breaks) > max(qf)) cdf[breaks > max(qf)] <- 1
  if (sum(cdf>1)) cdf[cdf>1] <- 1
  density = (cdf[-1] - cdf[-length(cdf)])
  counts = as.integer(density * 1e5)
  density = density / (breaks[-1] - breaks[-length(breaks)])
  mids = (breaks[-1] + breaks[-length(breaks)]) / 2
  hist = list(breaks=breaks, counts=counts, density=density, mids=mids)
  return(hist)
}
