#'@title Plot of a single covariance matrix.
#'@noRd
#'@param mout A covariance or correlation matrix.
#' @param optns A list of optns control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control optns are
#' \describe{
#' \item{plot.type}{Character with two choices, "continuous" and "categorical".
#' The former plots the correlations in a continuous scale of colors by magnitude
#' while the latter categorizes the positive and negative entries into two different colors.
#' Default is "continuous"}
#' \item{plot.clust}{Character, the ordering method of the correlation matrix.
#' "original" for original order (default);
#' "AOE" for the angular order of the eigenvectors;
#' "FPC" for the first principal component order;
#' "hclust" for the hierarchical clustering order, drawing 4 rectangles on the graph according to the hierarchical cluster;
#' "alphabet" for the alphabetical order.}
#' \item{plot.method}{Character, the visualization method of correlation matrix to be used.
#' Currently, it supports seven methods, named "circle" (default), "square", "ellipse", "number", "pie", "shade" and "color". }
#' \item{CorrOut}{Logical, indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{plot.display}{Character, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.}
#'}
#'@examples
#'\donttest{
#'yy = matrix(c(rnorm(100)),nrow =10)
#'mm = cov(yy)
#'covplot(mm)
#'covplot(mm, optns= list(plot.type = "categorical", plot.clust = "hclust"))
#'covplot(mm, optns= list(plot.clust = "hclust"))
#'}
#'@importFrom corrplot corrplot


covplot <- function(mout,optns = list()){
  if(is.null(optns$plot.type)){
    plot.type = "continuous"
  }
  else{
    plot.type = optns$plot.type
  }
  if(is.null(optns$plot.clust)){
    plot.clust = NULL
    addrect = 4
  }
  else{
    plot.clust = optns$plot.clust
    if(plot.clust == "hclust"){
      addrect = 4
    }
  }
  if(is.null(optns$plot.method)){
    plot.method = "circle"
  }
  else{
    plot.method = optns$plot.method
  }
  if(is.null(optns$CorrOut)){
    CorrOut = FALSE
  }
  else{
    CorrOut = optns$CorrOut
  }
  if(is.null(optns$plot.display)){
    plot.display = "full"
  }
  else{
    plot.display = optns$plot.display
  }


  if(plot.type == "continuous"){
    col <- colorRampPalette(c("blue","white","red"))(200)
    if((dim(mout)[1] >10)){
      corrplot::corrplot(mout, method = plot.method, col = col, bg='lightblue',
                         type = plot.display, outline = TRUE,
                         order = plot.clust, addrect = addrect,
                         diag = FALSE, tl.pos = "n", is.corr = CorrOut)
    }
    else{
      corrplot::corrplot(mout, method = plot.method, col = col, bg='lightblue',
                         type = plot.display, outline = TRUE,
                         order = plot.clust, addrect = addrect,
                         diag = FALSE, is.corr = CorrOut)
    }
  }
  if(plot.type == "categorical"){
    col =c("blue","red")
    if((dim(mout)[1] >10)){
      corrplot::corrplot(mout, method = plot.method, col = col, bg='lightblue',
                         type = plot.display, outline = TRUE,
                         order = plot.clust, addrect = addrect,
                         diag = FALSE, tl.pos = "n",is.corr = CorrOut)
    }
    else{
      corrplot::corrplot(mout, method = plot.method, col = col, bg='lightblue',
                         type = plot.display, outline = TRUE,
                         order = plot.clust, addrect = addrect,
                         diag = FALSE, is.corr = CorrOut)
    }
  }
}
