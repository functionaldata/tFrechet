#' @title Plots for Fréchet regression for univariate densities.
#' @param x A \code{denReg} object, result of \code{\link{DenFMean}}, \code{\link{GloDenReg}} or \code{\link{LocDenReg}}.
#' @param obj An integer indicating which output to be plotted; 1, 2, 3, 4, and 5 for \code{dout}, \code{qout}, \code{din}, \code{qin}, and reference chart for \code{qout}, respectively - default: 1.
#' @param prob A vector specifying the probability levels for reference chart if \code{obj} is set to 5. Default: \code{c(0.05,0.25,0.5,0.75,0.95)}.
#' @param xlab Character holding the label for x-axis; default: \code{"Probability"} when \code{obj} is 2 or 4, \code{""} when \code{obj} is 1 or 3, \code{"x"} when \code{obj} is 5.
#' @param ylab Character holding the label for y-axis; default: \code{"Quantile"} when \code{obj} is 2, 4, or 5, and \code{"Density"} when \code{obj} is 1 or 3.
#' @param main Character holding the plot title; default: \code{NULL}.
#' @param xlim A numeric vector of length 2 holding the range of the x-axis to be drawn; default: automatically determined by the input \code{x}.
#' @param ylim A numeric vector of length 2 holding the range of the y-axis to be drawn; default: automatically determined by the input \code{x}.
#' @param col.bar A logical variable indicating whether a color bar is presented on the right of the plot - default: \code{TRUE}.
#' @param widrt A scalar giving the width ratio between the main plot and the color bar - default: 4.
#' @param col.lab A character giving the color bar label.
#' @param nticks An integer giving the number of ticks used in the axis of color bar.
#' @param ticks A numeric vector giving the locations of ticks used in the axis of color bar; it overrides \code{nticks}.
#' @param add Logical; only works when \code{obj} is 5. If \code{TRUE} add to an already existing plot. Taken as \code{FALSE} (with a warning if a different value is supplied) if no graphics device is open.
#' @param pos.prob \code{FALSE} or a scalar less than 0 or larger than 1. FALSE: no probability levels will be labeled on the quantile curves; a scalar between 0 and 1: indicating where to put the probability levels along the curves on growth charts: 0 and 1 correspond to left and right ends, respectively. Default: 0.9.
#' @param colPalette A function that takes an integer argument (the required number of colors) and returns a character vector of colors interpolating the given sequence
#' (e.g., \code{\link{heat.colors}}, \code{\link{terrain.colors}} and functions created by \code{\link{colorRampPalette}}).
#' Default is \code{colorRampPalette(colors = c("pink","royalblue"))} for more than one curves and \code{"black"} otherwise.
#' @param ... Can set up \code{lty}, \code{lwd}, etc.
#' @return No return value.
#' @note see \code{\link{DenFMean}}, \code{\link{GloDenReg}} and \code{\link{LocDenReg}} for example code.
#' @export

plot.denReg <- function(x, obj = NULL, prob = NULL,
                        xlab = NULL, ylab = NULL, main = NULL,
                        ylim = NULL, xlim = NULL,
                        col.bar = TRUE, widrt = 4, col.lab = NULL,
                        nticks = 5, ticks = NULL,
                        add = FALSE, pos.prob = 0.9,
                        colPalette = NULL,...) {
  if (is.null(obj)) obj <- 1
  if(! obj %in% 1:5) stop("obj is mis-specified.")

  if (is.null(colPalette))
    colPalette <- colorRampPalette(colors = c("pink","royalblue"))
  if (obj == 5) {
    # growth chart
    if (is.null(prob))
      prob <- c(0.05,0.25,0.5,0.75,0.95)
    plotObj <- sapply(apply(x$qout, 1, approx, x=x$qSup, xout=prob), with, y)
    plotGrid <- x$xout
    if (!is.vector(plotGrid)) {
      plotGrid <- plotGrid[,1]
      warning("x$xout is not a vector. Only the first column is used.")
    }
    if (is.null(ylim)) ylim <- range(plotObj)
    if (is.null(xlim)) xlim <- range(plotGrid)
    if(is.null(ylab)) ylab <- "Quantile"
    if(is.null(xlab)) xlab <- "x"

    n <- nrow(plotObj)
    if (n == 1)
      colPalette <- function(num) rep("black",num)
    i <- 1
    if (add) {
      lines(plotGrid, plotObj[i,], col = colPalette(n)[i], ...)
    } else {
      if (add != FALSE) warning("add is mis-specified and is taken to be FALSE.")
      plot(plotGrid, plotObj[i,], type='l', col = colPalette(n)[i],
           ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, main = main, ...)
    }
    if (is.numeric(pos.prob) & pos.prob >= 0 & pos.prob <= 1) {
      textidx <- which.min(abs(plotGrid - (xlim[1]*(1-pos.prob)+xlim[2]*pos.prob)))
      text(x = plotGrid[textidx], y = plotObj[i,textidx], labels = prob[i], col = colPalette(n)[i], cex=0.7)
    }
    if (n > 1) {
      for (i in 2:n) {
        lines(plotGrid, plotObj[i,], col = colPalette(n)[i], ...)
        if (is.numeric(pos.prob) & pos.prob >= 0 & pos.prob <= 1) {
          textidx <- which.min(abs(plotGrid - (xlim[1]*(1-pos.prob)+xlim[2]*pos.prob)))
          text(x = plotGrid[textidx], y = plotObj[i,textidx], labels = prob[i], col = colPalette(n)[i], cex=0.7)
        }
      }
    }

  } else {

    if(is.null(ylab)) ylab <- ifelse(obj %in% c(1,3), "Density", "Quantile")
    if(is.null(xlab)) xlab <- ifelse(obj %in% c(1,3), "", "Probability")

    if ((obj == 1 & is.list(x$dout)) | obj == 3) {
      if (obj == 1) {
        tmp <- x$dout
      } else tmp <- x$din
      n <- length(tmp)

      if (n == 1) {
        colPalette <- function(num) rep("black",num)
      } else {
        if (col.bar) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          
          #def.par <- par(no.readonly = TRUE)
          layout(t(1:2), widths = c(widrt,1))
        }
      }

      if (is.null(ylim)) {
        ylim <- sapply(tmp, function(d) range(d$y))
        ylim <- c(min(ylim[1,]), max(ylim[2,]))
      }
      if (is.null(xlim)) {
        xlim <- sapply(tmp, function(d) range(d$x))
        xlim <- c(min(xlim[1,]), max(xlim[2,]))
      }

      i <- 1
      plot(tmp[[i]]$x, tmp[[i]]$y, col=colPalette(n)[i], type='l',
           xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
      if (n > 1) {
        for (i in 2:n) {
          lines(tmp[[i]]$x, tmp[[i]]$y, col=colPalette(n)[i], ...)
        }
      }
    } else {
      if (obj == 1) {
        plotObj <- x$dout
        plotGrid <- x$dSup
      } else if (obj == 2) {
        plotObj <- x$qout
        plotGrid <- x$qSup
      } else {
        plotObj <- x$qin
        plotGrid <- x$qSup
      }
      if (is.null(ylim)) ylim <- range(plotObj)
      if (is.null(xlim)) xlim <- range(plotGrid)

      if (is.vector(plotObj)) {
        plotObj <- matrix(plotObj, nrow = 1)
      }
      n <- nrow(plotObj)
      
      if (n == 1) {
        colPalette <- function(num) rep("black",num)
      } else {
        if (col.bar) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          
          layout(t(1:2), widths = c(widrt,1))
        }
      }

      i <- 1
      plot(plotGrid, plotObj[i,], type='l', col = colPalette(n)[i],
           ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, main = main, ...)
      if (n > 1) {
        for (i in 2:n) {
          lines(plotGrid, plotObj[i,], col = colPalette(n)[i], ...)
        }
      }
    }
    if (obj %in% c(1,2)) {
      colVal <- x$xout
    } else colVal <- x$xin
    if (!is.vector(colVal)) {
      colVal <- colVal[,1]
      warning(ifelse(obj %in% c(1,2), "x$xout","x$xin")," is not a vector. Only the first column is used.")
    }

    if (n > 1) {
      if (col.bar) {
        color.bar(colVal = colVal, lut = colPalette(length(colVal)), nticks = nticks, ticks = ticks, title = col.lab)
        #par(def.par)
      }
    }
  }
}
