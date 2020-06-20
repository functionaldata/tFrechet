#' Generate color bar/scale.
#' @param colVal A numeric vector giving the variable values to which each color is corresponding. It overrides \code{min} (and \code{max}) if \code{min > min(colVal)} (\code{max < max(colVal)}).
#' @param colBreaks A numeric vector giving the breaks dividing the range of variable into different colors. It overrides \code{min} and \code{max}.
#' @param min A scalar giving the minimum value of the variable represented by colors.
#' @param max A scalar giving the maximum value of the variable represented by colors.
#' @param lut Color vector. Default is \cr
#' \code{colorRampPalette(colors = c("pink","royalblue"))(length(colBreaks)-1)}.
#' @param nticks An integer giving the number of ticks used in the axis of color bar.
#' @param ticks A numeric vector giving the locations of ticks used in the axis of color bar; it overrides \code{nticks}.
#' @param title A character giving the label of the variable according to which the color bar is generated.
#' @return No return value.
#' @export

color.bar <- function(colVal = NULL, colBreaks = NULL, min = NULL, max = NULL,
                      lut = NULL, nticks = 5, ticks = NULL, title = NULL) {
  if (is.null(colVal) & is.null(min) & is.null(colBreaks))
    stop("At least one of colVal, colBreaks and min must be specified.")

  if (!is.null(colBreaks)) {
    if (!is.null(lut)) {
      if (length(colBreaks) - length(lut) != 1)
        stop("length(colBreaks) - length(lut) must be equal to 1.")
    }
    min <- min(colBreaks)
    max <- max(colBreaks)
  } else {
    if (!is.null(min)) {
      if (is.null(colVal)) {
        if (is.null(max)) stop("max must be given when min is specified but neither are colVal and colBreaks.")
        colVal <- seq(min, max, length.out = length(lut))
      } else {
        if (min > min(colVal))
          min <- min(colVal)
        if (is.null(max)) {
          max <- max(colVal)
        } else if(max < max(colVal)) {
          max <- max(colVal)
        }
      }
    } else if (!is.null(colVal)) {
      if (!is.null(lut)) {
        if (length(colVal) != length(lut))
          stop("colVal and lut must have the same length.")
      }
      min <- min(colVal)
      max <- max(colVal)
    }
    colBreaks <- c(colVal[1], colVal[-1] - diff(colVal)/2, colVal[length(colVal)])
    colBreaks[1] <- min
    colBreaks[length(colBreaks)] <- max
  }
  if (is.null(lut))
    lut <- colorRampPalette(colors = c("pink","royalblue"))(length(colBreaks)-1)
  if (is.null(ticks)) ticks <- seq(min, max, len=nticks)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  title(xlab = title, line = 1)
  #plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:length(lut)) {
    rect(0,colBreaks[i],10,colBreaks[i+1], col=lut[i], border=NA)
  }
}
