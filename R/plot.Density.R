#' Plot a Density Object
#' 
#' Plot the computed density and points for a density object.
#' 
#' @param x object of class "Density."
#' @param xlab the x-axis label. If missing, then the label is the variable name.
#' @param ylab the y-axis label. If missing, then the label is automatically 
#'generated from information in \code{x}. 
#' @param type the type of plot, default is "l," which draws the fit as a line.
#' @param \dots additional arguments modifying the plot. 
#' 
#' @importFrom graphics abline segments par text polygon mtext
#' 
#' @return The object is returned invisibly.
#' 
#' @method plot Density
#' @export
plot.Density <- function(x, xlab, ylab, type="l", ...) {
  if(missing(xlab)) {
    if(x$Pzero > 0) {
      xlab <- paste0("Nonzero ", x$xnam)
    } else
      xlab <- x$xnam
  }    
  # Get y values to adjust for nice numbers
  y <- x$y
  ymax <- max(y)
  if(missing(ylab)) {
    if(ymax < 1e-5) {
      ylab <-  "Density * 1000000"
      y <- y * 1000000
    } else if(ymax < 1e-2) {
      ylab <-  "Density * 1000"
      y <- y * 1000
    } else
      ylab <- "Density"
  }
  # Show truncated distributions by adding droplines
  xx <- x$x
  if(y[1L]/max(y) > 0.001) {
    xx <- c(xx[1L], xx)
    y <- c(0, y)
  }
  if(y[length(y)]/max(y) > 0.001) {
    xx <- c(xx, xx[length(xx)])
    y <- c(y, 0)
  }
  # Get ylimits if not specified
  if("ylim" %in% names(list(...)))
    plot(xx, y, xlab=xlab, ylab=ylab, type=type, ...)
  else {
    if(is.null(x$boot)) {
      ylim <- range(pretty(y))
    } else
      ylim <- range(pretty(c(0, x$boot$upper)))
    plot(xx, y, xlab=xlab, ylab=ylab, type=type, ylim=ylim, ...)
  }
  if(!is.null(x$boot)) {
    xxx <- x$boot$x
    if(is.null(xxx))
      xxx <- x$x
    polygon(c(xxx, rev(xxx)), c(x$boot$upper, rev(x$boot$lower)), col='gray90',
	        border = NA)
    lines(xxx, x$boot$mean, lwd=2, col='white')
    lines(xx, y) # recover the obliterated line
	mtext(paste0(x$bootCI*100, "% Confidence Band"), line=-1.2, adj=0.05)
  }
  abline(h=0, col='gray')
  if(!is.null(x$xin)) {
	xin <- x$xin
    dupes <- duplicated(xin)
    xin[dupes] <- jitter(xin[dupes])
    segments(xin, 0, xin, par("usr")[3L]/3)
  }
  invisible(x)
}