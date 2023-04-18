#' Print Density Object
#' 
#' Print a general summary of a Density object.
#' 
#' @param x the Density object.
#' @param digits the number of significant digits to print in the numeric values.
#' @param \dots further arguments for other methods.
#' 
#' @method print Density
#' @export
print.Density <- function(x, digits=4, ...) {
  cat("Method:",x$method, "\nCall:\n")
  dput(x$call)
  cat("\nTuning arguments:\n")
  print(x$tuning, digits=digits)
  cat("\nSummary statisics:\n")
  print(x$stats, digits=digits)
  if(!is.null(x$logstats))
    print(x$logstats, digits=digits)
  if(x$Pzero > 0) {
    cat("Density estimates adjusted by ", signif(x$Pzero*100, digits), 
        " percent zero values.\n", sep="")
  } else {
    cat("\nGoodness of Fit:\n")
    cat("Normalized Za Score: ", round(x$ZaNorm, digits), "\n", sep="")
    print(x$loglik, digits=digits)
    cat("Variance Ratio: ", round(x$VarRat, digits), "\n", sep="")
    cat("Variability of Density Estimate: ", round(x$VarDen, digits), "\n", sep="")
  }
  if(!is.null(x$boot)) {
    area <- sum(x$boot$upper) - sum(x$boot$lower)
	if(is.null(x$boot$x)) {
	  area <- area * diff(x$x[1:2])
	} else
	  area <- area * diff(x$boot$x[1:2])
	cat("\nArea within the ", 100*x$bootCI, " percent confidence band: ",
	    round(area, digits), "\n", sep="")
  }
  invisible(x)
}
