#' Density Distribution
#' 
#' Density, cumulative distribution, and quantile functions for objects of
#'class "Density."
#'
#' @aliases Density_dist dDensity pDensity qDensity
#' 
#' @param x vector of quantities.
#' @param q vector of quantities.
#' @param p vector of probabilities.
#' @param object an object of class "Density."
#' 
#' @return The function \code{ddensity} gives the density, \code{pdensity} gives
#'the cumulative distribution, and \code{qdensity} gives the quantiles of the
#'empirical density estimate. One value is returned for each element of the first
#'argument.
#'
#' @rdname Density_dist
#' @export
dDensity <- function(x, object) {
  # Protect against NAs in y
  if(all(is.na(object$y))) {
    out <- rep(NA_real_, length(x))
  } else
    out <- approx(object$x, object$y, xout=x)$y
  if(object$Pzero > 0) {
    message("Conditional densities adjusted for zero values.")
    out[x == 0] <- NA
  }
  return(out)
}

#' @rdname Density_dist
#' @export
pDensity <- function(q, object) {
  xdiff <- diff(object$x[1:2])
  if(all(is.na(object$y))) {
    out <- rep(NA_real_, length(q))
  } else
    out <- approx(object$x, cumsum(object$y)*xdiff, xout=q)$y
  if(object$Pzero > 0) {
    out <- out + object$Pzero
    out[q == 0] <- NA
  }
  out <- pmin(out, 1) # Just in case
  return(out)
}

#' @rdname Density_dist
#' @export
qDensity <- function(p, object) {
  xdiff <- diff(object$x[1:2])
  if(all(is.na(object$y))) {
    out <- rep(NA_real_, length(p))
  } else {
    if(object$Pzero > 0){
      out <- approx(cumsum(object$y)*xdiff + object$Pzero, object$x, xout=p, 
                    ties="ordered")$y
      out[p < object$Pzero] <- 0
    } else
      out <- approx(cumsum(object$y)*xdiff, object$x, xout=p, 
                    ties="ordered")$y
  }
  return(out)
}