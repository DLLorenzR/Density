#' Density Conversion
#'
#' Functions to convert selected kernel density objects to objects of class "Density."
#'
#' @aliases as.Density as.Density.default as.Density.Density as.Density.density
#'as.Density.sm as.Density.ssden as.Density.kde as.Density.logspline
#'as.Density.npdensity as.Density.pendensity as.Density.ash
#'
#' @param object the output from one of the selected kernel density methods.
#' @param \dots optional additional arguments required by certain methods.
#' @param xin the input data for the selected kernel density method (required).
#' @param xnam the name of the variable.
#' @param augment logical, if \code{TRUE}, then augment the truncated ends of
#'a spline density estimate to make the area as close to 1 as possible. If \code{FALSE},
#'then do not adjust the density estimate.
#' @param smooth logical, if \code{TRUE}, then smooth the original estimate. If
#'\code{FALSE}, then do not smooth the original estimate.
#' @param n the number of points to compute if \code{smooth} is \code{TRUE}.
#'
#' @note The default method assumes at least a list with components \code{x},
#'which contains the ordinates that are uniformly distributed and \code{y},
#'which are the density estimates. The integral of y by x must be close to 1.\cr
#'The default method also traps objects created by \code{sm} that have a specific
#'structure, but have no class that ties them to that function.\cr
#' The method for ash objects extrapolates the ends to a zero density value based
#'on the slope at the endpoints if "ash estimate nonzero outside interval ab."
#'Otherwise, the density at adjacent points is set to 0.
#'
#' @importFrom gss dssden
#' @importFrom logspline qlogspline dlogspline
#' @importFrom np npudens
#' @importFrom pendensity dpendensity
#' @importFrom stats var median spline approx dnorm sd quantile
#' 
#' @return An object of class "Density" with \code{object} included in the
#'object component to facilitate exploring the original object.
#'
#' @rdname as.Density
#' @export
as.Density <- function(object, ...) UseMethod("as.Density")

#' @rdname as.Density
#' @method as.Density default
#' @export
as.Density.default <- function(object, xin, ...) {
  # Trap unclassed objects created by sm or ash
  if(!is.null(object$eval.points)) {
    class(object) <- "sm"
    xnam <- deparse(substitute(xin))
    return(as.Density.sm(object, xin, xnam = xnam, ...)) 
  } else if(!is.null(object$m) && !is.null(object$kopt)) {
    if(missing(xin))
      stop("xin required for Density object")
    # Pack object with xin and xnam
    object$xnam <- deparse(substitute(xin)) 
    object$xin <- xin
    object$call <- match.call()
    class(object) <- "ash"
    return(as.Density.ash(object, ...)) # xin and xnam packed
  }
  call <- match.call()
  if(is.null(object$x) || is.null(object$y))
    stop("Input object does not contain either x or y")
  x <- object$x
  y <- object$y
  NOBS <- length(x)
  if(length(y) != NOBS)
    stop("The lengths of x and y do not match")
  xdiff <- diff(x[1:2])
  # Check uniformity of x
  if(max(abs(xdiff - range(diff(x))))/xdiff > 1.e-9)
    stop("The values of x are not uniform")
  # Check integral
  outcdf <- cumsum(y) * xdiff
  if(abs(outcdf[NOBS] - 1) > 0.01) # A pretty generous criterion
    stop("The integral of y by x is not close to 1")
  # OK, proceed as far as possible
  # Is input data available?
  if(missing(xin)) {
    xin <- object$xin
    if(is.null(xin))
      stop("xin required for Density object")
    xnam <- object$xnam
    if(is.null(xnam))
      xnam <- "X"
  } else {
    xnam <- deparse(substitute(xin))
    xin
  }
  # Compute stats for xin
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  # Compute stats for xout
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  if(is.null(xstats)) {
    stats <- as.data.frame(t(c(mean=xoutmean, var=xoutvar, median=xoutmedian)))
    row.names(stats) <- "Computed"
  } else {
    stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
    row.names(stats) <- c(xnam, "Computed")
  }
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=NA_real_, xin=xin, 
                 fit=NA_real_, call=call, method="Unknown", Pzero=0, 
                 object=object)
  class(retval) <- "Density"
  loglik <- sum(log(dDensity(xin, retval)))
  # Any kernel density methods without class include h or bw in object?
  # Goodness of fit scores
  attr(loglik, "nobs") <- NOBS
  attr(loglik, "df") <- NA_real_
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density Density
#' @export
as.Density.Density <- function(object, ...) {
  return(object)
}

#' @rdname as.Density
#' @method as.Density density
#' @export
as.Density.density <- function(object, xin, ...) {
  x <- object$x
  y <- object$y
  # Try to get input data
  # Is input data available?
  if(missing(xin)) {
    stop("xin required for Density object")
  } 
  xin
  NOBS <- length(xin)
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  xnam <- object$data.name
  # Compute stats for xout
  xdiff <- diff(x[1:2])
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  if(is.null(xstats)) {
    stats <- as.data.frame(t(c(mean=xoutmean, var=xoutvar, median=xoutmedian)))
    row.names(stats) <- "Computed"
  } else {
    stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
    row.names(stats) <- c(xnam, "Computed")
  }
  tuning <- c("bw" = object$bw)
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 call=object$call, method="Kernel Density", Pzero=0,
                 object=object)
  class(retval) <- "Density"
  # Construct log-likelihood and other goodnes of fit scores
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- NOBS
  # Check that kernel was not specified in the call
  if(is.null(object$call$kernel)) { # the default is gaussian
    # Compute ENP-- in this case tmp is the transpose of the true form
    tmp <- matrix(0, ncol=NOBS, nrow=NOBS)
    for(i in seq(NOBS))
      tmp[, i] <- dnorm((xin - xin[i])/object$bw)/NOBS/object$bw
    attr(loglik, "df") <- sum(diag(tmp)/colSums(tmp))
  } else 
    attr(loglik, "df") <- NA_real_
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density sm
#' @export
as.Density.sm <- function(object, xin, ...) {
  x <- object$eval.points
  y <- object$estimate
  if(missing(xin)) {
    stop("xin required for Density object")
  } 
  dots <- list(...)
  if(length(dots) && "xnam" %in% names(dots)) {
    xnam <- dots[["xnam"]]
  } else
    xnam <- deparse(substitute(xin)) 
  xin
  NOBS <- length(xin)
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  # Compute stats for xout
  xdiff <- diff(x[1:2])
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- c("h" = object$h)
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 call=object$call, method="Kernel Density", Pzero=0,
                 object=object)
  class(retval) <- "Density"
  # Construct log-likelihood and other goodness of fit scores
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- NOBS
  tmp <- matrix(0, ncol=NOBS, nrow=NOBS)
  for(i in seq(NOBS))
    tmp[, i] <- dnorm((xin - xin[i])/object$h)/NOBS/object$h
  attr(loglik, "df") <- sum(diag(tmp)/colSums(tmp))
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density ssden
#' @export
as.Density.ssden <- function(object, n=301, ...) {
  # Extract the necessary info from the object
  xnam <- object$terms$labels
  xin <- object$mf[[xnam]]
  x <- seq(object$terms[[xnam]]$phi$env$env$min, object$terms[[xnam]]$phi$env$env$max,
           length=n)
  # Proceed with computations
  y <- dssden(object, x)
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  # Compute stats for xout
  xdiff <- diff(x[1:2])
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- c("alpha" = object$alpha)
  # And the output object
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 fit=NA_real_, call=object$call, 
                 method="Smoothing Splines Density", Pzero=0,
                 object=object)
  class(retval) <- "Density"
  # Goodness of fit scores
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- length(xin)
  attr(loglik, "df") <- NA_real_
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density kde
#' @export
as.Density.kde <- function(object, xnam, ...) {
  # Extract the necessary info from the object
  xin <- object$x
  NOBS <- length(xin)
  if(missing(xnam))
    xnam <- object$names # This seems to always be "x"
  x <- object$eval.points
  y <- object$estimate
  # Proceed with computations
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  # Compute stats for xout
  xdiff <- diff(x[1:2])
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- c("h" = object$h) # 1-D value
  # And the output object
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 fit=NA_real_, call=object$call, 
                 method="Kernel Density", Pzero=0, object=object)
  class(retval) <- "Density"
  # Construct log-likelihood and other goodness of fit
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- NOBS
  tmp <- matrix(0, ncol=NOBS, nrow=NOBS)
  for(i in seq(NOBS))
    tmp[, i] <- dnorm((xin - xin[i])/object$h)/NOBS/object$h
  attr(loglik, "df") <- sum(diag(tmp)/colSums(tmp))
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density logspline
#' @export
as.Density.logspline <- function(object, xin, augment=TRUE, n=301, ...) {
  xnam <- deparse(object$call$x)
  # Compute x range
  NOBS <- object$samples
  # Compute range such that CDF is .995 (actually greater)
  xmin <- qlogspline(.0025, object)
  xmax <- qlogspline(.9975, object)
  x <- seq(xmin, xmax, length=n)
  y <- dlogspline(x, object)
  method <- "Logspline Density"
  xdiff <- diff(x[1:2])
  if(augment) {
    # Calculate area under the curve to compute the needed augmentation
    AreaNeeded <- 1 - sum(y) * xdiff
    # First and last nonzero densities are 1, n
    NumDiffs <- floor(2*AreaNeeded/(y[1L] + y[n])/xdiff)
    ypre <- seq(1, NumDiffs)/(NumDiffs + 1) * y[1L]
    xpre <- x[1L] - seq(NumDiffs, 1)*xdiff
    ypost <- seq(NumDiffs, 1)/(NumDiffs + 1) * y[n]
    xpost <- x[n] + seq(1, NumDiffs)*xdiff
    y <- c(ypre, y, ypost)
    x <- c(xpre, x, xpost)
    method <- "Augmented Logspline Density"
  }
  # Proceed with computations
  # Compute stats for xout
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  # Get xin if possible and complete stats
  if(missing(xin)) {
    stop("xin required for Density objects")
  } else {
    xin 
    xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
    stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
    row.names(stats) <- c(xnam, "Computed")
  }
  tuning <- c("nknots" = object$nknots) # 1-D value
  # Spline log-likelihoods are inflated, use computed loglike and df
  pick <- object$logl[,1] == object$nknots
  loglik <- object$logl[pick, 3]
  attr(loglik, "nobs") <- length(xin)
  attr(loglik, "df") <- object$nknots - 1
  class(loglik) <- "logLik"
  # And the output object
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 fit=NA_real_, call=object$call, 
                 method=method, Pzero=0, object=object, loglik=loglik)
  class(retval) <- "Density"
  # Other goodness of fit scores
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}
#' @rdname as.Density
#' @method as.Density npdensity
#' @export
as.Density.npdensity <- function(object, n=301, ...) {
  xnam <- object$bws$xnames
  xin <- object$eval[[xnam]]
  NOBS <- length(xin)
  # Compute x range
  xoutmin <- 1.5*min(xin) - 0.5*mean(xin)
  xoutmax <- 1.5*max(xin) - 0.5*mean(xin)
  x <- seq(xoutmin, xoutmax, length=n) # Adjusted in conversion to density
  edat <- data.frame(x=x)
  names(edat) <- xnam
  # Now reconstruct the call and append edat
  call <- object$bws$call
  call[[1]] <- as.name("npudens")
  m <- call
  m$edat <- edat
  m <- eval(m)
  y <- m$dens
  # Proceed with computations
  # Compute stats for xout
  xdiff <- diff(x[1:2])
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  # Compute stats
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- c("bw" = object$bw) # 1-D value
  # And the output object
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 fit=NA_real_, call=call, 
                 method="Kernel Density", Pzero=0, object=object)
  # Construct log-likelihood and other goodness of fits
  NPBS <- length(xin)
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- NOBS
  tmp <- matrix(0, ncol=NOBS, nrow=NOBS)
  for(i in seq(NOBS))
    tmp[, i] <- dnorm((xin - xin[i])/object$bw)/NOBS/object$bw
  attr(loglik, "df") <- sum(diag(tmp)/colSums(tmp))
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  class(retval) <- "Density"
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density pendensity
#' @export
as.Density.pendensity <- function(object, n=301, augment=TRUE, ...) {
  xnam <- deparse(object$call[[2]])
  xin <- object$values$y
  # Compute x range from range of help knots
  xrng <- range(object$splines$knots.val$help)
  x <- seq(xrng[1L], xrng[2L], length=n) 
  # Now reconstruct the call and append edat
  call <- paste0("Formula: ", deparse(object$call))
  y <- dpendensity(object, x)
  # Need this now
  xdiff <- diff(x[1:2])
  method <- "Penalized Density"
  if(augment) {
    # Calculate area under the curve to compute the needed augmentation
    AreaNeeded <- 1 - sum(y) * xdiff
    # Find first and last nonzero densities
    FirstNZ <- which(y > 0)[1L]
    LastNZ <- rev(which(y > 0))[1L]
    NumDiffs <- floor(2*AreaNeeded/(y[FirstNZ] + y[LastNZ])/xdiff)
    if(NumDiffs >= FirstNZ) {
      message("Can't augment in current version.")
    } else {
      y[seq(FirstNZ - NumDiffs, FirstNZ - 1)] <- 
        seq(1, NumDiffs)/(NumDiffs + 1) * y[FirstNZ]
      y[seq(LastNZ + NumDiffs, LastNZ + 1)] <- 
        seq(1, NumDiffs)/(NumDiffs + 1) * y[LastNZ]
      method <- "Augmented Penalized Density"
    }
  }
  # Proceed with computations
  # Compute stats for xout

  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  # Compute stats
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- data.frame("no.base" = (object$splines$K - 1)/2, 
                       "base" = object$splines$base,
                       "q" = object$splines$q, row.names="")
  # And the output object
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 fit=object$results$fitted, call=call, 
                 method=method, Pzero=0, object=object)
  # Goodness of fit
  loglik <- sum(log(dDensity(xin, retval)))
  attr(loglik, "nobs") <- length(xin)
  attr(loglik, "df") <- NA_real_
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  class(retval) <- "Density"
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname as.Density
#' @method as.Density ash
#' @export
as.Density.ash <- function(object, xin, smooth=TRUE, n=301, ...) {
  x <- object$x
  y <- object$y
  xdiff <- diff(x[1:2])
  # Is input data available? Should be packed from the redirected call from default
  if(missing(xin)) {
    xin <- object$xin
    if(is.null(xin))
      stop("xin required for Density object")
    xnam <- object$xnam
    if(is.null(xnam))
      xnam <- "X"
  } else {
    xnam <- deparse(substitute(xin))
    xin
  }
  if(is.null(object$call))
    object$call <- match.call()
  # Compute stats for xin
  xstats <- c(mean=mean(xin), var=var(xin), median=median(xin))
  # Adjust x and y for ends to 0 density
  if(object$ier == 0) { # ab values are 0, but represent 1/2 the spacing of x
    x <- c(x[1L] - xdiff, x, x[length(x)] + xdiff)
    y <- c(0, y, 0)
  } else { # Extrapolate
    dy <- diff(y[1:2])
    steps <- ceiling(y[1L]/dy)
    prex <- x[1L] - seq(steps, 1) * xdiff
    prey <- pmax(0, y[1L] - seq(steps, 1) * dy)
    nxy <- length(x)
    dy <- diff(y[c(nxy, nxy-1)]) # reverse slope
    steps <- ceiling(y[nxy]/dy)
    postx <- x[nxy] + seq(steps) * xdiff
    posty <- pmax(0, y[nxy] - seq(steps) * dy)
    x <- c(prex, x, postx)
    y <- c(prey, y, posty)
  }
  if(smooth) { # Interpolate
    xy <- spline(x, y, n=n)
    x <- xy$x
    y <- pmax(xy$y, 0)
    xdiff <- diff(x[1:2])
  }
  # Compute stats for xout
  xoutmean <- sum(x * y) * xdiff
  xoutvar <- sum((x - xoutmean)^2 * y) * xdiff
  xoutcdf <- cumsum(y) * xdiff
  xoutmedian <- approx(xoutcdf, x, xout=0.5, ties="ordered")$y
  if(is.null(xstats)) {
    stats <- as.data.frame(t(c(mean=xoutmean, var=xoutvar, median=xoutmedian)))
    row.names(stats) <- "Computed"
  } else {
    stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
    row.names(stats) <- c(xnam, "Computed")
  }
  # Check if ash output
  tuning <- c("m" = object$m, "kopt.1" = object$kopt[1L], 
              "kopt.2" = object$kopt[2L])
  retval <- list(x=x, y=y, xnam=xnam, stats=stats, tuning=tuning, xin=xin, 
                 call=call, method="ASH", Pzero=0, object=object)
  class(retval) <- "Density"
  loglik <- sum(log(dDensity(xin, retval)))
  # Goodness of fit
  attr(loglik, "nobs") <- length(xin)
  attr(loglik, "df") <- NA_real_
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}
