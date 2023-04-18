#' Density Estimate for Log-Transformed Data
#' 
#' Some highly skewed data or data that have a minimum of 0 benefit from a 
#'log-transformation and possible probability adjustment for 0 values.
#' 
#' @aliases logDensity logDensity.default logDensity.Surv
#' 
#' @param x the variable for which density estimates are needed.
#' @param method the name of the density estimation method to use on the 
#'log-transformed data. The variable for which density estimates are needed 
#'must be \code{x}.
#' @param n the number of output density points for the-back-transformed values.
#' @param \dots any number of additional arguments for \code{method}
#' 
#' @details Zero values are extracted from \code{x} and the remaining values are
#'log-transformed. The back-transformed probability distribution is adjusted for
#'the proportion of zero values.
#'
#' @return An object of class "Density."
#' @rdname logDensity
#' @export
logDensity <- function(x, method, n=301, ...) UseMethod("logDensity")

#' @rdname logDensity
#' @method logDensity default
#' @export
logDensity.default <- function(x, method, n=301, ...) {
  if(!typeof(method) == "character")
    method <- deparse(substitute(method))
  xnam <- deparse(substitute(x))
  dots <- list(...)
  call <- match.call()
  if(any(x < 0))
    stop("x values cannot be less than 0")
  NOBS <- length(x)
  Nzero <- sum(x == 0)
  xgt0 <- x[x > 0]
  logxgt0 <- log(xgt0)
  dots$x <- logxgt0
  retval <- as.Density(do.call(method, dots), xin=logxgt0) # xin required for some methods
  # Fix some stuff in retval
  retval$call <- call
  stats <- retval$stats
  if(nrow(stats) == 2L) {
    row.names(stats) <- c(paste0("Log ", xnam), "Log Computed")
  } else { # Need to compute the stats
    stats <- rbind(stats, c(mean=mean(logxgt0), var=var(logxgt0),
                            median=median(logxgt0)))
    # And we need to flip the order.
    stats <- stats[c(2,1), ]
    row.names(stats) <- c(paste0("Log ", xnam), "Log Computed")
  }
  retval$logstats <- stats
  retval$Pzero <- Pzero <- Nzero/NOBS
  retval$xin <- x
  retval$xnam <- xnam
  # Back transform the density
  xt <- exp(retval$x)
  yt <- retval$y/xt
  # Calculate a reasonable upper limit for xout based on slope at peak value
  dMax <- suppressMessages(dDensity(max(logxgt0), retval))
  dMax1 <- suppressMessages(dDensity(max(logxgt0) - .1, retval))
  Fact <- max(exp(dMax/(dMax1 - dMax)/9), 1.25)
  # Prepend 0 if needed
  if(Pzero > 0 || mean(logxgt0)/sd(logxgt0) < 4.5) {
    xt <- c(0, xt)
    yt <- c(0, yt)
    xout <- seq(0, max(x)*Fact, length=n)
  } else {
    xout <- seq(min(xt), max(x)*Fact, length=n)
  }
  yout <- approx(xt, yt, xout=xout)$y
  # Adjust by Pzero
  retval$x <- xout 
  retval$y <- yout*(1-Pzero)
  # Protect against quirky value for y[2]
  if(retval$y[2L] > retval$y[3L])
    retval$y[2L] <- retval$y[3L]/2
  # Compute stats on raw data and density estimate, if Pzero = 0
  xstats <- c(mean=mean(x), var=var(x), median=median(x))
  # Compute stats for xout
  xdiff <- diff(xout[1:2])
  if(Pzero == 0) {
    xoutmean <- sum(xout * yout) * xdiff
    xoutvar <- sum((xout - xoutmean)^2 * yout) * xdiff
    xoutcdf <- cumsum(yout) * xdiff
    xoutmedian <- approx(xoutcdf, xout, xout=0.5)$y
  } else {
    xoutmean <- xoutvar <- xoutmedian <- NA_real_
  }
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  retval$stats <- stats
  # Goodness of Fit Scores:
  loglik <- sum(log(dDensity(x, retval)))
  attr(loglik, "nobs") <- NOBS
  attr(loglik, "df") <- NA_real_ # Until we figure out how to capture info from method
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}

#' @rdname logDensity
#' @method logDensity Surv 
#' @export
logDensity.Surv <- function(x, method, n=301, ...) {
  call <- match.call()
  xnam <- deparse(substitute(x))[1L]
  xtype <- attr(x, "type")
  # Construct a dataset of the censored data
  if(xtype == "left") {
    xcen <- data.frame(Obs=x[,1], Cens=1-x[,2]) # reverse sense
  } else # right censored
    xcen <- data.frame(Obs=-x[,1], Cens=1-x[,2]) # flip data
  # Construct the simple density object, just x and y components
  # Start with normal distribution 
  xmdl <- survreg(x ~ 1, dist="loggaussian")
  xmn <- xmdl$coefficients
  if(xtype == "right")
    xmn <- -xmn
  xsd <- xmdl$scale
  NOBS <- nrow(xcen)
  xdens <- list(x=seq(from=xmn - 4*xsd, to=xmn + 4*xsd, length=351))
  xdens$y <- dnorm(xdens$x, xmn, xsd)
  xdens$Pzero <- 0
  # The estimated values (order not important!)
  xcen$Est <- with(xcen, rosEst(log(Obs), Cens, xdens))
  # First estimate using the default method
  xdens <- as.Density(do.call(method, list(x=xcen$Est, n=n, ...)), xin=x)
  # Update estimated Values
  xcen$Est <- with(xcen, rosEst(log(Obs), Cens, xdens))
  # Second estimate using the default method
  xdens <- as.Density(do.call(method, list(x=xcen$Est, n=n, ...)), xin=x)
  # Update estimated Values # this time in back-transformed units
  xcen$Est <- exp(with(xcen, rosEst(log(Obs), Cens, xdens)))
  # And the final density estimate
  retval <- logDensity(xcen$Est, method=method, n=n, ...)
  # Update stats
  retval$stats[1L,1L] <- xmn
  retval$stats[1L,2L] <- xsd^2
  row.names(retval$stats)[1] <- xnam
  if(xtype == "right") { # flip sense of  density and xin
    retval$x <- rev(-retval$x)
    retval$y <- rev(retval$y)
    retval$xin <- -retval$xin
    # And the stats
    retval$stats[,1L] <- -retval$stats[,1L]
    retval$stats[,3L] <- -retval$stats[,3L]
  }
  retval$xnam <- xnam
  retval$Pzero <- 0
  retval$Surv <- x
  retval$call <- call
  return(retval)
}
