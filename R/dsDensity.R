#' Double Smoothing Density
#' 
#' An implementation of the double-smoothing density described by Yao.
#' 
#' @aliases dsDensity dsDensity.default dsDensity.Surv
#' 
#' @param x the variable for which density estimates are needed.
#' @param h the bandwidth, if missing, then Silverman’s rule of thumb
#'is used. See \code{Details}.
#' @param n the number of output density points.
#' @param from,to the left and right limits of the range of desired output points
#'where the density is to be estimated. If missing, then reasonable limits are
#'computed beyond the range of the data.
#'
#' @details Silverman’s rule of thumb, computes the bandwidth as 0.9 times the 
#'minimum of the standard deviation and the interquartile range divided by 1.34 
#'times the sample size to the negative one-fifth power. (From the documentation
#'for \code{bw.nrd0} in the stats package.) This tends to work better for double
#'smoothing than the normal density estimate, because it can handle smaller values
#'of bandwidth than traditional methods.
#'
#' @return A list of class "Density" with \code{n} values of x, the ordinate; 
#'and y the density values; summary statistics and tuning parameters; and the call. 
#'
#' @useDynLib Density dsden
#' @importFrom survival survreg
#' 
#' @rdname dsDensity
#' @export
dsDensity <- function(x, h, n=301, from, to) UseMethod("dsDensity")

#' @rdname dsDensity
#' @method dsDensity default
#' @export
dsDensity.default <- function(x, h, n=301, from, to) {
  call <- match.call()
  xnam <- deparse(substitute(x))[1L] # Just to prevent quirky behavior in text calls
  NOBS <- length(x)
  # Compute stats for x
  xstats <- c(mean=mean(x), var=var(x), median=median(x))
  # normalize x to make h consistently 1.0
  if(missing(h)) {
    h <- 0.9*min(c(sd(x), diff(quantile(x, c(.25, .75)))/1.34))*NOBS^-0.2
  }
  x <- x/h
  # Generate range if not specified
  if(missing(from)) {
    xoutmin <- 1.5*min(x) - 0.5*mean(x)
  } else
    xoutmin <- from/h
  if(missing(to)) {
    xoutmax <- 1.5*max(x) - 0.5*mean(x)
  } else
    xoutmax <- to/h
  xout <- seq(xoutmin, xoutmax, length=n)
  # the fortran interface
  hmm <- .Fortran("dsden",
                  x=x,
                  nobs=as.integer(NOBS),
                  wj=double(NOBS),
                  xout=xout,
                  yout=double(n),
                  nout=as.integer(n))
  # Normalize yout to sum to 1
  # I cannot get any reasonable theoretical value for c(h)
  # The divisor of 0.5908 is empirical 
  yout <- hmm$yout/NOBS^2/.5908
  retval <- list(x=xout*h, y=yout/h, xnam=xnam)
  # Compute stats for xout
  xdiff <- diff(retval$x[1:2])
  xoutmean <- sum(retval$x * retval$y) * xdiff
  xoutvar <- sum((retval$x - xoutmean)^2 * retval$y) * xdiff
  xoutcdf <- cumsum(retval$y) * xdiff
  xoutmedian <- approx(xoutcdf, retval$x, xout=0.5, ties="ordered")$y
  # Recover x
  x <- x*h
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  retval$stats <- stats
  retval$tuning <- c("h" = h)
  retval$xin=x
  retval$call <- call
  retval$Pzero <- 0
  retval$method <- "Double Smooth"
  class(retval) <- "Density"
  # Goodness of Fit Scores:
  loglik <- sum(log(dDensity(x, retval)))
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

#' @rdname dsDensity
#' @method dsDensity Surv 
#' @export
dsDensity.Surv <- function(x, h, n = 301, from, to) {
  call <- match.call()
  xnam <- deparse(substitute(x))[1L]
  xtype <- attr(x, "type")
  # Construct a dataset of the censored data
  if(xtype == "left") {
    xcen <- data.frame(Obs=x[,1], Cens=1-x[,2]) # reverse sense
  } else { # right censored
    stop("Right-censored data not supported in this version.")
    xcen <- data.frame(Obs=-x[,1], Cens=1-x[,2]) # flip data
  }
  # Construct the simple density object, just x and y components
  # Start with normal distribution 
  xmdl <- survreg(x ~ 1, dist="gaussian")
  xmn <- xmdl$coefficients
  if(xtype == "right")
    xmn <- -xmn
  xsd <- xmdl$scale
  NOBS <- nrow(xcen)
  xdens <- list(x=seq(from=xmn - 4*xsd, to=xmn + 4*xsd, length=351))
  xdens$y <- dnorm(xdens$x, xmn, xsd)
  xdens$Pzero <- 0
  # The estimated values (order not important!)
  xcen$Est <- with(xcen, rosEst(Obs, Cens, xdens))
  if(missing(h)) { # should be good enough
    h <- 0.9*min(c(sd(xcen$Est), diff(quantile(xcen$Est, c(.25, .75)))/1.34))*NOBS^-0.2
  }
  # First estimate using the default method
  xdens <- dsDensity(xcen$Est, h=h, n=n)
  # Update estimated Values
  xcen$Est <- with(xcen, rosEst(Obs, Cens, xdens))
  # Second estimate using the default method
  xdens <- dsDensity(xcen$Est, h=h, n=n)
  # Update estimated Values
  xcen$Est <- with(xcen, rosEst(Obs, Cens, xdens))
  # And the final density estimate
  # Generate range if not specified
  if(missing(from)) {
    xoutmin <- 1.5*min(xcen$Est) - 0.5*mean(xcen$Est)
  } else
    xoutmin <- from
  if(missing(to)) {
    xoutmax <- 1.5*max(xcen$Est) - 0.5*mean(xcen$Est)
  } else
    xoutmax <- to
  retval <- dsDensity(xcen$Est, h=h, n=n, from=xoutmin, to=xoutmax)
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
