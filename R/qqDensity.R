#' Q-Q Density Estimate
#' 
#' A proposed density estimate computed by fitting a smooth, monotonically 
#'increasing curve through the quantiles of a distribution related to 
#'the observed values to compute a cumulative probability estimate (CDF). The 
#'density estimate is the derivative of the CDF. The purpose is to better 
#'describe the observed data.
#' 
#' @aliases qqDensity qqDensity.default qqDensity.Surv
#' 
#' @param x the variable for which density estimates are needed.
#' @param alpha the plotting position tuning parameter. See \bold{Details}.
#' @param k the number of knots in the monotonic cubic spline.  See \bold{Details}.
#' @param n the number of output density points.
#' @param from,to the left and right limits of the range of desired output points
#'where the density is to be estimated. If missing, then reasonable limits are
#'computed beyond the range of the data.
#' @param \dots Additional arguments for specific methods.
#' @param shape the general description of the shape of the distribution of 
#'\code{x}. Must be one of "auto" (default), "normal," "left," "right," "platy,"
#'"lepto," or "multimodal." See \bold{Details}.
#' @param plot.fit logical controlling a plot of the Q-Q fit. If \code{TRUE}, then
#'the plot is created on the current plotting device. Otherwise, no plot is 
#'created.
#' 
#' @details The value of \code{alpha} can affect the tails of the
#'density estimate---small values, near 0, have more narrow tails than those 
#'near 0.5, which produces wider tails. The default is 0.5\cr
#' The argument \code{k} determines the number of knots in the cubic spline fit.
#'The minimum is 4 larger the values can decrease the smoothness of the fit.
#'The default is 5, which appears to be a good compromise.\cr
#' The default \code{shape} method selects a distribution based on the 
#'statistics of \code{x}. If the skewness is greater than 0.5, then a gamma 
#'distribution is selected, specified by "right." If the skewness is less than 
#'-0.5, then a flipped gamma distribution is selected, specified by "left."
#'If the kurtosis is greater than 3.67, then a t distribution is selected,
#'specified by "lepto." If the kurtosis is less than 2.5, then a generalized 
#'normal distribution is selected, specified by "platy." If the statistics of \code{x} 
#'meets none of those criteria, then a normal distribution is used, specified by
#'"normal." If the shape is multimodal, then the user must specify "multimodal."
#'Can be abbreviated. The Surv method always uses a normal distribution.
#' 
#' @return A list of class "Density" with \code{n} values of x, the ordinate; 
#'and y the density values; stats, the input and derived (computed) summary 
#'statistics; h the value of h used in the estimate; and call, the call. 
#'
#' @import mgcv
#' @import gnorm
#' @importFrom stats pgamma qgamma qt pt qnorm pnorm ppoints 
#' @importFrom graphics lines
#' 
#' @rdname qqDensity
#' @export
qqDensity <- function(x, alpha=0.5, k=5, n=301, from, to, ...)
  UseMethod("qqDensity")

#' @rdname qqDensity
#' @method qqDensity default
#' @export
qqDensity.default <- function(x, alpha=0.5, k=5, n=301, from, to, 
                              shape="auto", plot.fit=TRUE, ...) {
  # simple logistic distribution functions
  # dslogis <- function(x) exp(-x)/(1 + exp(-x))^2 # Not needed
  qslogis <- function(p) log(p/(1-p))
  pslogis <- function(q) 1/(1 + exp(-q))
  # begin execution
  call <- match.call()
  shape <- match.arg(shape, c("auto", "normal", "left", "right", "platy", "lepto",
                     "multimodal"))
  xnam <- deparse(substitute(x))[1L] # Just to prevent quirky behavior in text calls
  xin <- sort(x)
  NOBS <- length(x)
  # Compute stats for x
  xstats <- c(mean=mean(x), var=var(x), median=median(x))
  # Skewness and kurtosis and append to xstats:
  xskew <- sum((x - xstats[1L])^3)/NOBS/mean((x - xstats[1L])^2)^1.5
  xkurt <- sum((x - xstats[1L])^4)/NOBS/mean((x - xstats[1L])^2)^2
  xstats <- c(xstats, skew=xskew, kurtosis=xkurt)
  # Establish distribution
  if(shape == "multimodal") {
    qDist <- qslogis
    pDist <- pslogis
    method <- "Quantile-Q-logistic Fit"
  } else if(shape == "right" || (shape == "auto" && xskew > 0.5)) {
    qDist <- function(p) qgamma(p, 9)
    pDist <- function(q) pgamma(q, 9)
    method <- "Quantile-Q-gamma Fit"
    shape <- "right"
  } else if(shape == "left" || (shape == "auto" && xskew < -0.5)) {
    qDist <- function(p) -qgamma(1 - p, 9)
    pDist <- function(q) 1 - pgamma(-q, 9)
    method <- "Quantile-Q-flipped-gamma Fit"
    shape <- "left"
  } else if(shape == "lepto" || (shape == "auto" && xkurt > 3.67)) {
    qDist <- function(p) qt(p, 10)
    pDist <- function(q) pt(q, 10)
    method <- "Quantile-Q-t Fit"
    shape <- "lepto"
  } else if(shape == "platy" || (shape == "auto" && xkurt < 2.5)) { 
    # The generalized normal distribution with beta > 2 (about 3) should work
    # It is available in the gnorm package
    qDist <- function(p) qgnorm(p, beta=2.5)
    pDist <- function(q) pgnorm(q, beta=2.5)
    method <- "Quantile-Q-generalized Fit"
    shape <- "platy"
  } else { # default to normal
    qDist <- qnorm
    pDist <- pnorm
    method <- "Quantile-Q-normal Fit"
    shape <- "normal"
  }
  ppx <- ppoints(xin, a=alpha)
  qqx <- qDist(ppx)
  if(plot.fit) { # adjust visual Quantiles to match data
    mnqqx <- mean(qqx)
    sdqqx <- sd(qqx)
    plot(xin, (qqx - mnqqx)/sdqqx*sqrt(xstats[2L]) + xstats[1L], 
         xlab=xnam, ylab="Theoretical Quantiles")
  }
  # Generate range if not specified
  if(missing(from)) {
    xoutmin <- 1.5*min(x) - 0.5*mean(x)
  } else
    xoutmin <- from
  if(missing(to)) {
    xoutmax <- 1.5*max(x) - 0.5*mean(x)
  } else
    xoutmax <- to
  xout <- seq(xoutmin, xoutmax, length=n+1) # Adjusted in conversion to density
  df <- data.frame(x=xin, y=qqx)
  # Generate unconstrained model for parameters needed later
  # Number of knows needed may depend on the variability of the data
  # May need some measure of the uniformity of the increases
  unc <- gam(y ~ s(x, k=k, bs="cr"), data=df)
  # Get the penalty matrix
  sm <- smoothCon(s(x, k=k, bs="cr"), data=df, knots=NULL)[[1]]
  #  Create the constraint matrix
  Fmat <- mono.con(sm$xp, up=TRUE)
  # Now construct then list needed for pcsl:
  Glist <- list(X=sm$X, C=matrix(0,0,0), sp=unc$sp, p=sm$xp, y=df$y, 
                w=rep(1, length(df$y)), Ain=Fmat$A, bin=Fmat$b, S=sm$S, off=0)
  # Fit the spline producing the coefficients
  Coef <- pcls(Glist)
  # Compute the  monotonic fitted values
  fit <- (Predict.matrix(sm, data.frame(x=xout)) %*% Coef)[,1]
  if(plot.fit)
    lines(xout, (fit - mnqqx)/sdqqx*sqrt(xstats[2L]) + xstats[1L])
  # Convert to probability
  pfit <- pDist(fit)
  # And to density
  dfit <- diff(pfit)/diff(xout[1:2])
  # Adjust xout
  xout <- (xout[-1L] + xout[-n])/2
  # Compute stats for xout
  xdiff <- diff(xout[1:2])
  xoutmean <- sum(xout * dfit) * xdiff
  xoutvar <- sum((xout - xoutmean)^2 * dfit) * xdiff
  xoutcdf <- cumsum(dfit) * xdiff
  xoutmedian <- approx(xoutcdf, xout, xout=0.5, ties="ordered")$y
  xoutskew <- (sum((xout - xoutmean)^3 * dfit) * xdiff)/xoutvar^1.5
  xoutkurt <- (sum((xout - xoutmean)^4 * dfit) * xdiff)/xoutvar^2
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian, 
                                         xoutskew, xoutkurt)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- data.frame("alpha" = alpha, "k" = k, "shape" = shape, row.names="")
  retval <- list(x=xout, y=dfit, xnam=xnam, stats=stats, tuning=tuning, xin=x, 
                 call=call, method=method, Pzero=0, fit=fit, qqx=qqx)
  class(retval) <- "Density"
  # Goodness of fit scores
  loglik <- sum(log(dDensity(x, retval)))
  attr(loglik, "nobs") <- NOBS
  attr(loglik, "df") <- k
  class(loglik) <- "logLik"
  retval$loglik <- loglik
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}  


#' @rdname qqDensity
#' @method qqDensity Surv 
#' @export
qqDensity.Surv <- function(x, alpha = 0.5, k = 5, n = 301, from, to,
                           plot.fit = TRUE, ...) {
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
  # First estimate using the default method
  xdens <- qqDensity(xcen$Est, alpha = alpha, k = k, n = n, shape="normal", 
                     plot.fit = FALSE)
  # Update estimated Values
  xcen$Est <- with(xcen, rosEst(Obs, Cens, xdens))
  # Second estimate using the default method
  # Generate range if not specified
  if(missing(from)) {
    from <- xdens$x[1L]
  } 
  if(missing(to)) {
    to <- xdens$x[n]
  }
  xdens <- qqDensity(xcen$Est, alpha = alpha, k = k, n = n, from=from, to=to, 
                     shape="normal", plot.fit = FALSE)
  # Update estimated Values
  xcen$Est <- with(xcen, rosEst(Obs, Cens, xdens))
  # And the final density estimate
  retval <- qqDensity(xcen$Est, alpha = alpha, k = k, n = n, from=from, to=to,
                      shape="normal", plot.fit = plot.fit)
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
