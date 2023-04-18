#' Bootstrap Resampling 
#' 
#' Generates bootstrap analysis of selected density estimates. The mean and 
#'confidence intervals are computed.
#'
#' @param x the numeric variable for which density estimates are needed.
#' @param method character the name of the density function to use. See 
#'\bold{Details}.
#' @param bootN the number of bootstrap iterations. The default is 200.
#' @param CI the confidence interval to estimate. The default is 0.90.
#' @param \dots additional arguments for \code{method}. The arguments must be fully
#'specified, abbreviations are not allowed.
#' 
#' @details \code{Method} must be one of "qqDensity," "dsDensity," "density,"
#'"ash," "ssden," "local," "logspline," "bootdata," or "merge." In
#'general, \code{method} is the name of the function, but "ash" indicates 
#'\code{ash1} and "local" indicates \code{density.lf}. For specialized density
#'estimates, "bootdata" resamples the data so that consistent bootstrap analyses
#'can be computed and "merge" merges the bootstrap estimates into a Density object.
#'An example is provided in the Overview vignette.\cr
#' The only option for kernel density estimation is "density." The results of the 
#'various methods are nearly identical and the bandwidth can be set from one of 
#'the other methods.\cr
#' For \code{method} = "merge," \code{x} must be the output matrix of the resampled
#'density estimates and the Density object must be supplied in an argument called
#'\code{Density}.\cr
#' The \code{mixDensity} and \code{pendensity} methods are not implemented 
#'because the fit cannot be specified to ensure that the method assumes that the
#'resampled data come from the same population as the original. The same kind of 
#'issues can also affect \code{logDensity}, but the user can set up to ensure
#'that assumption. See the example in the Overview vignette.
#' 
#' @note In general, \code{set.seed} should be set before calling 
#'\code{bootDensity} to ensure reproducibility and consistency across different
#'functions.
#' 
#' @return An object of class "Density" with the bootstrap analysis in the 
#'components \code{boot} and \code{bootCI}.
#'
#' @importFrom ash ash1 bin1
#' @importFrom gss ssden
#' @importFrom locfit density.lf
#' @importFrom logspline logspline dlogspline
#' @importFrom utils capture.output
#' 
#' @export
bootDensity <- function(x, method, bootN=200, CI=0.90, ...) {
  method <- match.arg(method, c("qqDensity", "dsDensity", "density", "ash",
                                "ssden", "local", "logspline",
                                "bootdata", "merge"))
  xnam <- deparse(substitute(x))
  NOBS <- length(x)
  dots <- list(...)
  # Create the matrix of bootstrap samples 
  # (required because some methods may generate random numbers)
  bootSamp <- matrix(0, ncol=bootN, nrow=NOBS)
  for(i in seq(bootN))
    bootSamp[,i] <- sample(x, NOBS, replace=TRUE)
  if(method == "bootdata") {
    return(bootSamp)
  }
  if(method == "merge") {
    if(is.null(retval <- dots$Density))
      stop("The Density object must be supplied in an argument called Density.")
    # Extract the stats
    ymean <- rowMeans(x)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(x, 1, quantile, probs=ci)
    yll <- apply(x, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "dsDensity") {
    cat("Take a break. Bootstrapping dsDensity is extremely slow!\n")
    # First build the Density object
    retval <- dsDensity(x, ...)
    # Fix xnam
    retval$xnam <- xnam
    # Extract the necessary info
    n <- length(retval$x)
    h <- retval$tuning
    from <- retval$x[1L]
    to <- retval$x[n]
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=n)
    for(i in seq(bootN))
      bootY[,i] <- dsDensity(bootSamp[,i], h=h, n=n, from=from, to=to)$y
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "qqDensity") {
    # First build the Density object
    retval <- qqDensity(x, ...)
    # Fix xnam
    retval$xnam <- xnam
    # Extract/update the necessary info
    n <- length(retval$x)
    from <- retval$x[1L]
    to <- retval$x[n]
    alpha <- retval$tuning[["alpha"]]
    k <- retval$tuning[["k"]]
    shape <- retval$tuning[["shape"]]
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=n)
    for(i in seq(bootN))
      bootY[,i] <- qqDensity(bootSamp[,i], alpha=alpha, k=k, n=n, from=from, to=to,
                             shape=shape, plot.fit=FALSE)$y
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "density") {
    # First execute the function and build the Density object
    object <- density(x, ...)
    retval <- as.Density(object, xin=x)
    # Fix xnam
    retval$xnam <- xnam
    # Extract/update the necessary info
    n <- length(retval$x)
    dots$from <- retval$x[1L]
    dots$to <- retval$x[n]
    dots$bw <- retval$tuning
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=n)
    for(i in seq(bootN)) {
      dots$x <- bootSamp[,i]
      bootY[,i] <- do.call(density, dots)$y
    }
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "ash") {
    # First extract dots arguments to build the call
    # Any bin1 args?
    ab <- FALSE
    nbin <- 50
    if("ab" %in% names(dots))
      ab <- TRUE
    if("nbin"  %in% names(dots))
      nbin <- dots$nbin
    if(ab) {
      object <- bin1(x, ab=dots$ab, nbin=nbin)
    } else 
      object <- bin1(x, nbin=nbin)
    m <- 5
    kopt <- c(2,2)
    if("m" %in% names(dots))
      m <- dots$m
    if("kopt"  %in% names(dots))
      kopt <- dots$kopt
    object <- ash1(object, m=m, kopt=kopt)
    retval <- as.Density(object, xin=x)
    # Fix xnam
    retval$xnam <- xnam
    # Extract/update the necessary info
    ab <- object$ab
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=nbin)
    # Discard any  "ash estimate nonzero outside interval ab"
    capture.output(
      for(i in seq(bootN)) {
        bootY[,i] <- ash1(bin1(bootSamp[,i], ab=ab, nbin=nbin), m=m, kopt=kopt)$y
      },
      file=NULL)
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll, x=object$x)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "ssden") {
    # First extract dots arguments to build the call
    # Is data in dots? If so, trick into local dataset
    if("data" %in% names(dots)) {
      x <- dots$data[[xnam]]
      dots$data <- list(x=x)
    }
    # How about domain?
    if("domain" %in% names(dots)) {
      xrng <- dots$domain[[xnam]]
      dots$domain <- list(x=xrng)
    }
    dots$formula <- ~ x
    object <- do.call(ssden, dots)
    retval <- as.Density(object)
    # Fix xnam
    retval$xnam <- xnam
    # Extract/update the necessary info
    dots$domain <- object$domain
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=301) # The default for as.Density
    for(i in seq(bootN)) {
      x <- bootSamp[,i]
      dots$data <- list(x=x)
      bootY[,i] <- as.Density(do.call(ssden, dots))$y
    }
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "local") {
    # Check on width and use the modified width, not default
    if(is.null(dots$width)) {
      dots$width <- diff(range(x))/(logb(length(x), base = 2) + 1)/0.5
    }
    if(is.null(dots$n)) {
      message("The default number of output points for density.lf may be too small, increasing to 301.")
      n <- dots$n <- 301
    } else
      n <- dots$n
    ft <- fromTo(x)
    if(is.null(dots$from))
      dots$from <- ft[1L]
    if(is.null(dots$to))
      dots$to <- ft[2L]
    dots$x <- x
    retval <- as.Density(do.call(density.lf, dots), x)
    retval$xnam <- xnam
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=n) # The default for as.Density
    for(i in seq(bootN)) {
      dots$x <- bootSamp[,i]
      bootY[,i] <- do.call(density.lf, dots)$y
    }
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
  if(method == "logspline") {
    # No special concerns with dots
    object <- logspline(x, ...)
    n <- 301
    # Fix knots or number of knots?
    # The actual knots seem to produce fewer warnings
    dots$knots <- object$knots
    retval <- as.Density(object, xin=x, augment=TRUE, n=n)
    retval$xnam <- xnam
    # The matrix for density estimates and compute
    bootY <- matrix(0, ncol=bootN, nrow=length(retval$x)) # The default for as.Density
    for(i in seq(bootN)) {
      dots$x <- bootSamp[,i]
      bootY[,i] <- dlogspline(retval$x, do.call(logspline, dots))
    }
    # Extract the stats
    ymean <- rowMeans(bootY)
    # Adjust CI to reflect 2 sides
    ci <- 1 - (1 - CI)/2
    yul <- apply(bootY, 1, quantile, probs=ci)
    yll <- apply(bootY, 1, quantile, probs=1-ci)
    retval$boot <- data.frame(mean=ymean, upper=yul, lower=yll)
    retval$bootCI <- CI
    return(retval)
  }
}