#' Mixed Gaussian Density Estimate
#' 
#' A density estimation method that moels the data as a mixture of two or more
#'Gaussian distributions. This method does not support censored data.
#' 
#' @param x the variable for which density estimates are needed.
#' @param Nmix the number of Gaussian mixtures in \code{x}. If 1, then the density
#'is computed assuming the data are normally distributed. The default is 2.
#' @param n the number of output density points.
#' @param \dots Additional arguments for specific methods. Not used currently.
#' 
#' @details The value of \code{alpha} can affect the peakedness of the
#'density estimate---small values, near 0, have higher peaks than those 
#'near 0.5, which produces wider tails. The default is 0.5\cr
#' The argument \code{k} determines the number of knots in the cubic spline fit.
#'The minimum is 4 larger the values can decrease the smoothness of the fit.
#'The default is 5, which appears to be a good compromise.\cr
#' The default \code{qqDensity} method selects a distribution based on the 
#'statistics of \code{x}. If the skewness is greater than 0.5, then a gamma 
#'distribution is selected. If the skewness is less than -0.5, then a flipped 
#'gamma distribution is selected. If the kurtosis is greater than 3.67, then a 
#'t distribution is selected. If the kurtosis is less than 2.5, then a 
#'generalized normal distribution is selected. If the statistics of \code{x} 
#'meets non of those criteria, then a normal distribution is used. The Surv
#'method always uses a logistic distribution.
#' 
#' @return A list of class "Density" with \code{n} values of x, the ordinate; 
#'and y the density values; stats, the input and derived (computed) summary 
#'statistics; h the value of h used in the estimate; and call, the call. In 
#'addition to the standard components of a Density object, \code{posterir},
#'the posterior probabilities; \code{logLike}, the log likelihood; and 
#'\code{mix}, the mean, standard deviation, and \code{posterior} the proportion 
#' of each Gaussian distribution for every point.
#' 
#' @rdname mixDensity
#' @export
mixDensity <- function(x, Nmix=2, n=301, ...) {
  # Define the expectation and maximization functions.
  eStep <- function(x, Mix) {
    NOBS <- length(x)
    NMIX <- nrow(Mix)
    denX <- matrix(0, nrow=NOBS, ncol=NMIX)
    for(i in seq(NMIX)) 
      denX[,i] <- dnorm(x, Mix[i, 1L], sqrt(Mix[i, 2L])) * Mix[i,3L]
    denSumX <- rowSums(denX)
    return(list("loglik" = sum(log(denSumX)),
                "posterior" = denX/denSumX))
  }
  mStep <- function(x, posterior) {
    NOBS <- length(x)
    NMIX <- ncol(posterior)
    postSum <- colSums(posterior)
    means <- colSums(posterior * x)/postSum
    vars <- colSums(posterior * (x - rep(means, each=NOBS))^2)/postSum
    alphas <- postSum/NOBS
    retval <- cbind(mean=means, var=vars, alpha=alphas)
  }
  # Begin execution
  call <- match.call()
  xnam <- deparse(substitute(x))[1L] # Just to prevent quirky behavior in text calls
  xin <- sort(x)
  NOBS <- length(x)
  # Compute stats for x
  xstats <- c(mean=mean(x), var=var(x), median=median(x))
  xoutmin <- 1.5*min(x) - 0.5*mean(x)
  xoutmax <- 1.5*max(x) - 0.5*mean(x)
  xout <- seq(xoutmin, xoutmax, length=n)
  if(Nmix > 1) {
    # Compute arbitrary groups and stats for each group
    xgrp <- cut(x, Nmix, labels=FALSE)
    Mix <- cbind(mean=tapply(x, xgrp, mean), var=tapply(x, xgrp, var), alpha= 1/Nmix)
    loop <- 0
    lastLL <- -Inf
    repeat {
      e <- eStep(x, Mix)
      Mix <- mStep(x, e$posterior)
      # Trap collapse of one of the Gaussian mixtures
      if(any(is.na(Mix))) {
        e$loglik <- NA_real_
        warning("Failed to converge, try reducing Nmix")
        break
      }
      if(abs(lastLL - e$loglik) < 1e-6) break
      lastLL <- e$loglik
    }
  } else { # Single normal distribution
    Mix <- matrix(c(xstats[1L], xstats[2L], 1), nrow=1)
    colnames(Mix) <- c("mean", "var", "alpha")
    e <- list(loglik=sum(log(dnorm(x, Mix[1L], sqrt(Mix[2L])))),
              posterior=matrix(1, nrow=NOBS, ncol=1))
  }
  # Compute yout positions
  yout <- numeric(n)
  for(i in seq(Nmix))
    yout <- yout + dnorm(xout, Mix[i,1L], sqrt(Mix[i, 2L]))*Mix[i, 3L]
  xdiff <- diff(xout[1:2])
  xoutmean <- sum(xout * yout) * xdiff
  xoutvar <- sum((xout - xoutmean)^2 * yout) * xdiff
  xoutcdf <- cumsum(yout) * xdiff
  if(is.na(e$loglik)) {
    xoutmedian <- NA_real_
  } else
    xoutmedian <- approx(xoutcdf, xout, xout=0.5, ties="ordered")$y
  stats <- as.data.frame(rbind(xstats, c(xoutmean, xoutvar, xoutmedian)))
  row.names(stats) <- c(xnam, "Computed")
  tuning <- c(Nmix=Nmix)
  # Fix LogLik
  loglik <- e$loglik
  attr(loglik, "nobs") <- NOBS
  if(is.na(loglik)) {
    attr(loglik, "df") <- NA_real_
  } else
    attr(loglik, "df") <- length(Mix) - 1
  class(loglik) <- "logLik"
  retval <- list(x=xout, y=yout, xnam=xnam, stats=stats, tuning=tuning, xin=x, 
                 call=call, method="Gaussian Mixture", Pzero=0, loglik=loglik,
                 mix=Mix, posterior=e$posterior)
  class(retval) <- "Density"
  # More goodness of fit scores
  GOF <- gofScores(x, retval)
  retval$ZaNorm <- GOF$ZaNorm
  retval$VarRat <- GOF$VarRat
  retval$VarDen <- GOF$VarDen
  return(retval)
}
