#' Goodness of Fit Scores
#' 
#' Computes the Zhang modified goodness of fit test type A, without computing the 
#'p-level. The value is normalized such that 0 indicates a perfect fit (ZaNorm). A 
#'convolution measure of the variability of the density estimate is 
#'also computed as VarDen. VarRat is the ratio of the variance of the fit to that 
#'of the data. Used primarily to asses the relative goodness of fits between 
#'several density estimates.
#' 
#' @details The convolution measure of the variability of the density estimate
#'(VarDen) is computed as the sum of the squared differences between 101 uniformly 
#'distributed points in the range of \code{x} (100 differences). The square-root
#'of the value divided by the peak value of density and is adjusted so that a
#'normal distribution gives a value of 1.0.
#' 
#' @param x the original sample data
#' @param plot.scores logical, if \code{TRUE}, then plot the ZaNorm and LogLik scores 
#'with a reference line for the kernel density estimate relation. If \code{FALSE},
#'the default, then no plot is created. The points are labelled if the Density
#'objects are named in the call.
#' @param \dots any number of objects of class "Density." To create abbreviations
#'for the Density objects, they can be named in the call.
#' 
#' @note Smaller values of Za indicate closer matches to the distribution.
#' 
#' @return a data frame of the ZaNorm, VarD, and VarRat scores, and the log-likelihood
#'(LogLik) and degrees of freedom (DF) with one row per 
#'density object. The data are sorted by ZaNorm.
#'
#' @importFrom stats bw.nrd0 density
#' @export
gofScores <- function(x, ..., plot.scores=FALSE) {
  dots <- list(...)
  dotnames <- names(dots)
  Labels <- TRUE
  if(is.null(dotnames)) {
    dotnames <- unlist(lapply(substitute(list(...))[-1], deparse))
    Labels <- FALSE
  }
  if(length(dots) == 0)
    stop("No Density objects provided.")
  x <- sort(x)
  n <- length(x)
  seqn <- seq(n)
  N <- length(dots)
  Za <- VarD <- VarRat <- LogLik <- DF <- numeric(N)
  Za0 <- 3.28978 - 0.67428/n
  VarSamp <- var(x)
  for(i in seq(N)) {
    p1 <- pDensity(x, dots[[i]])
    Za[i] <- -sum(log(p1)/(n - seqn + .5) + log(1 - p1)/(seqn - .5)) - Za0
    VarD[i] <- sum(diff(dDensity(seq(min(x), max(x), length=101),dots[[i]]))^2)
    VarD[i] <- sqrt(VarD[i])/max(dots[[i]]$y) * 5.044505
    VarRat[i] <- dots[[i]]$stats[nrow(dots[[i]]$stats), 2L]/VarSamp
    LogLik[i] <- as.numeric(dots[[i]]$loglik)
    DF[i] <- attr(dots[[i]]$loglik, "df")
  }
  retval <- data.frame(Name = dotnames, ZaNorm=Za, VarDen=VarD, VarRat=VarRat, 
                       LogLik=LogLik, DF=DF)
  retval <- retval[order(Za), ]
  row.names(retval) <- as.character(seq(N))
  if(plot.scores) {
    plot(retval$ZaNorm, retval$LogLik, log='x', xlab="Normalized Za Score",
         ylab="Log Likelihood", pch=16, cex=0.6)
    if(Labels)
      text(retval$ZaNorm, retval$LogLik, labels=dotnames, pos=4)
    # Compute the reference line if Pzero  is 0
    if(dots[[1]]$Pzero == 0) {
      bw1 <- bw.nrd0(x)
      Fs <- c(0.70, 0.72, 0.74, 0.77, 0.81, 0.87, 0.95, 1.05, 1.15, 1.25, 1.40, 
              1.60, 1.8, 2.0)
      Zas <- LLs <- numeric(length(Fs))
      for(i in seq(along=Fs)) {
        d1 <- density(x, bw=bw1*Fs[i])
        d1$Pzero <- 0
        p1 <- pDensity(x, d1)
        Zas[i] <- -sum(log(p1)/(n - seqn + .5) + log(1 - p1)/(seqn - .5)) - Za0
        LLs[i] <- sum(log(dDensity(x, d1)))
      }
      lines(Zas,LLs)
    }
    
  }
  return(retval)
}
