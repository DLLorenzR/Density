#' Censored Estimation
#' 
#' Internal function  to compute the regression on order statistics (ROS)-
#'estimation values for censored values.
#'
#' @param Obs the reported values. Missing values are not permitted.
#' @param Cens the left-censor codes, if \code{TRUE}, then the corresponding
#'value in \code{Obs} is left-censored, otherwise the corresponding value is
#'the detected value. Missing values are not permitted.
#' @param dens a list containing x (ordinate) and y (density) values 
#'describing the distribution.
#' @param alpha the plotting position offset factor, default is 0.4.
#'
#' @note ROS-substitution works on single- or multiple-detection limit data. 
#'
#' @importFrom survival Surv
#' @importFrom stats lm predict
#' 
#' @return A vector like \code{Obs} with the left-censored values replaced with
#'the estimated values.
#'
rosEst <- function (Obs, Cens, dens, alpha = 0.4) {
  # Plotting points for censored data function
  censpp <- function (x, a = 0.4) {
    x <- sort(x)
    x.ndflag <- !x[,2L] # for now, only left-censored
    dl <- c(-Inf, unique(x[x.ndflag, 1L]), Inf)
    m <- length(dl)
    aj <- bj <- cj <- numeric(m)
    xunc <- x[, 1L][!x.ndflag]
    xcen <- x[, 1L][x.ndflag]
    for (j in 1:m) {
      aj[j] <- sum(dl[j] <= xunc & xunc < dl[j + 1])
      bj[j] <- sum(xunc < dl[j]) + sum(xcen <= dl[j])
      cj[j] <- sum(xcen == dl[j])
    }
    if (m > 3L) {
      inew <- 1L
      for (i in 2L:(m - 1L)) {
        if (aj[i] == 0) {
          cj[i + 1L] <- cj[i] + cj[i + 1L]
        }
        else {
          inew <- inew + 1L
          aj[inew] <- aj[i]
          bj[inew] <- bj[i]
          cj[inew] <- cj[i]
          dl[inew] <- dl[i]
        }
      }
      inew <- inew + 1L
      aj[inew] <- aj[m]
      bj[inew] <- bj[m]
      cj[inew] <- cj[m]
      dl[inew] <- dl[m]
      m <- inew
      length(aj) <- length(bj) <- length(cj) <- length(dl) <- m
    }
    if (sum(aj) != sum(!x.ndflag)) 
      warning("sum(aj) != # of hits = sum(!x.ndflag)\n  sum(aj) =", 
              sum(aj), "\n  sum(!x.ndflag) =", sum(!x.ndflag), 
              "\n")
    if (bj[m] != length(x)) 
      warning("bj[m] =", bj[m], "!= # of obs =", length(x), 
              "\n")
    pe <- c(1, numeric(m))
    for (j in m:2L) {
      pe[j] <- pe[j + 1L] + (1 - pe[j + 1L]) * aj[j]/(aj[j] +  bj[j])
    }
    ppunc <- xunc * 0
    ppcen <- numeric(length(xcen))
    isum <- 0
    jsum <- 0
    for (i in seq(m)) {
      if (aj[i] > 0) {
        for (j in seq(aj[i])) {
          ppt <- (j - a)/(aj[i] + 1 - 2 * a)
          ppunc[isum + j] <- (1 - pe[i]) + (pe[i] - pe[i + 
                                                         1]) * ppt
        }
        isum <- isum + aj[i]
      }
      if (cj[i] > 0) {
        for (j in seq(cj[i])) {
          ppt <- (j - a)/(cj[i] + 1 - 2 * a)
          ppcen[jsum + j] <- (1 - pe[i]) * ppt
        }
        jsum <- jsum + cj[i]
      }
    }
    list(x = xunc, pp = ppunc, xcen = xcen, ppcen = ppcen)
  }
  # End of plotting point function
  # Begin execution
  tmp.x <- Surv(Obs, event=1 - Cens, type="left")
  step1 <- censpp(tmp.x, a = alpha)
  qden <- qDensity(step1$pp, dens)
  step2 <- with(step1, lm(x ~ qden))
  # Step3 are the censored estimates
  qden <- qDensity(step1$ppcen, dens)
  step3 <- predict(step2, newdata = data.frame(qden=qden))
  # Arbitrary replacement of censored values
  Obs[Cens == 1] <- step3
  return(Obs)
}
