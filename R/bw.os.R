#' Bandwidth for Kernel Density Estimation
#' 
#' Computes the bandwidth from the variance of x, specifically the oversmoothed 
#'bandwidth selector of Wand and Jones (1995, page 61) for \code{bw.os} or
#'the modified bandwidth for \code{density.lf} (\code{bw.lf}).
#' 
#' @aliases bw.os bw.lf
#' 
#' @param x numeric vector of observations from the distribution whose density 
#'is to be estimated.
#' @param kernel character string indicating the type of kernel. Must be one of
#'"normal," "box," "epanech," "biweight," or "triweight." The default is
#'"normal." See \bold{Details}.
#'
#' @details The value for \code{kernel} is required for functions that do not
#'automatically adjust the bandwidth for the variance of the kernel, like 
#'\code{bkde}. The function \code{density}, for example, does and the kernel
#'should always be "normal." 
#' 
#' @return A single numeric value representing the value.
#' 
#' @rdname bw.os
#' @export
bw.os <- function(x, kernel="normal") {
  kernel <- match.arg(kernel, c("normal", "box", "epanech", 
                                "biweight", "triweight"))
  del0 <- switch(kernel, normal = (1/(4 * pi))^(1/10), box = (9/2)^(1/5), 
                 epanech = 15^(1/5), biweight = 35^(1/5), triweight = (9450/143)^(1/5))
  bandwidth <- del0 * (243/(35 * length(x)))^(1/5) * sd(x)
  return(bandwidth)
}

#' @rdname bw.os
#' @export
bw.lf <- function(x) {
  return(diff(range(x))/(logb(length(x), base = 2) + 1)/0.5)
}
