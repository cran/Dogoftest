#' @title One-Sample Cramér–von Mises Goodness-of-Fit Test
#'
#' @description
#' Performs the one-sample Cramér–von Mises goodness-of-fit (GoF) test to assess whether
#' a sample comes from a specified distribution using asymptotic p-value approximations.
#'
#' @param x A numeric vector of observations.
#' @param dist A character string specifying the theoretical distribution. Must be one of
#'   \code{"norm"}, \code{"exp"}, \code{"unif"}, \code{"lnorm"}, \code{"weibull"},
#'   \code{"gamma"}, \code{"t"}, or \code{"chisq"}.
#' @param ... Distribution parameters passed to the corresponding \code{p} functions
#'   (e.g., \code{mean}, \code{sd}, \code{rate}, \code{df}, etc.).
#' @param eps A small value to truncate extreme p-values (default is \code{1e-15}).
#'
#' @return An object of class \code{"htest"} with the following components:
#' \describe{
#'   \item{statistic}{The computed Cramér–von Mises test statistic.}
#'   \item{p.value}{The asymptotic p-value.}
#'   \item{method}{A description of the test and distribution.}
#'   \item{data.name}{The name of the data vector.}
#' }
#'
#' @details
#' The test uses the Cramér–von Mises statistic to assess how well the empirical
#' distribution function (EDF) of the sample agrees with the cumulative distribution
#' function (CDF) of the specified theoretical distribution. The p-value is computed
#' using approximation formulas derived from the asymptotic distribution of the test statistic.
#'
#' @examples
#' set.seed(123)
#' x1 <- rnorm(500, mean = 0, sd = 1)
#' CVMgof2(x1, dist = "norm", mean = 0, sd = 1)
#'
#' x2 <- rexp(500, rate = 2)
#' CVMgof2(x2, dist = "exp", rate = 2)
#'
#' x3 <- runif(200, min = -1, max = 3)
#' CVMgof2(x3, dist = "unif", min = -1, max = 3)
#'
#' @export

CVMgof2 <- function(x,
                   dist = c("norm", "exp", "unif", "lnorm", "weibull", "gamma", "t", "chisq"),
                   ...,
                   eps = 1e-15) {
  dist   <- match.arg(dist)
  DNAME  <- deparse(substitute(x))
  x      <- na.omit(x)
  n      <- length(x)
  if (n < 1) stop("At least one observation is required.")
  pfuns <- list(
    norm    = function(...) function(u) pnorm((u - list(...)[["mean"]]) / list(...)[["sd"]]),
    exp     = function(...) function(u) pexp(u, ...),
    unif    = function(...) function(u) punif(u, ...),
    lnorm   = function(...) function(u) plnorm(u, ...),
    weibull = function(...) function(u) pweibull(u, ...),
    gamma   = function(...) function(u) pgamma(u, ...),
    t       = function(...) function(u) pt(u, ...),
    chisq   = function(...) function(u) pchisq(u, ...)
  )
  pfun   <- pfuns[[dist]](...)
  x_s    <- sort(x)
  p      <- pfun(x_s)
  p[p < eps]       <- eps
  p[p > 1 - eps]   <- 1 - eps
  i      <- seq_len(n)
  W      <- (1 / (12 * n)) + sum((p - (2 * i - 1) / (2 * n))^2)
  WW     <- (1 + 0.5 / n) * W
  if (WW < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
  } else if (WW < 0.051) {
    pval <- 1 - exp(-5.903  + 179.546 * WW - 1515.29  * WW^2)
  } else if (WW < 0.092) {
    pval <-      exp( 0.886  - 31.62  * WW +   10.897 * WW^2)
  } else if (WW < 1.1) {
    pval <-      exp( 1.111  - 34.242 * WW +   12.832 * WW^2)
  } else {
    pval <- 7.37e-10
  }
  RVAL <- list(
    statistic = c(W = WW),
    p.value   = pval,
    method    = paste("Cramer von Mises gof test for", dist, "distribution"),
    data.name = DNAME
  )
  class(RVAL) <- "htest"
  RVAL
}
