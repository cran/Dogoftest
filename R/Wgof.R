#' Watson goodness-of-fit test
#' Performs the Watson test for goodness-of-fit to a specified distribution.
#'
#' @param x Numeric vector of observations.
#' @param dist Character string specifying the distribution to test against.
#'   One of \code{"norm"}, \code{"exp"}, \code{"unif"}, \code{"lnorm"}, or \code{"gamma"}.
#' @param ... Additional parameters passed to the distribution's cumulative distribution function (CDF).
#'   For example, \code{mean} and \code{sd} for the normal distribution.
#' @param eps Numeric tolerance for probability bounds to avoid extremes (default: 1e-15).
#'
#' @return An object of class \code{"htest"} containing the test statistic, p-value, method description, data name,
#'   and any distribution parameters used.
#'
#' @details
#' The Watson test is a modification of the Cramér–von Mises test, adjusting for mean deviations.
#' It measures the squared distance between the empirical distribution function of the data and the specified
#' theoretical cumulative distribution function, with a correction for location.
#'
#' @examples
#' set.seed(123)
#' x_norm <- rnorm(1000, mean = 5, sd = 2)
#' Wgof(x_norm, dist = "norm", mean = 5, sd = 2)
#'
#' x_exp <- rexp(500, rate = 0.5)
#' Wgof(x_exp, dist = "exp", rate = 0.5)
#'
#' x_unif <- runif(300, min = 0, max = 10)
#' Wgof(x_unif, dist = "unif", min = 0, max = 10)
#'
#' x_lnorm <- rlnorm(200, meanlog = 0, sdlog = 1)
#' Wgof(x_lnorm, dist = "lnorm", meanlog = 0, sdlog = 1)
#'
#' x_gamma <- rgamma(400, shape = 1, rate = 1)
#' Wgof(x_gamma, dist = "gamma", shape = 1, rate = 1)
#'
#' @export

Wgof <- function(x,
                 dist = c("norm", "exp", "unif", "lnorm", "gamma"),
                 ...,
                 eps = 1e-15) {
  dist <- match.arg(dist)
  DNAME <- deparse(substitute(x))
  x <- na.omit(x)
  n <- length(x)
  if (n < 1) stop("At least one observation is required.")
  params <- list(...)

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
  pfun <- pfuns[[dist]](...)
  x_s   <- sort(x)
  p     <- pfun(x_s)
  p[p < eps]       <- eps
  p[p > 1 - eps]   <- 1 - eps
  i     <- seq_len(n)
  W <- 1/(12 * n) + sum((p - (2 * i - 1)/(2 * n))^2)
  U <- W - n * (mean(p) - 0.5)^2
  UU <- (1 + 0.5 / n) * U
  if (UU < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * UU - 12542.61 * UU^2)
  } else if (UU < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * UU - 1515.29 * UU^2)
  } else if (UU < 0.092) {
    pval <- exp(0.886 - 31.62 * UU + 10.897 * UU^2)
  } else if (UU < 1.1) {
    pval <- exp(1.111 - 34.242 * UU + 12.832 * UU^2)
  } else {
    pval <- 7.37e-10
  }
  RVAL <- list(
    statistic = c(UU = UU),
    p.value = pval,
    method = paste("Watson test of goodness-of-fit for", dist, "distribution"),
    data.name = DNAME,
    parameter = if (length(params) > 0) unlist(params) else NULL
  )
  class(RVAL) <- "htest"
  return(RVAL)
}
