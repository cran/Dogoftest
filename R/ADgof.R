#' @title Anderson-Darling Goodness-of-Fit Test for a Specified Distribution
#'
#' @description
#' Performs the Anderson-Darling (AD) goodness-of-fit test for a given univariate distribution.
#' The function computes the AD statistic and returns an approximate p-value based on adjusted formulas.
#'
#' @param x A numeric vector of sample observations.
#' @param dist A character string specifying the null distribution. Options are
#'   \code{"norm"}, \code{"exp"}, \code{"unif"}, \code{"lnorm"}, \code{"weibull"},
#'   \code{"gamma"}, \code{"t"}, and \code{"chisq"}.
#' @param ... Additional named parameters passed to the corresponding distribution functions
#'   (e.g., \code{mean}, \code{sd}, \code{rate}, \code{df}, etc.).
#' @param eps A small positive constant to avoid log(0) during computation (default: \code{1e-15}).
#'
#' @return A list of class \code{"htest"} with components:
#' \describe{
#'   \item{statistic}{The value of the Anderson-Darling test statistic.}
#'   \item{p.value}{The approximate p-value computed using adjustment formulas.}
#'   \item{method}{A description of the test performed.}
#'   \item{data.name}{The name of the input data.}
#' }
#'
#' @details
#' This implementation supports several common distributions. Parameters of the null distribution
#' must be supplied via \code{...}. The p-value is calculated using the approximations suggested
#' by Stephens (1986) and other refinements. For small samples or custom distributions, a bootstrap
#' version may be preferred.
#'
#' @examples
#' set.seed(123)
#' x1 <- rnorm(500, mean = 5, sd = 2)
#' ADgof(x1, dist = "norm", mean = 5, sd = 2)
#'
#' x2 <- rexp(400, rate = 1.5)
#' ADgof(x2, dist = "exp")
#' ADgof(x2, dist = "exp", rate = 1.5)
#'
#' x3 <- runif(300, min = -2, max = 4)
#' ADgof(x3, dist = "unif", min = -2, max = 4)
#' @importFrom stats ecdf na.omit pchisq pexp pgamma plnorm pt punif pweibull quantile runif
#' @export

ADgof <- function(x,
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
  pfun <- pfuns[[dist]](...)
  x_s  <- sort(x)
  p    <- pfun(x_s)
  p[p < eps]       <- eps
  p[p > 1 - eps]   <- 1 - eps
  i    <- seq_len(n)
  Wp   <- (2 * i - 1) * (log(p) + log(1 - rev(p)))
  A    <- -n - (1 / n) * sum(Wp)
  AA   <- A * (1 + 0.75 / n + 2.25 / n^2)
  if (AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14  * AA - 223.73  * AA^2)
  } else if (AA < 0.34) {
    pval <- 1 - exp(-8.318  + 42.796  * AA - 59.938  * AA^2)
  } else if (AA < 0.6) {
    pval <-      exp( 0.9177 - 4.279   * AA - 1.38    * AA^2)
  } else if (AA < 100) {
    pval <-      exp( 1.2937 - 5.709   * AA + 0.0186  * AA^2)
  } else {
    pval <- exp(-1.733 * sqrt(AA))
  }
  RVAL <- list(
    statistic = c(AD = A),
    p.value   = pval,
    method    = paste("Anderson-Darling gof test for", dist, "distribution"),
    data.name = DNAME
  )
  class(RVAL) <- "htest"
  RVAL
}
