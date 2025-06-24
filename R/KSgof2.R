#' One-sample Kolmogorov-Smirnov goodness-of-fit test
#'
#' Performs the one-sample Kolmogorov-Smirnov test for a specified theoretical distribution.
#'
#' @param x Numeric vector of observations.
#' @param dist Character string specifying the distribution to test against.
#'   One of \code{"norm"}, \code{"exp"}, \code{"unif"}, \code{"lnorm"}, \code{"weibull"},
#'   \code{"gamma"}, \code{"t"}, or \code{"chisq"}.
#' @param ... Additional parameters passed to the distribution’s cumulative distribution function (CDF).
#'   For example, \code{mean} and \code{sd} for the normal distribution.
#' @param eps Numeric lower and upper bound for tail probabilities to avoid numerical issues (default: \code{1e-15}).
#'
#' @return An object of class \code{"htest"} containing the test statistic, p-value, method description, and data name.
#'
#' @details
#' The test compares the empirical distribution function of \code{x} with the cumulative distribution function
#' of a specified theoretical distribution using the Kolmogorov-Smirnov statistic.
#' For large sample sizes, a p-value approximation based on the asymptotic distribution is used.
#'
#' A correction is applied when sample size exceeds 100, adjusting the test statistic to approximate a fixed sample size.
#' For very small or very large statistics, piecewise polynomial approximations are used to compute the p-value.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(1000, mean = 5, sd = 2)
#' KSgof2(x, dist = "norm", mean = 5, sd = 2)
#'
#' y <- rexp(500, rate = 0.5)
#' KSgof2(y, dist = "exp", rate = 0.5)
#'
#' u <- runif(300, min = 0, max = 10)
#' KSgof2(u, dist = "unif", min = 0, max = 10)
#'
#' @export

KSgof2 <- function(x,
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
  x_s   <- sort(x)
  p     <- pfun(x_s)
  p[p < eps]       <- eps
  p[p > 1 - eps]   <- 1 - eps
  i     <- seq_len(n)
  Dplus  <- max(i/n - p)
  Dminus <- max(p - (i - 1)/n)
  K      <- max(Dplus, Dminus)
  if (n <= 100) {
    Kd <- K; nd <- n
  } else {
    Kd <- K * (n/100)^0.49; nd <- 100
  }
  pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) +
                  2.99587 * Kd * sqrt(nd + 2.78019) -
                  0.122119 + 0.974598/sqrt(nd) + 1.67997/nd)

  if (pvalue > 0.1) {
    KK <- (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
    pvalue <- if      (KK <= 0.302) 1
    else if (KK <= 0.5)   2.76773 - 19.828315*KK + 80.709644*KK^2 -
      138.55152*KK^3 + 81.218052*KK^4
    else if (KK <= 0.9)  -4.901232 + 40.662806*KK - 97.490286*KK^2 +
      94.029866*KK^3 - 32.355711*KK^4
    else if (KK <= 1.31)  6.198765 - 19.558097*KK + 23.186922*KK^2 -
      12.234627*KK^3 + 2.423045*KK^4
    else                 0
  }
  RVAL <- list(
    statistic = c(D = (if (exists("KK")) KK else K)),
    p.value   = pvalue,
    method    = paste("Kolmogorov-Smirnov gof test for", dist, "distribution"),
    data.name = DNAME
  )
  class(RVAL) <- "htest"
  RVAL
}
