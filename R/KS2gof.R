#' @title Two-Sample Kolmogorov–Smirnov Test with Bootstrap
#'
#' @description
#' Performs a two-sample Kolmogorov–Smirnov (KS) test using a bootstrap method
#' to assess whether two independent samples come from the same distribution.
#'
#' @param x,y Numeric vectors of data values for the two independent samples.
#' @param alternative Character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}.
#' @param nboots Number of bootstrap resamples used to approximate the null distribution (default: 5000).
#' @param keep.boots Logical; if \code{TRUE}, return the vector of bootstrap statistics.
#'
#' @return An object of class \code{"htest"} with the following components:
#' \describe{
#'   \item{statistic}{The observed KS statistic.}
#'   \item{p.value}{The p-value based on the bootstrap distribution.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{method}{Description of the test used.}
#' }
#'
#' @details
#' This implementation performs a nonparametric KS test for equality of distributions
#' by resampling under the null hypothesis. It supports one-sided and two-sided alternatives.
#'
#' If \code{keep.boots = TRUE}, the function returns all bootstrap statistics,
#' which can be used for further analysis (e.g., plotting).
#'
#' If the p-value is zero due to no bootstrap statistic exceeding the observed value,
#' it is adjusted to \code{1 / (2 * nboots)} to avoid a zero p-value.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100, mean = 0, sd = 4)
#' y <- rnorm(100, mean = 2, sd = 4)
#' KS2gof(x, y)
#'
#' @export

KS2gof <- function(x, y,
                   alternative = c("two.sided", "less", "greater"),
                   nboots = 5000,
                   keep.boots = FALSE) {
  alternative <- match.arg(alternative)
  x <- na.omit(x); y <- na.omit(y)
  n.x <- length(x); n.y <- length(y)
  if (n.x < 1 || n.y < 1)
    stop("Samples x and y must each contain at least one observation")
  combined   <- c(x, y)
  labels     <- c(rep(TRUE, n.x), rep(FALSE, n.y))
  ord        <- order(combined)
  combined   <- combined[ord]
  labels     <- labels[ord]
  ks_stat <- function(values, labels) {
    ecdf.x <- ecdf(values[labels])
    ecdf.y <- ecdf(values[!labels])
    t      <- sort(unique(values))
    Fxd    <- ecdf.x(t)
    Fyd    <- ecdf.y(t)
    switch(alternative,
           two.sided = max(abs(Fxd - Fyd)),
           greater   = max(Fxd - Fyd),
           less      = max(Fyd - Fxd))
  }
  obs_stat <- ks_stat(combined, labels)
  boots    <- if (keep.boots) numeric(nboots) else NULL
  count    <- 0L
  for (i in seq_len(nboots)) {
    lab.b <- rep(FALSE, length(combined))
    lab.b[sample.int(length(combined), n.x)] <- TRUE
    bstat <- ks_stat(combined, lab.b)
    if (bstat >= obs_stat) count <- count + 1L
    if (keep.boots) boots[i] <- bstat
  }
  p.value <- count / nboots
  if (p.value == 0) p.value <- 1 / (2 * nboots)
  result <- list(
    statistic   = obs_stat,
    p.value     = p.value,
    alternative = alternative,
    method      = "Two samples Kolmogorov-Smirnov test (Bootstrap)"
  )
  if (keep.boots) result$bootstraps <- boots
  class(result) <- "htest"
  return(result)
}
