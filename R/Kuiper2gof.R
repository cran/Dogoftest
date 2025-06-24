#' @title Two-Sample Kuiper Test with Bootstrap
#'
#' @description
#' Performs a two-sample Kuiper test using bootstrap resampling to test whether two independent samples come from the same distribution.
#'
#' @param x,y Numeric vectors of data values for the two samples.
#' @param alternative Character string indicating the alternative hypothesis. Must be one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}.
#' @param nboots Integer. Number of bootstrap resamples to compute the empirical null distribution (default: 2000).
#' @param keep.boots Logical. If \code{TRUE}, returns all bootstrap test statistics.
#'
#' @return An object of class \code{"htest"} containing:
#' \describe{
#'   \item{statistic}{The observed Kuiper statistic.}
#'   \item{p.value}{The p-value computed from the bootstrap distribution.}
#'   \item{alternative}{The specified alternative hypothesis.}
#'   \item{method}{A character string describing the test.}
#'   \item{bootstraps}{(If requested) A numeric vector of bootstrap statistics.}
#' }
#'
#' @details
#' The Kuiper test is a nonparametric test similar to the Kolmogorov–Smirnov test,
#' but sensitive to discrepancies in both location and shape between two distributions.
#' This implementation uses bootstrap resampling to estimate the p-value.
#'
#' The \code{two.sided} test uses the sum of maximum positive and negative ECDF differences.
#' The \code{greater} and \code{less} options use one-sided variations.
#'
#' If the observed test statistic exceeds all bootstrap values, the p-value is set to \code{1 / (2 * nboots)} to avoid zero.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100, 0, 4)
#' y <- rnorm(100, 2, 4)
#' Kuiper2gof(x, y)
#'
#' @export

Kuiper2gof <- function(x, y,
                       alternative = c("two.sided", "less", "greater"),
                       nboots = 2000,
                       keep.boots = FALSE) {
  alternative <- match.arg(alternative)
  x <- na.omit(x); y <- na.omit(y)
  n.x <- length(x); n.y <- length(y)
  if (n.x < 1 || n.y < 1)
    stop("Samples x and y must each contain at least one observation")
  comb   <- c(x, y)
  labels <- c(rep(TRUE, n.x), rep(FALSE, n.y))
  ord    <- order(comb)
  comb   <- comb[ord]
  labels <- labels[ord]
  kuiper_stat <- function(comb, labels) {
    ecdf.x <- ecdf(comb[labels])
    ecdf.y <- ecdf(comb[!labels])
    t <- sort(unique(comb))
    Fxd <- ecdf.x(t)
    Fyd <- ecdf.y(t)

    Dplus  <- max(Fxd - Fyd)
    Dminus <- max(Fyd - Fxd)

    switch(alternative,
           two.sided = Dplus + Dminus,
           greater   = Dplus,
           less      = Dminus)
  }
  obs_stat <- kuiper_stat(comb, labels)
  boots <- if (keep.boots) numeric(nboots) else NULL
  count <- 0L
  for (i in seq_len(nboots)) {
    lab.b <- rep(FALSE, length(comb))
    lab.b[sample.int(length(comb), n.x)] <- TRUE
    bstat <- kuiper_stat(comb, lab.b)
    if (bstat >= obs_stat) count <- count + 1L
    if (keep.boots) boots[i] <- bstat
  }
  p.value <- count / nboots
  if (p.value == 0) p.value <- 1 / (2 * nboots)
  result <- list(
    statistic   = obs_stat,
    p.value     = p.value,
    alternative = alternative,
    method      = "Two samples Kuiper test (Bootstrap)"
  )
  if (keep.boots) result$bootstraps <- boots
  class(result) <- "htest"
  return(result)
}
