#' @title Two-Sample Anderson-Darling Test (Bootstrap Version)
#'
#' @description
#' Performs a two-sample Anderson-Darling (AD) goodness-of-fit test using bootstrap resampling
#' to compare whether two samples come from the same distribution. This test is sensitive to
#' differences in both location and shape between the two distributions.
#'
#' @param x A numeric vector of data values from the first sample.
#' @param y A numeric vector of data values from the second sample.
#' @param alternative Character string specifying the alternative hypothesis.
#'        One of `"two.sided"` (default), `"less"`, or `"greater"`.
#' @param nboots Integer. Number of bootstrap replicates to compute the null distribution (default: 2000).
#' @param keep.boots Logical. If `TRUE`, returns the full vector of bootstrap statistics (default: `FALSE`).
#'
#' @return A list of class `"htest"` containing:
#' \describe{
#'   \item{statistic}{The observed Anderson-Darling test statistic.}
#'   \item{p.value}{The estimated bootstrap p-value.}
#'   \item{alternative}{The alternative hypothesis used.}
#'   \item{method}{A character string describing the test.}
#'   \item{bootstraps}{(Optional) A numeric vector of bootstrap statistics if `keep.boots = TRUE`.}
#' }
#'
#' @details
#' The test computes the Anderson-Darling statistic using the pooled empirical distribution
#' functions (ECDFs) of the two samples. A bootstrap procedure resamples the group labels
#' to approximate the null distribution and compute a p-value. If `p.value = 0`, it is adjusted
#' to `1 / (2 * nboots)` for stability.
#'
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100, mean = 0, sd = 4)
#' y <- rnorm(100, mean = 2, sd = 4)
#' AD2gof(x, y)
#'
#' @export
AD2gof <- function(x, y,
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
  ad_stat <- function(comb, labels) {
    ecdf.x <- ecdf(comb[labels])
    ecdf.y <- ecdf(comb[!labels])
    t <- comb
    Fxd <- ecdf.x(t)
    Fyd <- ecdf.y(t)
    diff2 <- switch(alternative,
                    two.sided = (Fxd - Fyd)^2,
                    greater   = pmax(Fxd - Fyd, 0)^2,
                    less      = pmax(Fyd - Fxd, 0)^2)
    Fmid <- (Fxd + Fyd) / 2
    weight <- 1 / (Fmid * (1 - Fmid))
    weight[!is.finite(weight)] <- 0
    dt   <- diff(t)
    mids <- (diff2[-length(diff2)] * weight[-length(weight)] +
               diff2[-1] * weight[-1]) / 2
    sum(mids * dt)
  }
  obs_stat <- ad_stat(comb, labels)
  boots <- if (keep.boots) numeric(nboots) else NULL
  count <- 0L
  for (i in seq_len(nboots)) {
    lab.b <- rep(FALSE, length(comb))
    lab.b[sample.int(length(comb), n.x)] <- TRUE
    bstat <- ad_stat(comb, lab.b)
    if (bstat >= obs_stat) count <- count + 1L
    if (keep.boots) boots[i] <- bstat
  }
  p.value <- count / nboots
  if (p.value == 0) p.value <- 1 / (2 * nboots)
  result <- list(
    statistic   = obs_stat,
    p.value     = p.value,
    alternative = alternative,
    method      = "Two samples Anderson-Darling test (Bootstrap)"
  )
  if (keep.boots) result$bootstraps <- boots
  class(result) <- "htest"
  return(result)
}
