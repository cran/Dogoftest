#' @title Two-Sample Cramér–von Mises Test (Bootstrap Version)
#'
#' @description
#' Performs a nonparametric two-sample Cramér–von Mises test using a permutation-based
#' bootstrap method to assess whether two samples come from the same distribution.
#'
#' @param x Numeric vector of observations from the first sample.
#' @param y Numeric vector of observations from the second sample.
#' @param alternative Character string specifying the alternative hypothesis.
#'        Must be one of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param nboots Number of bootstrap replicates to approximate the null distribution (default: 2000).
#' @param keep.boots Logical. If \code{TRUE}, the bootstrap statistics will be returned (default: \code{FALSE}).
#'
#' @return An object of class \code{"htest"} with elements:
#' \describe{
#'   \item{statistic}{Observed Cramér–von Mises test statistic.}
#'   \item{p.value}{Bootstrap-based p-value.}
#'   \item{alternative}{The alternative hypothesis used.}
#'   \item{method}{A description of the test.}
#'   \item{bootstraps}{(Optional) Vector of bootstrap test statistics if \code{keep.boots = TRUE}.}
#' }
#'
#' @details
#' The test compares two empirical cumulative distribution functions (ECDFs).
#' The bootstrap procedure permutes group labels to generate the null distribution.
#' Tailored one-sided tests use one-sided squared differences of ECDFs.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100, mean = 0, sd = 4)
#' y <- rnorm(100, mean = 2, sd = 4)
#' CVM2gof(x, y)
#'
#' # One-sided test
#' CVM2gof(x, y, alternative = "greater")
#'
#' # Store bootstrap replicates
#' res <- CVM2gof(x, y, keep.boots = TRUE)
#' hist(res$bootstraps, main = "Bootstrap Distribution", xlab = "Test Statistic")
#'
#' @export

CVM2gof <- function(x, y,
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
  cvm_stat <- function(comb, labels) {
    ecdf.x <- ecdf(comb[labels])
    ecdf.y <- ecdf(comb[!labels])
    t    <- comb
    Fxd  <- ecdf.x(t)
    Fyd  <- ecdf.y(t)
    diff2 <- switch(alternative,
                    two.sided = (Fxd - Fyd)^2,
                    greater   = pmax(Fxd - Fyd, 0)^2,
                    less      = pmax(Fyd - Fxd, 0)^2)

    sum(diff2)
  }
  obs_stat <- cvm_stat(comb, labels)
  boots <- if (keep.boots) numeric(nboots) else NULL
  count <- 0L
  for (i in seq_len(nboots)) {
    lab.b <- rep(FALSE, length(comb))
    lab.b[sample.int(length(comb), n.x)] <- TRUE
    bstat <- cvm_stat(comb, lab.b)
    if (bstat >= obs_stat) count <- count + 1L
    if (keep.boots) boots[i] <- bstat
  }
  p.value <- count / nboots
  if (p.value == 0) p.value <- 1 / (2 * nboots)
  result <- list(
    statistic   = obs_stat,
    p.value     = p.value,
    alternative = alternative,
    method      = "Two samples Cramer von Mises test (Bootstrap)"
  )
  if (keep.boots) result$bootstraps <- boots
  class(result) <- "htest"
  result
}
