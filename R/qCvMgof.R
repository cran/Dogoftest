#' Calculate the Quantile of the Cramer-von Mises Goodness-of-Fit Statistic
#'
#' This function calculates the quantile of the Cramer-von Mises goodness-of-fit
#' statistic using the `uniroot` function to find the root of the given function.
#'
#' @param X A numeric vector containing the sample data.
#' @param p A numeric value representing the desired quantile probability.
#'
#' @return
#' \item{root}{The quantile value corresponding to the given probability.}
#'
#' @export
#' @importFrom stats uniroot
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- rnorm(100)  # Generate a sample from a normal distribution
#' p <- 0.95        # Desired quantile probability
#' result <- qCvMgof(X, p)
#' print(result)

qCvMgof <- function(X, p) {
  X <- sort(X)
  n <- length(X)
  k <- seq_len(n)
  W2 <- 1/(12 * n) + sum((X - (2*k - 1)/(2*n))^2)
  f <- function(q) {1 - exp(-q^2 / 2) - p}
  res <- uniroot(f, c(0, 10))
  return(res$root)
}
