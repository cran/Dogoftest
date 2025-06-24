#' Perform a Simple Cramer-von Mises Goodness-of-Fit Test
#'
#' This function performs a simple Cramer-von Mises goodness-of-fit test to assess
#' whether a given sample comes from a uniform distribution. The test statistic and
#' p-value are calculated based on the sorted sample data.
#'
#' @param X A numeric vector containing the sample data.
#'
#' @return
#' \item{statistic}{The value of the Cramer-von Mises test statistic.}
#' \item{pvalue}{The p-value for the test.}
#' \item{statname}{The name of the test statistic.}
#'
#' @export
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- runif(100)  # Generate a sample from a uniform distribution
#' result <- simpleCvMgof(X)
#' print(result)
#'
#' # Example with non-uniform data:
#' Y <- rnorm(100)  # Generate a sample from a normal distribution
#' result <- simpleCvMgof(Y)
#' print(result)

simpleCvMgof <- function(X) {
  X <- sort(X)
  n <- length(X)
  k <- seq_len(n)
  W2 <- 1/(12 * n) + sum((X - (2*k - 1)/(2*n))^2)
  pvalue <- 1 - exp(-W2^2 / 2)
  return(list(statistic = W2, pvalue = pvalue, statname = "CvM"))
}
