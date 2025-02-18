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

qCvMgof <- function(X, p) { # 辅助函数，用于uniroot寻找分位数
  X <- sort(X)
  n <- length(X)      # 获取样本大小
  k <- seq_len(n)     # 创建一个从1到n的序列
  W2 <- 1/(12 * n) + sum((X - (2*k - 1)/(2*n))^2) # 计算Cramer-von Mises统计量
  f <- function(q) {1 - exp(-q^2 / 2) - p}
  res <- uniroot(f, c(0, 10)) # 使用uniroot函数找到分位数
  return(res$root)
}
