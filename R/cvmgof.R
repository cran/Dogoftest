#' Perform the Cramer-von Mises Goodness-of-Fit Test for Normality
#'
#' @param x A numeric vector containing the sample data.
#'
#' @return
#' \item{statistic}{The value of the Cramer-von Mises test statistic.}
#' \item{p.value}{The p-value for the test.}
#' \item{method}{A character string describing the test.}
#'
#' @export
#' @importFrom stats pnorm sd
#' @examples
#' # Example usage:
#' set.seed(123)
#' x <- rnorm(100)  # Generate a sample from a normal distribution
#' result <- cvmgof(x)
#' print(result)
#'
#' # Example with non-normal data:
#' y <- rexp(100)  # Generate a sample from an exponential distribution
#' result <- cvmgof(y)
#' print(result)

cvmgof <- function(x) {
  x <- sort(x) # 对数据进行排序
  n <- length(x) # 计算样本大小
  mean.x <- mean(x) # 计算样本均值和标准差
  sd.x <- sd(x)
  p <- pnorm((x - mean.x) / sd.x) # 计算正态分布的累积分布函数（CDF）值
  W <- (1/(12 * n) + sum((p - (2 * seq(1:n) - 1)/(2 * n))^2)) # 计算W统计量
  WW <- (1 + 0.5/n) * W # 调整W统计量
  if (WW < 0.0275) { # 根据WW的值计算p值
    pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
  } else if (WW < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
  } else if (WW < 0.092) {
    pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)
  } else if (WW < 1.1) {
    pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)
  } else {
    pval <- 7.37e-10
  }
  cvmt <- list(statistic = c(W = W), p.value = pval, method = "Cramer-von Mises normality test")
  return(cvmt)
}
