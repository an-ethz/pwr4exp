#' Determine the sample size required to achieve a target power
#'
#' This function finds the sample size required to achieve a target power for a given design quote.
#'
#' @param design.quote The design quote to use
#' @param alpha significance level, default 0.05
#' @param target.power The target power to achieve, default 0.8
#' @param n_init The initial sample size to start with, default 2
#' @param n_max The maximum sample size to try, default 99
#' @param ... Additional arguments to pass to pwr.anova
#'
#' @return A data frame with the alpha, power, and best sample size for each sample size tried, NA if the power is not achieved within the maximum sample size
#' @export
find_sample_size <- function(design.quote, alpha = 0.05, target.power = 0.8, n_init = 2, n_max = 99, ...){
  n = n_init
  res <- data.frame(pwr.anova(eval(design.quote), alpha = alpha, ...))
  res <- res[, c("alpha", "power")]
  best_n <- rep(NA, nrow(res))

  for (n in n_init:n_max) {
    power_result <- pwr.anova(eval(design.quote), alpha = alpha, ...)
    power_diff <- power_result$power - target.power
    res[is.na(best_n) & power_diff > 0, "power"] <- power_result[is.na(best_n) & power_diff > 0, "power"]
    best_n[is.na(best_n) & power_diff > 0] <- n
    if (all(power_diff > 0)) break
  }
  res$best_n <- best_n
  return(res)
}
