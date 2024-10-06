#' Determine the sample size required to achieve the target power
#'
#' This function finds the minimum sample size needed to achieve the target power
#' for a given design. It uses an iterative approach to determine the minimum
#' number of replications by traversing through a series of integers.
#'
#' @param design.quote a quoted design object with unknown and unevaluated
#' replications to be evaluated with varying values
#' @param alpha type I error rate, default is 0.05
#' @param target.power the target power can be a single value for all factors or
#' a vector of containing individual values for different factors, default is 0.8
#' @param n_init the initial replications for the iterative process, default is 2
#' @param n_max the maximum number of replications for the iterative process, default is 99
#' @param ... additional arguments passed to \code{\link{pwr.anova}}
#'
#' @return A data frame with type I error rate (alpha), realized power (power),
#' and minimum sample size (best_n).
#' @examples
#' # create a LSD object with unknown replications (\code{squares = n})
#' # simply \code{\link{quote}} the design generating function with
#' lsd_quote <- quote(
#'   designLSD(
#'     treatments = 4,
#'     squares = n,
#'     reuse = "row",
#'     beta = c(10, 2, 3, 4),
#'     VarCov = list(5, 2),
#'     sigma2 = 10
#'   )
#' )
#'
#' # find the minimum number of squares required to achieve the target power of 0.8
#' find_sample_size(lsd_quote)
#' @export
find_sample_size <- function(design.quote, alpha = 0.05, target.power = 0.8, n_init = 2, n_max = 99, ...){
  n = n_init
  res <- data.frame(pwr.anova(eval(design.quote), alpha = alpha, ...))
  res <- res[, c("alpha", "power")]
  res$power <- NA
  best_n <- rep(NA, nrow(res))
  for (n in n_init:n_max) {
    power_result <- pwr.anova(eval(design.quote), alpha = alpha, ...)
    power_diff <- power_result$power - target.power
    res[is.na(best_n) & power_diff > 0, "power"] <- power_result[is.na(best_n) & power_diff > 0, "power"]
    best_n[is.na(best_n) & power_diff > 0] <- n
    res[is.na(best_n), "power"] <- power_result[is.na(best_n), "power"]
    if (all(power_diff > 0)) break
  }
  best_n[is.na(best_n)] <- n_max
  res$best_n <- best_n
  return(res)
}
