#' Extend lmerModLmerTest class
#'
#' This class extends \code{lmerModLmerTest} by adding a \code{DenDF} slot.
#' @importClassesFrom lmerTest lmerModLmerTest
#' @slot DenDF Numeric vector of denominator degrees of freedom.
#' @export
methods::setClass(
  Class = "customLmerMod",
  contains = "lmerModLmerTest",
  slots = c(DenDF = "numeric")
)

#' Calculate variance covariance parameters
#'
#' Scale variance-covariance matrices as the relative Cholesky factors of each random effect term.
#'
#' @param VarCov variance-covariance matrices. If there are multiple random effect groups,
#' supply the variance-covariance matrix of each group as an element in a list.
#' @param sigma standard deviation of random errors.
#' @return theta
#' @seealso "theta" in \code{\link[lme4]{getME}}, \code{\link[lme4]{lmer}}
calc.theta <- function(VarCov, sigma) {
  UseMethod('calc.theta', VarCov)
}

#' @exportS3Method calc.theta default
#' @noRd
calc.theta.default <- function(VarCov, sigma) {
  L <- chol(VarCov)
  theta <- L[upper.tri(L, diag = TRUE)] / sigma
  return(theta)
}

#' @exportS3Method calc.theta list
#' @noRd
calc.theta.list <- function(VarCov, sigma) {
  unname(unlist(
    lapply(
      VarCov,
      function(x) {
        calc.theta.default(x, sigma)
      }
    )
  ))
}

#' Naming theta
#' Naming the vector in the order of model specification and in the actual order used in the model
#' @param data data frame
#' @param formula model formula
theta.names <- function(data, formula) {
  # random effect terms in the order of fitting process
  lmod <- lme4::lFormula(formula, data)
  reTrms <- lmod$reTrms$cnms
  # random effect terms in the order of model formula
  grp <- as.character(sapply(lme4::findbars(formula), '[[', 3))
  reTrms0 <- reTrms[grp]
  # generate theta terms from random effects terms
  get_thetaTrms <- function(reTrms){
    thetaTrms <- c()
    for (group in names(reTrms)) {
      terms <- reTrms[[group]]
      n_terms <- length(terms)
      # 1. Add the variance term for intercept first
      if ("(Intercept)" %in% terms) {
        thetaTrms <- c(thetaTrms, paste0(group, ".(Intercept)"))
      }
      # 2. Add covariance terms if there is more than one random effect term (random slope + intercept)
      if (n_terms > 1) {
        for (term in terms) {
          if (term != "(Intercept)") {
            # Covariance between slope and intercept
            thetaTrms <- c(thetaTrms, paste0(group, ".", term, ".(Intercept)"))
          }
        }
      }
      # 3. Add the variance terms for slopes after the covariances
      for (term in terms) {
        if (term != "(Intercept)") {
          thetaTrms <- c(thetaTrms, paste0(group, ".", term))
        }
      }
    }
    return(thetaTrms)
  }
  thetaNames <- list(input = get_thetaTrms(reTrms0), output = get_thetaTrms(reTrms))
  return(thetaNames)
}

#' Create an artificial model object
#'
#' Create a pseudo-model object with the response variable being simulated
#' according to the fixed and random effects. Model coefficients are replaced
#' by the expectations specified in the argument \code{beta}. Variance-covariance
#' components of random effects are replaced by the values specified in argument
#' \code{VarCov}. The standard deviation of random error is replaced by the
#' argument \code{sigma}. Creating such a pseudo-model facilitates power
#' calculations by leveraging the \code{anova} function in \code{lmerTest} and
#' the \code{Anova} function in \code{car}.
#'
#' @param formula an object of class \code{\link{formula}}
#' @param data a data frame with the independent variables of the design as
#' columns, e.g., treatment factors and block factors.
#' @param beta a vector of the expectations of model coefficients.
#' @param VarCov variance-covariance matrices. If there are multiple random effect groups,
#' supply the variance-covariance matrix of each group as an element in a list.
#' @param sigma standard deviation of error
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of sum of squares, with Type III as the default, and the method for
#' computing the denominator degrees of freedom, with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#' @return a pseudo model object.
fit.pseu.model <- function(formula, data, beta, VarCov, sigma, ...) {
  has_random_effects <- any(grepl("\\(|\\)", as.character(formula)))
  if (has_random_effects) {
    return(fit.pseu.model.lmer(formula, data, beta, VarCov, sigma, ...))
  } else {
    return(fit.pseu.model.lm(formula, data, beta, sigma))
  }
}

#' @noRd
fit.pseu.model.lmer <- function(formula, data, beta, VarCov, sigma, ...) {
  # pseudo variable for construct model structure
  data$y <- 0
  # Calculate covariance parameters
  theta <- calc.theta(VarCov = VarCov, sigma = sigma)
  thetaNames <- theta.names(data = data, formula = formula)
  names(theta) <- thetaNames[["input"]]
  theta <- theta[thetaNames[["output"]]]
  # Internal function to fit a pseudo-model
  pseu.model <- function() {
    suppressMessages(
      data$y <- lme4::simulate.formula(
        object = formula[-2],
        nsim = 1,
        newparams = list(beta = beta, theta = theta, sigma = sigma),
        newdata = data
      )[[1]]
    )
    pseu_model <- try(lmerTest::lmer(formula, data))
    if (inherits(pseu_model, "try-error") || lme4::isSingular(pseu_model)) {
      return(NULL)
    } else {
      return(pseu_model)
    }
  }
  # Repeat fitting until a non-singular model is obtained
  repeat {
    suppressMessages(model <- pseu.model())
    if (!is.null(model)) break
  }
  # Modify model parameters
  DenDF <- round(stats::anova(model, ...)$DenDF, 0)
  model@beta <- beta
  model@theta <- theta
  model@devcomp$cmp[["sigmaREML"]] <- sigma
  model@sigma <- sigma
  model@pp <- lme4::merPredD$new(
    model@pp$X,
    model@pp$Zt,
    model@pp$Lambdat,
    model@pp$Lind,
    theta,
    n = nrow(model@pp$X)
  )
  model@vcov_beta <- as.matrix(lme4::vcov.merMod(model))
  model <- methods::new("customLmerMod", model, DenDF = DenDF)
  return(model)
}

#' @noRd
fit.pseu.model.lm <- function(formula, data, beta, sigma) {
  model_matrix <- stats::model.matrix(formula[-2], data)
  data$y <- model_matrix %*% beta + stats::rnorm(nrow(data), 0, sigma)
  pseu_model <- stats::lm(formula, data)
  # modify model parameters
  pseu_model$coefficients[] <- beta
  # expected RSS
  desired_rss <- (sigma^2) * pseu_model$df.residual
  # Adjust the residuals to achieve the expected RSS
  original_residuals <- stats::residuals(pseu_model)
  current_rss <- sum(original_residuals^2)
  adjustment_factor <- sqrt(desired_rss / current_rss)
  adjusted_residuals <- original_residuals * adjustment_factor
  # Create new fitted values using the new coefficients
  new_fitted_values <- as.vector(model_matrix %*% beta)
  # Update the pseu_model components manually
  pseu_model$fitted.values[] <- new_fitted_values
  pseu_model$residuals <- adjusted_residuals
  pseu_model$model$y[, 1] <- new_fitted_values + adjusted_residuals
  return(pseu_model)
}
