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

#' Calculate covariance parameter vector (theta)
#'
#' This function calculates the covariance parameters
#'
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma error standard deviation.
#' @return covariance parameter vector
#' @details
#' For more details on the structure and estimation of `theta`, refer to the documentation
#' of the \code{\link[lme4]{getME}} function in the \code{lme4} package.
#' @seealso \code{\link[lme4]{getME}}, \code{\link[lme4]{lmer}}
#' @export
calc.theta <- function(VarCov, sigma) {
  UseMethod('calc.theta', VarCov)
}

#' @exportS3Method calc.theta default
calc.theta.default <- function(VarCov, sigma) {
  L <- chol(VarCov)
  theta <- L[upper.tri(L, diag = TRUE)] / sigma
  return(theta)
}

#' @exportS3Method calc.theta list
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

#' Create an artificial model object
#'
#' Create a pseudo-model object with the response variable being simulated according to the fixed and random effects.
#' Model coefficients are replaced by the expectations specified in the argument \code{beta}.
#' Variance-covariance components of random effects are replaced by the values specified in argument \code{VarCov}.
#' The standard deviation of random error is replaced by the argument \code{sigma}.
#' Creating such a pseudo-model facilitates power calculations by leveraging the \code{anova} function in \code{lmerTest} and the \code{Anova} function in \code{car}.
#'
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model that would be used to test effects in post-experimental data analysis.
#' For details on model specification, see \code{\link{lm}} for fixed effects and \code{\link{lmer}} for random effects.
#' @param data a data frame with the independent variables of a design as columns, e.g., treatment factors and block factors.
#' @param beta a vector of the expectations of model coefficients.
#' The coefficients for categorical variables are the coefficients of dummy variables created using \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma standard deviation of error
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default, and the method for computing the denominator degrees of freedom, with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#' @return a pseudo model object.
#' @export
fit.pseu.model <- function(formula, data, beta, VarCov, sigma, ...) {
  # Check if the formula contains random effects terms
  has_random_effects <- any(grepl("\\(|\\)", as.character(formula)))

  if (has_random_effects) {
    return(fit.pseu.model.lmer(formula, data, beta, VarCov, sigma, ...))
  } else {
    return(fit.pseu.model.lm(formula, data, beta, sigma))
  }
}

#' @noRd
fit.pseu.model.lmer <- function(formula, data, beta, VarCov, sigma, ...) {
  # Fill up missing beta values with zeros
  fixed_formula <- stats::reformulate(
    attr(stats::terms(formula), "term.labels")[!grepl("\\|", attr(stats::terms(formula), "term.labels"))],
    response = "y"
  )
  model_matrix <- stats::model.matrix(fixed_formula[-2], data)
  fixed_terms <- dimnames(model_matrix)[[2]]
  if (length(beta) != length(fixed_terms)) {
    beta <- c(beta, rep(0, length(fixed_terms) - length(beta)))
  }

  # Calculate covariance parameters
  theta <- calc.theta(VarCov = VarCov, sigma = sigma)

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
  terms <- dimnames(model_matrix)[[2]]

  if (length(beta) < length(terms)) {
    beta <- c(beta, rep(0, length(terms) - length(beta)))
  }

  data$y <- model_matrix %*% beta + stats::rnorm(nrow(data), 0, sigma)
  pseu_model <- stats::lm(formula, data)

  # Modify the pseu_model coefficients
  pseu_model$coefficients[] <- beta

  # Calculate the desired RSS based on the desired sigma
  desired_rss <- (sigma^2) * pseu_model$df.residual

  # Adjust the residuals to achieve the desired RSS
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
