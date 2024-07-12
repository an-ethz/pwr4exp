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
#' This function calculates the covariance parameter vector
#'
#' @param VarCov variance-covariance components of random effects. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma error standard deviation
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
#' Create a pseudo model object with the same structure as the model object fitted using the exact data. Model coefficients are replaced by the effect size specified in the argument \code{beta}. The variance-covariance matrix of random effects is replaced by the argument \code{VarCov}. The standard deviation of error is replaced by the argument \code{sigma}.
#'
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param data a data frame with the independent variables of a design, e.g., treatment factors, block factors.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma standard deviation of error
#' @param ... other arguments passed to \code{\link{anova}}.
#' @return a pseudo model object \code{lmerMod} or \code{\link{lm}}.
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
  # fill up missing beta values with zeros
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

  # modify model parameters
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

#' Power calculation for main effects and interactions
#'
#' Calculate power for testing overall effects of factors and interactions in the form of ANOVA table.
#'
#' @param design a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @param alpha significance level, type I error. default 0.05
#' @param ...	other arguments passed to anova. By default, type III sum of squares and Satterthwaite’s approximation of denominator degrees of freedom are applied.
#' @return a data frame
#' @export
#' @examples
#' mydesign = designCRD(treatments = c(2, 2),
#'                      beta = c(10, 9, 8, 7),
#'                      sigma2 = 9,
#'                      replications = 10)
#' pwr.anova(mydesign)
#'
#' mydesign = designRCB(treatments = c(2, 2),
#'                       blocks = 10,
#'                       beta = c(10, 9, 8, 7),
#'                       VarCov = 8,
#'                       sigma2 = 9)
#' pwr.anova(mydesign)
pwr.anova <- function(design, alpha, ...) {
  UseMethod("pwr.anova", design)
}

#' @exportS3Method pwr.anova lmmDesign
pwr.anova.lmmDesign <- function(design, alpha = 0.05, ...){

  args <- list(...)
  if (!("type" %in% names(args))) {
    args$type <- "3"
  }

  if (!("ddf" %in% names(args))) {
    args$ddf <- "Satterthwaite"
  }

  res <- as.data.frame(stats::anova(design$pseu_model, ...))
  res$DenDF <- design$pseu_model@DenDF
  res$`non-centrality` <- res$`F value` * res$NumDF
  res$Fcrt <- stats::qf(p = alpha, df1 = res$NumDF, df2 = res$DenDF, lower.tail = FALSE)
  res$power <- 1 - stats::pf(q = res$Fcrt, df1 = res$NumDF, df2 = res$DenDF, ncp = res$`non-centrality`)
  res <- res[, c("NumDF", "DenDF", "non-centrality", "power")]
  res$alpha <- alpha
  res <- res[, c("NumDF", "DenDF", "non-centrality", "alpha", "power")]
  attr(res, "heading") <- paste("Power analysis of", design$design, "design")
  class(res) <- c("anova", "data.frame")
  return(res)
}

#' @exportS3Method pwr.anova lmDesign
pwr.anova.lmDesign <- function(design, alpha = 0.05, ...){
  formula <- design$formula
  # change contrast to contr.sum for type 3 anova
  data2 <- design$pseu_model$model
  for (factor_name in names(data2)) {
    if (is.factor(data2[[factor_name]])) {
      stats::contrasts(data2[[factor_name]]) <- stats::contr.sum
    }
  }
  design$pseu_model <- stats::update(design$pseu_model, data = data2)
  args <- list(...)
  if (!("type" %in% names(args))) {
    args$type <- "3"
  }
  res <- as.data.frame(car::Anova(design$pseu_model, ...))
  res$NumDF <- res$Df
  res$DenDF <- design$pseu_model$df.residual
  res$`non-centrality` <- res$`F value` * res$NumDF
  res$Fcrt <- stats::qf(p = alpha, df1 = res$NumDF, df2 = res$DenDF, lower.tail = FALSE)
  res$power <- 1 - stats::pf(q = res$Fcrt, df1 = res$NumDF, df2 = res$DenDF, ncp = res$`non-centrality`)
  res$alpha <- alpha
  res <- res[, c("NumDF", "DenDF", "non-centrality", "alpha", "power")]
  res <- res[!rownames(res) %in% c("(Intercept)", "Residuals"), ]
  attr(res, "heading") <- paste("Power analysis of", design$design, "design")
  class(res) <- c("anova", "data.frame")
  return(res)
}

#' Power calculation for contrasts
#'
#' Calculate power for testing various contrasts of fixed effects. The same syntax of \code{\link{emmeans}} package is employed to specify contrast types.
#'
#' @param design a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @param specs see \code{\link[emmeans]{emmeans}}
#' @param contrast see \code{\link[emmeans]{contrast}}
#' @param alpha significance level, type I error. default 0.05
#' @param ... other arguments passed to \code{\link[emmeans]{contrast}}. By default, kenward-roger approximation for degrees of freedom is applied.
#' @return a data frame
#' @export
#' @examples
#' mydesign = designRCB(treatments = c(2, 2),
#'                      blocks = 10,
#'                      beta = c(10, 9, 8, 7),
#'                      VarCov = 8, sigma2 = 9)
#' pwr.contrast(mydesign, ~ facA|facB, contrast = "pairwise")
pwr.contrast <- function(design, specs, contrast, alpha = 0.05, ...){
  args <- list(...)
  if (!("lmer.df" %in% names(args))) {
    args$lmer.df <- "kenward-roger"
  }
  if (inherits(design, "lmmDesign")) {
    emm <- do.call(
      emmeans::emmeans,
      c(list(object = design$pseu_model, specs = specs), args)
    )
  } else if (inherits(design, "lmDesign")) {
    emm <- do.call(
      emmeans::emmeans,
      c(list(object = design$pseu_model, specs = specs), args)
    )
  } else {
    stop("design is neither an lmmDesign nor an lmDesign object.")
  }
  emm.contrast <- as.data.frame(emmeans::contrast(emm, method = contrast))
  emm.contrast$f.ratio <- emm.contrast$t.ratio^2
  emm.contrast$df <- round(emm.contrast$df)
  emm.contrast$`non-centrality` <- emm.contrast$f.ratio
  emm.contrast$Fcrt <- stats::qf(p = alpha, df1 = 1, df2 = emm.contrast$df, lower.tail = FALSE)
  emm.contrast$alpha <- alpha
  emm.contrast$power <- 1 - stats::pf(
    q = emm.contrast$Fcrt,
    df1 = 1,
    df2 = emm.contrast$df,
    ncp = emm.contrast$`non-centrality`
  )
  emm.contrast <- emm.contrast[, !names(emm.contrast) %in% c("SE", "t.ratio", "p.value", "f.ratio", "Fcrt")]
  return(emm.contrast)
}
