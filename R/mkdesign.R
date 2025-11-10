###############################################################################
# This file is part of the pwr4exp package.
#
# pwr4exp is licensed under the GNU General Public License (GPL);
# see the LICENSE file for details.
#
# pwr4exp is provided in the hope that it will be useful, but WITHOUT ANY WARRANTY,
# including the implied warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
# For more information, please refer to the GNU General Public License.
###############################################################################

#' Create a Design Object for Power Calculation
#'
#' Generate a design object for power analysis by specifying a model formula and data frame.
#' This object is not a true experimental design as created by design generation procedures,
#' where randomization and unit allocation are performed.
#' Instead, it serves as an object containing all necessary information for power analysis,
#' including design matrices, assumed values of model effects, and other internally
#' calculated parameters.
#'
#' @param formula A right-hand-side [formula] specifying the model for testing treatment effects,
#' with terms on the right of [~] , following [lme4::lmer] syntax for random effects.
#' @param data A data frame with all independent variables specified in the model,
#' matching the design's structure.
#' @param beta One of the optional inputs for fixed effects.
#' A vector of model coefficients where factor variable coefficients correspond
#' to dummy variables created using treatment contrast ([stats::contr.treatment]).
#' @param means One of the optional inputs for fixed effects.
#' A vector of marginal or conditioned means (if factors have interactions).
#' Regression coefficients are required for numerical variables.
#' Either `beta` or `means` must be provided, and their values must strictly follow a specific order.
#' A template can be created to indicate the required input values and their order.
#' See "Details" for more information.
#' @param vcomp A vector of variance-covariance components for random effects, if present.
#' The values must follow a strict order. See "Details".
#' @param sigma2 error variance.
#' @param correlation Specifies residual (R-side) correlation structures using [nlme::corClasses] functions.
#' See "Details" for more information.
#' @param template Default is `FALSE`.
#' If `TRUE` or when only the formula and data are provided, a template for `beta`,
#' `means`, and `vcomp` is generated to indicate the required input order.
#' @param REML Specifies whether to use REML or ML information matrix. Default is `TRUE`.
#' @details
#' - **data**: A long-format data frame is required, as typically used in R for fitting linear models.
#'   This data frame can be created manually or with the help of design creation packages
#'   such as \pkg{agricolae}, \pkg{AlgDesign}, \pkg{crossdes}, or \pkg{FrF2}.
#'   It should include all independent variables specified in the model (e.g., treatments, blocks, subjects).
#'   All the irrelevant variables not specified in the model are ignored.
#'
#' - **template**: Templates are automatically generated when only the formula and data are supplied,
#'   or explicitly if `template = TRUE`. Templates serve as guides for specifying inputs:
#'
#'   - **Template for `beta`**: Represents the sequence of model coefficients.
#'
#'   - **Template for `means`**: Specifies the order of means (for categorical variables)
#'     and/or regression coefficients (for continuous variables), depending on the scenario:
#'
#'     - *Categorical variables without interactions*: Requires marginal means for each level of the categorical variable(s).
#'     - *Interactions among categorical variables*: Requires conditional (cell) means for all level combinations.
#'     - *Numerical variables without interactions*: Requires regression coefficients.
#'       The intercept must also be included if there are no categorical variables in the model.
#'     - *Interactions among numerical variables*: Requires regression coefficients for both main effects and interaction terms.
#'       The intercept must also be included if there are no categorical variables in the model.
#'     - *Categorical-by-numerical interactions*: Requires regression coefficients for the numerical variable
#'       at each level of the categorical variable, as well as marginal means for the levels of the categorical variable.
#'
#'     Note: For models containing only numerical variables, the inputs for `means` and `beta` are identical.
#'     See the "Examples" for illustrative scenarios.
#'
#'   - **Template for `vcomp`**: Represents a variance-covariance matrix, where integers indicate
#'     the order of variance components in the input vector.
#'
#' - **correlation**: Various residual correlation structures can be specified following the instructions from \code{\link[nlme]{corClasses}} constructors.
#'
#'   Note: In original [nlme::corAR1()] and [nlme::corARMA()] with `p=1` and `q=0`, the time variable must be an integer class.
#'   In `pwr4exp`, the time variable can also be a factor class.
#' @return A list object containing all essential components for power calculation.
#' This includes:
#' - Structural components (deStruct): including design matrices for fixed and random effects,
#' variance-covariance matrices for random effects and residuals, etc.
#' - Internally calculated higher-level parameters (deParam), including variance-covariance
#' matrix of beta coefficients (vcov_beta), variance-covariance matrix of variance parameters (vcov_varpar),
#' gradient matrices (Jac_list), etc.
#' @export
#' @examples
#' # Using templates for specifying "means"
#'
#' # Create an example data frame with four categorical variables (factors)
#' # and two numerical variables
#' df1 <- expand.grid(
#'   fA = factor(1:2),
#'   fB = factor(1:2),
#'   fC = factor(1:3),
#'   fD = factor(1:3),
#'   subject = factor(1:10)
#' )
#' df1$x <- rnorm(nrow(df1))  # Numerical variable x
#' df1$z <- rnorm(nrow(df1))  # Numerical variable z
#'
#' ## Categorical variables without interactions
#' # Means of each level of fA and fB are required in sequence.
#' mkdesign(~ fA + fB, df1)$fixeff$means
#'
#' ## Interactions among categorical variables
#' # Cell means for all combinations of levels of fA and fB are required.
#' mkdesign(~ fA * fB, df1)$fixeff$means
#'
#' ## Numerical variables without and with interactions, identical to beta.
#' # Without interactions: Regression coefficients are required
#' mkdesign(~ x + z, df1)$fixeff$means
#'
#' # With interactions: Coefficients for main effects and interaction terms are required.
#' mkdesign(~ x * z, df1)$fixeff$means
#'
#' ## Categorical-by-numerical interactions
#' # Marginal means for each level of fA, and regression coefficients for x
#' # at each level of fA are required.
#' mkdesign(~ fA * x, df1)$fixeff$means
#'
#' ## Three factors with interactions:
#' # Cell means for all 12 combinations of the levels of fA, fB, and fC are required.
#' mkdesign(~ fA * fB * fC, df1)
#'
#' # A design that mixes the above-mentioned scenarios:
#' # - Interactions among three categorical variables (fA, fB, fC)
#' # - A categorical-by-numerical interaction (fD * x)
#' # - Main effects for another numerical variable (z)
#' # The required inputs are:
#' # - Cell means for all combinations of levels of fA, fB, and fC
#' # - Means for each level of fD
#' # - Regression coefficients for x at each level of fD
#' # - Regression coefficients for z
#' mkdesign(~ fA * fB * fC + fD * x + z, df1)$fixeff$means
#'
#' # Using templates for specifying "vcomp"
#'
#' # Assume df1 represents an RCBD with "subject" as a random blocking factor.
#' ## Variance of the random effect "subject" (intercept) is required.
#' mkdesign(~ fA * fB * fC * fD + (1 | subject), df1)$varcov
#'
#' # Demonstration of templates for more complex random effects
#' ## Note: This example is a demo and statistically incorrect for this data
#' ## (no replicates under subject*fA). It only illustrates variance-covariance templates.
#' ## Inputs required:
#' ## - Variance of the random intercept (1st)
#' ## - Covariance between the intercept and "fA2" (2nd)
#' ## - Variance of "fA2" (3rd)
#' mkdesign(~ fA * fB * fC * fD + (1 + fA | subject), df1)$varcov
#'
#' # Power analysis for repeated measures
#'
#' ## Create a data frame for a CRD with repeated measures
#' n_subject <- 6
#' n_trt <- 3
#' n_hour <- 8
#' trt <- c("CON", "TRT1", "TRT2")
#' df2 <- data.frame(
#'   subject = as.factor(rep(seq_len(n_trt * n_subject), each = n_hour)), # Subject as a factor
#'   hour = as.factor(rep(seq_len(n_hour), n_subject * n_trt)),           # Hour as a factor
#'   trt = rep(trt, each = n_subject * n_hour)                           # Treatment as a factor
#' )
#'
#' ## Templates
#' temp <- mkdesign(formula = ~ trt * hour, data = df2)
#' temp$fixeff$means  # Fixed effects means template
#'
#' ## Create a design object
#' # Assume repeated measures within a subject follow an AR1 process with a correlation of 0.6
#' design <- mkdesign(
#'   formula = ~ trt * hour,
#'   data = df2,
#'   means = c(1, 2.50, 3.50,
#'             1, 3.50, 4.54,
#'             1, 3.98, 5.80,
#'             1, 4.03, 5.40,
#'             1, 3.68, 5.49,
#'             1, 3.35, 4.71,
#'             1, 3.02, 4.08,
#'             1, 2.94, 3.78),
#'   sigma2 = 2,
#'   correlation = corAR1(value = 0.6, form = ~ hour | subject)
#' )
#'
#' pwr.anova(design)  # Perform power analysis
#'
#' ## When time is treated as a numeric variable
#' # Means of treatments and regression coefficients for hour at each treatment level are required
#' df2$hour <- as.integer(df2$hour)
#' mkdesign(formula = ~ trt * hour, data = df2)$fixeff$means
#'
#' ## Polynomial terms of time in the model
#' mkdesign(formula = ~ trt + hour + I(hour^2) + trt:hour + trt:I(hour^2), data = df2)$fixeff$means
mkdesign <- function(formula, data,
                     beta = NULL, means = NULL, vcomp = NULL,
                     sigma2 = NULL, correlation = NULL, template = FALSE, REML = TRUE) {
  rownames(data) <- 1:nrow(data) # re-index rownames
  mc <- match.call()
  if (all(names(mc[-1]) %in% c("formula", "data")))
    template <- TRUE
  if (length(formula) > 2) {
    formula <- update(formula, NULL ~ .)
    warning(sprintf("Only rhs formula is required\nformula coerced to: %s.", deparse(formula)),  call. = F)
  }

  deStruct <- mkStruct(formula, data, correlation)
  fixedfr <- deStruct$fxTrms$fixedfr
  L <- means2beta_contrasts(fixedfr)

  if (template) {
    beta <- numeric(ncol(deStruct$fxTrms$X))
    fixeff <- means2beta(L, beta = beta)
    fixeff <- lapply(fixeff, function(x) setNames(seq_along(x), names(x)))
    varcov = deStruct$reTrms$G_temp
    return(list(fixeff = fixeff, varcov = varcov))
  }

  if (!is.null(deStruct$reTrms)) {
    deStruct$reTrms$G@x <- vcomp[deStruct$reTrms$Gind]
  }

  deStruct$rTrms$R <- sigma2 * deStruct$rTrms$R
  fixeff <- means2beta(L, means, beta)
  deStruct$fxTrms$fixeff <- fixeff

  corr <- coef(deStruct$rTrms$corStruct, unconstrained = FALSE)
  varpar <- c(vcomp, corr, sigma2)

  attr(varpar, "vcomp") <- length(vcomp)
  attr(varpar, "corr") <- length(corr)
  attr(varpar, "sigma2") <- length(sigma2)

  vcov_beta <- vcovbeta_vp(varpar, deStruct$fxTrms, deStruct$reTrms, deStruct$rTrms)
  info_mat <- informat(varpar = varpar, deStruct$fxTrms, deStruct$reTrms, deStruct$rTrms, REML)
  vcov_varpar <- solve(info_mat)
  Jac <- numDeriv::jacobian(vcovbeta_vp, x = varpar, fxTrms = deStruct$fxTrms, reTrms = deStruct$reTrms, rTrms = deStruct$rTrms)
  # gradient matrices of vcov(beta) w.r.t. variance parameters
  Jac_list <- lapply(1:ncol(Jac), function(i)
    array(Jac[, i], dim=rep(ncol(deStruct$fxTrms$X), 2)))

  deParam <- list(beta = fixeff$beta, vcov_beta = vcov_beta, vcov_varpar = vcov_varpar, Jac_list = Jac_list)
  object <- list(deStruct = deStruct, deParam = deParam)
  return(object)
}
