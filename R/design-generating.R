#' Creation of Experimental Designs
#'
#' These functions are used to create design objects for the further evaluation
#' of statistical power.
#'
#' @param treatments An integer-valued vector specifying the treatment structure,
#' in which the length of the vector indicates the number of treatment factors,
#' and each value represents the number of levels for each factor. A maximum of
#' two factors is allowed, and they are arranged in a factorial design.
#' For instance, \code{treatments = n} specifies one treatment factor with n
#' levels, and \code{treatments=c(2,3)} creates a "2x3" factorial design of
#' two treatment factors with 2 and 3 levels, respectively.
#' @param trt.main An integer-valued vector specifying the treatment structure at
#' main plot level for a split plot design, similar to `treatments`.
#' @param trt.sub An integer-valued vector specifying the treatment structure at
#' sub plot level for a split plot design, similar to `treatments`.
#' @param design.df Required input for creating a customized design. A data frame
#' with all independent variables of the design as columns, representing the actual
#' data structure (long format data frame) without response variables.
#' @param design.name Optional input for creating a customized design. A character.
#' @param label Optional. A list of character vectors specifying the names of
#' treatment factors and factor levels. Each vector in the list represents a
#' treatment factor, where the name of the vector specifies the name of the
#' factor, and the values in the vector are the labels for that factor's levels.
#' If not provided, factors and levels for one and two treatment factors are
#' labeled as \code{list(trt = c("1", "2", ...))} and
#' \code{list(facA = c("1", "2", ...), facB = c("1", "2", ...))}, respectively.
#' @param replicates The number of experimental units per treatment in a completely
#' randomized design or the number of experimental units (main plots) per treatment
#' of main plot factors.
#' @param blocks The number of blocks.
#' @param squares The number of replicated squares. By default, 1, i.e., no
#' replicated squares.
#' @param reuse A character string specifying how to replicate squares when
#' there are multiple squares. Options are: "row" for reusing row blocks, "col"
#' for reusing column blocks, or "both" for reusing both row and column blocks
#' to replicate a single square.
#' @param formula A model formula for testing treatment effects in
#' post-experimental data analysis. Use the syntax of \code{\link{lm}} for fixed
#' effects and \code{\link{lmer}} for random effects. The response variable is
#' always denoted as `y`. By default, all interaction terms between treatment
#' factors are included in the formula.
#' @param beta A numeric vector of expected model coefficients, representing the
#' effect sizes. The first element represents the intercept term, corresponding
#' to the mean of the reference level for categorical variables. Subsequent
#' elements correspond to the effect sizes of the independent variables in the
#' order they appear in the model matrix. For categorical variables, each coefficient
#' represents the difference between a non-reference level and the reference level
#' (intercept), as \code{\link{contr.treatment}} contrast coding is used for
#' constructing the model matrix. Ensure that \code{beta} aligns with the columns of the
#' model matrix, including any dummy variables created for categorical predictors.
#' @param VarCov Variance-covariance components of random effects. For multiple
#' random effect groups, supply the variance (for a single random effect term)
#' or variance-covariance matrix (for two or more random effect terms) of each
#' group in a list, following the order in the model formula.
#' @param sigma2 error variance.
#' @param ... Additional arguments passed to the \code{anova} function in
#' \code{lmerTest}. The type of ANOVA table (default is Type III) and the method
#' for computing denominator degrees of freedom (default is Satterthwaite's method)
#' can be modified. For balanced designs, the choice of sum of squares (SS) and
#' degrees of freedom (df) does not affect the results.
#' @details Each function creates a specific design as described below:
#' \describe{
#'   \item{\code{designCRD}}{Completely Randomized Design.
#' By default, the model formula is \code{y ~ trt} for one factor and
#' \code{y ~ facA*facB} for two factors, unless explicitly specified. If the
#' `label` argument is provided, the formula is automatically updated with the
#' specified treatment factor names.}
#'   \item{\code{designRCBD}}{Randomized Complete Block Design.
#' The default model formula is \code{y ~ trt + (1|block)} for one factor and
#' \code{y ~ facA*facB + (1|block)} for two factors. If `label` is provided, the
#' fixed effect parts of the formula are automatically updated with the specified
#' names. The label of block factor ("block") in the formula is not changeable.}
#'   \item{\code{designLSD}}{Latin Square Design.
#' The default formula is \code{y ~ trt + (1|row) + (1|col)} for one factor and
#' \code{y ~ facA*facB + (1|row) + (1|col)} for two factors. If `label` is provided,
#' the fixed effect parts of the formula are automatically updated with the specified
#' names. The labels of row ("row") and column ("col") block factors are not changeable.}
#'   \item{\code{designCOD}}{Crossover Design, which is a special case of LSD
#' with time periods and individuals as blocks. Period blocks are reused when
#' replicating squares.
#' The default formula is \code{y ~ trt + (1|subject) + (1|period)} for one factor
#' and \code{y ~ facA*facB + (1|subject) + (1|period)} for two factors. If `label`
#' is provided, the fixed effect parts of the formula are automatically updated
#' with the specified names. Note that "subject" and "period" are the labels for
#' the two blocking factors and cannot be changed.}
#'   \item{\code{designSPD}}{Split Plot Design.
#' The default formula includes the main effects of all treatment factors at
#' both the main and sub-plot levels, their interactions, and the random effects
#' of main plots: \code{y ~ . + (1|mainplot)}. If `label` is provided, the fixed
#' effect parts of the formula are automatically updated with the specified names.
#' The experimental unit at the main plot level (i.e., the block factor at the
#' subplot level) is always denoted as "mainplot".}
#' \item{\code{designCustom}}{Customized Design.}
#' }
#'
#' @rdname create_designs
#' @aliases design.CRD design.RCBD design.LSD design.COD design.SPD design.Custom
#' @return a list with the design name, data structure (data frame), model
#' formula, and a pseudo model object with the expected fixed and random effects.
#' @seealso [pwr.anova()], [pwr.contrast()]
#' @examples
#' # Example 1: Evaluate the power of a CRD with one treatment factor
#'
#' ## Create a design object
#'
#' crd <- designCRD(
#'   treatments = 4, # 4 levels of one treatment factor
#'   replicates = 12, # 12 units per level, 48 units totally
#'   # mean of level1, and the means of other levels minus level1, respectively
#'   beta = c(30, -2, 3, 5),
#'   sigma2 = 10 # error variance
#' )
#'
#' ## power of omnibus test
#' pwr.anova(crd)
#'
#' ## power of contrast
#' pwr.contrast(crd, specs = "trt", method = "pairwise") # pairwise comparisons
#' pwr.contrast(crd, specs = "trt", method = "poly") # polynomial contrasts
#'
#' # Example 2: Evaluate the power of an RCBD with 2 x 2 factorial treatments
#'
#' # Treatment factors are A (A1 vs. A2) and B (B1 vs. B2).
#' # To illustrate how to provide `beta`, treatment means are presented:
#' #     B1  B2
#' # A1  20  24
#' # A2  17  22
#' #
#' # From these means, we calculate:
#' # 1. the mean of reference level (A1B1): 20
#' # 2. the effect of A2 alone: Effect_A2 = A2B1 - A1B1 = 17 - 20 = -3
#' # 3. the effect of B2 alone: Effect_A2 = A1B2 - A1B1 = 24 - 20 = 4
#' # 4. the interaction effect of A2 and B2:
#' #    Interaction_A2B2 = A2B2 - A2B1 - A1B2 + A1B1 = 22 - 17 - 24 + 20 = 1, representing
#' #    the additional effect of combining A2B2 compared to what would be expected
#' #    from the sum of individual effects of A2 and B2.
#'
#' # The `beta` vector is constructed as:
#' # beta = c(mean_A1B1, Effect_A2, Effect_B2, Interaction_A2B2)
#' # beta = c(20, -3, 4, 1)
#'
#' ## Create a design object
#'
#' rcbd <- designRCBD(
#'   # 2x2 factorial design
#'   treatments = c(2, 2),
#'   # Specify treatment names
#'   label = list(A = c("A1", "A2"), B = c("B1", "B2")),
#'   # 12 blocks, totaling 48 experimental units
#'   blocks = 12,
#'   # Mean of the reference level and effect sizes as calculated above
#'   beta = c(20, -3, 4, 1),
#'   # Variance of block effects (between-block variance)
#'   VarCov = 30,
#'   # Error variance (within-block variance)
#'   sigma2 = 20
#' )
#'
#' ## power of omnibus test
#'
#' pwr.anova(rcbd)
#'
#' ## power of B2 vs. B1 at each level of A
#' pwr.contrast(rcbd, specs = ~B|A, method = "pairwise")
#'
#' # More examples are available in the package vignette("pwr4exp")
#' # and on the package website: https://an-ethz.github.io/pwr4exp/
#' @export
designCRD <- function(treatments, label, replicates, formula, beta, sigma2) {
  design_df <- df.crd(treatments = treatments, replicates = replicates, label = label)
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ trt
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", names(label)))
      }
    } else {
      formula <- y ~ facA*facB
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*")))
      }
    }
  }
  pseu_model <- fit.pseu.model(formula = formula, data = design_df, beta = beta, sigma = sqrt(sigma2))
  output <- list(design = "Completely Randomized Design", design_df = design_df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }
  return(output)
}

#' @rdname create_designs
#' @export
designRCBD <- function(treatments, label, blocks, formula, beta, VarCov, sigma2, ...) {
  design_df <- df.rcbd(treatments = treatments, blocks = blocks, label = label)
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ trt + (1|block)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", names(label), "+ (1|block)"))
      }
    } else {
      formula <- y ~ facA*facB + (1|block)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|block)"))
      }
    }
  }
  pseu_model <- fit.pseu.model(formula = formula, data = design_df, beta = beta, VarCov = VarCov, sigma = sqrt(sigma2), ...)
  output <- list(design = "Randomized Complete Block Design", design_df = design_df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}

#' @rdname create_designs
#' @export
designLSD <- function(treatments, label, squares = 1, reuse = c("row", "col", "both"), formula, beta, VarCov, sigma2, ...) {
  design_df <- df.lsd(treatments = treatments, squares = squares, reuse = reuse, label = label)
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ trt + (1|row) + (1|col)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", names(label), "+ (1|row) + (1|col)"))
      }
    } else {
      formula <- y ~ facA*facB + (1|row) + (1|col)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|row) + (1|col)"))
      }
    }
    # The model specification depends on the level coding of blocks.
    # if (squares > 1) {
    #   if (reuse == "row") {
    #     formula <- stats::update(formula, . ~ . - (1|col) + (1|square/col))
    #   }
    #   if (reuse == "col") {
    #     formula <- stats::update(formula, . ~ . - (1|row) + (1|square/row))
    #   }
    #   if (reuse == "both") {
    #     formula <- stats::update(formula, . ~ . - (1|row) - (1|col) + (1|square) + (1|square:row) + (1|square:col))
    #   }
    # }
  }
  pseu_model <- fit.pseu.model(formula = formula, data = design_df, beta = beta, VarCov = VarCov, sigma = sqrt(sigma2), ...)
  output <- list(design = "Latin Square Design", design_df = design_df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }
  return(output)
}

#' @rdname create_designs
#' @export
designCOD <- function(treatments, label, squares = 1, formula, beta, VarCov, sigma2, ...) {
  design_df <- df.cod(treatments = treatments, squares = squares, label = label)
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ trt + (1|subject) + (1|period)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", names(label), "+ (1|subject) + (1|period)"))
      }
    } else {
      formula <- y ~ facA*facB + (1|subject) + (1|period)
      if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|subject) + (1|period)"))
      }
    }
  }
  pseu_model <- fit.pseu.model(formula = formula, data = design_df, beta = beta, VarCov = VarCov, sigma = sqrt(sigma2), ...)
  output <- list(design = "Crossover Design", design_df = design_df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }
  return(output)
}

#' @rdname create_designs
#' @export
designSPD <- function(trt.main, trt.sub, label, replicates, formula, beta, VarCov, sigma2, ...) {
  design_df <- df.spd(trt.main, trt.sub, replicates, label = label)
  if (missing(formula)) {
    formula <- y ~ (1|mainplot)
    if (length(trt.main) == 1) {
      if (length(trt.sub) == 1) {
        formula <- stats::update(formula, . ~ . + trt.main*trt.sub)
        if (!missing(label)) {
          formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|mainplot)"))
        }
      } else {
        formula <- stats::update(formula, . ~ . + trt.main*facA.sub*facB.sub)
        if (!missing(label)) {
          formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|mainplot)"))
          }
        }
    }
    if (length(trt.main) == 2) {
      if (length(trt.sub) == 1) {
        formula <- stats::update(formula, . ~ . + facA.main*facB.main*trt.sub)
        if (!missing(label)) {
          formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|mainplot)"))
          }
      } else {
        formula <- stats::update(formula, . ~ . + facA.main*facB.main*facA.sub*facB.sub)
        if (!missing(label)) {
        formula <- stats::as.formula(paste("y ~", paste(names(label), collapse = "*"), "+ (1|mainplot)"))
        }}
    }
  }
  pseu_model <- fit.pseu.model(formula = formula, data = design_df, beta = beta, VarCov = VarCov, sigma = sqrt(sigma2), ...)
  output <- list(design = "Split Plot Design", design_df = design_df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}

#' @rdname create_designs
#' @export
designCustom <- function(design.df, formula, beta, VarCov, sigma2, design.name, ...) {
  if (missing(design.name)) {design.name = "Customized Design"}
  pseu_model <- fit.pseu.model(formula = formula, data = design.df, beta = beta, VarCov = VarCov, sigma = sqrt(sigma2), ...)
  output <- list(design = design.name, design_df = design.df, formula = formula, pseu_model = pseu_model)
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }
  return(output)
}
