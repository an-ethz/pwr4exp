#' Creation of Experimental Designs
#'
#' These functions are used to create some common designs in agricultural studies by specifying the design characteristics.
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
#' @param label A optional list of strings specifying the names of
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
#' @param formula a right-hand-side formula to be used for testing treatment effects
#' Use the syntax of \code{\link{lm}} for fixed effects and \link[lme4]{lmer} for random effects.
#' If not specified, a default formula with main effects and all interactions is used internally.
#' @param beta model coefficients, an optional input of fixed effects.
#' The coefficients of factor variables are the coefficients of dummy variables created
#' using "contr.treatment".
#' @param means marginal means or conditioned means if factors have interaction.
#' Regression coefficients should be provided for numerical variables. Use `template = TRUE` to create
#' a template for these inputs.
#' @param vcomp variance-covariance parameters of random effects provided in a vector.
#' Use `template = TRUE` to create a template for these inputs.
#' @param sigma2 error variance.
#' @param template if TRUE, a template for `beta`, `means`, and `vcomp` is
#' generated to indicate the orders of inputs
#'
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
#' @aliases design.CRD design.RCBD design.LSD design.COD design.SPD
#' @return a list with critical components for power analysis
#' @seealso [mkdesign], [pwr.anova()], [pwr.contrast()]
#' @examples
#' # Example 1: Evaluate the power of a CRD with one treatment factor
#'
#' ## Create a design object
#'
#' crd <- designCRD(
#'   treatments = 4, # 4 levels of one treatment factor
#'   replicates = 12, # 12 units per level, 48 units totally
#'   means = c(30, 28, 33, 35), # means of the 4 levels
#'   sigma2 = 10 # error variance
#' )
#'
#' ## power of omnibus test
#' pwr.anova(crd)
#'
#' ## power of contrast
#' pwr.contrast(crd, which = "trt", contrast = "pairwise") # pairwise comparisons
#' pwr.contrast(crd, which = "trt", contrast = "poly") # polynomial contrasts
#'
#' # More examples are available in the package vignette("pwr4exp")
#' # and on the package website: https://an-ethz.github.io/pwr4exp/
#'
#' @export
designCRD <- function(treatments,
                      label,
                      replicates,
                      formula,
                      beta = NULL,
                      means = NULL,
                      sigma2,
                      template = FALSE) {
  design_df <- df.crd(treatments = treatments, replicates = replicates, label = label)
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- ~ trt
      if (!missing(label)) {
        formula <- stats::reformulate(names(label))
      }
    } else {
      formula <- ~ facA*facB
      if (!missing(label)) {
        formula <- stats::reformulate(paste(names(label), collapse = "*"))
      }
    }
  }
  if (template) {
    return(mkdesign(formula, design_df, template = TRUE))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, sigma2 = sigma2)
  return(object)
}

#' @rdname create_designs
#' @export
designRCBD <- function(treatments, label, blocks, formula, beta = NULL, means = NULL, vcomp, sigma2, template = FALSE) {
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
  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2)
  return(object)
}


#' @rdname create_designs
#' @export
designLSD <- function(treatments,
                      label,
                      squares = 1,
                      reuse = c("row", "col", "both"),
                      formula,
                      beta = NULL,
                      means= NULL,
                      vcomp,
                      sigma2,
                      template=FALSE) {
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
  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2)
  return(object)
}

#' @rdname create_designs
#' @export
designCOD <- function(treatments, label, squares = 1, formula, beta =NULL, means =NULL, vcomp, sigma2, template=FALSE) {
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
  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2)
  return(object)
}

#' @rdname create_designs
#' @export
designSPD <- function(trt.main, trt.sub, label, replicates, formula, beta = NULL, means = NULL, vcomp, sigma2, template = FALSE) {
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
  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2)
  return(object)
}

