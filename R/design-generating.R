#' Creation of Standard Experimental Designs
#'
#' These functions facilitate the creation of standard experimental designs
#' commonly used in agricultural studies for power analysis. Instead of supplying
#' a data frame to [mkdesign], users can specify key design characteristics to generate the design.
#' Other design parameters are consistent with those in [mkdesign].
#'
#' @param treatments An integer vector where each element represents the number of levels
#' of the corresponding treatment factor. A single integer (e.g., \code{treatments = n})
#' specifies one treatment factor with \code{n} levels. When multiple factors are provided,
#' they are arranged in a factorial treatment factor design. For example,
#' \code{treatments = c(2, 3)} creates a 2x3 factorial design with the first factor having 2 levels
#' and the second factor having 3 levels.
#' @param trt.main An integer-valued vector specifying the treatment structure at
#' main plot level for a split plot design, similar to `treatments`.
#' @param trt.sub An integer-valued vector specifying the treatment structure at
#' sub plot level for a split plot design, similar to `treatments`.
#' @param label Optional. A list of character vectors, each corresponding to a treatment factor.
#' The name of each vector specifies the factor's name, and its elements provide the labels for that factor's levels.
#' If no labels are provided, default labels will be used. For a single treatment factor, the default is
#' \code{list(trt = c("1", "2", ...))}, and for two treatment factors, the default is
#' \code{list(facA = c("1", "2", ...), facB = c("1", "2", ...))}.
#' For split-plot designs, the defaults are similar but include the ".main" and ".sub" suffixes for main plot and subplot factors.
#' For example:
#' \code{list(trt.main = c("1", "2", ...), trt.sub = c("1", "2", ...))} and
#' \code{list(facA.main = c("1", "2", ...), facB.main = c("1", "2", ...),
#'       facA.sub = c("1", "2", ...), facB.sub = c("1", "2", ...))}.
#' Label sets should be arranged so that the main plot factors come first, followed by the subplot factors.
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
#' @param formula A right-hand-side [formula] specifying the model for testing treatment effects,
#' with terms on the right of [~] , following [lme4::lmer()] syntax for random effects.
#' If not specified, a default formula with main effects and all interactions is used internally.
#' @param beta One of the optional inputs for fixed effects.
#' A vector of model coefficients where factor variable coefficients correspond
#' to dummy variables created using "contr.treatment".
#' @param means One of the optional inputs for fixed effects.
#' A vector of marginal or conditioned means (if factors have interactions).
#' Regression coefficients are required for numerical variables.
#' Either `beta` or `means` must be provided, and their values must strictly follow a specific order.
#' A template can be created to indicate the required input values and their order.
#' See [mkdesign] for more information.
#' @param vcomp A vector of variance-covariance components for random effects, if present.
#' The values must follow a strict order. See [mkdesign].
#' @param sigma2 error variance.
#' @param template Default is `FALSE`.
#' If `TRUE`, a template for `beta`, `means`, and `vcomp` is generated to indicate the required input order.
#' @param REML Specifies whether to use REML information matrix for variance-covariance parameters instead of ML method.
#' Default is `TRUE`.
#' @details Each function creates a specific design as described below:
#' \describe{
#'   \item{\code{designCRD}}{Completely Randomized Design.
#' By default, the model formula is \code{~ trt} for one factor and
#' \code{~ facA*facB} for two factors, unless explicitly specified. If the
#' `label` argument is provided, the formula is automatically updated with the
#' specified treatment factor names.}
#'   \item{\code{designRCBD}}{Randomized Complete Block Design.
#' The default model formula is \code{~ trt + (1|block)} for one factor and
#' \code{~ facA*facB + (1|block)} for two factors. If `label` is provided, the
#' fixed effect parts of the formula are automatically updated with the specified
#' names. The label of block factor ("block") in the formula is not changeable.}
#'   \item{\code{designLSD}}{Latin Square Design.
#' The default formula is \code{~ trt + (1|row) + (1|col)} for one factor and
#' \code{~ facA*facB + (1|row) + (1|col)} for two factors. If `label` is provided,
#' the fixed effect parts of the formula are automatically updated with the specified
#' names. The labels of row ("row") and column ("col") block factors are not changeable.}
#'   \item{\code{designCOD}}{Crossover Design, which is a special case of LSD
#' with time periods and individuals as blocks. Period blocks are reused when
#' replicating squares.
#' The default formula is \code{~ trt + (1|subject) + (1|period)} for one factor
#' and \code{~ facA*facB + (1|subject) + (1|period)} for two factors. If `label`
#' is provided, the fixed effect parts of the formula are automatically updated
#' with the specified names. Note that "subject" and "period" are the labels for
#' the two blocking factors and cannot be changed.}
#'   \item{\code{designSPD}}{Split Plot Design.
#' The default formula includes the main effects of all treatment factors at
#' both the main and sub-plot levels, their interactions, and the random effects
#' of main plots: \code{~ . + (1|mainplot)}. If `label` is provided, the fixed
#' effect parts of the formula are automatically updated with the specified names.
#' The experimental unit at the main plot level (i.e., the block factor at the
#' subplot level) is always named as "mainplot".}}
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
                      template = FALSE,
                      REML = TRUE) {
  design_df <- df.crd(treatments = treatments, replicates = replicates, label = label)

  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  }
  if (!is.list(label)) label <- list(label)
  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  if (missing(formula))
    formula <- stats::reformulate(paste(names(label), collapse = "*"))

  if (template) {
    return(mkdesign(formula, design_df, template = TRUE))
  }

  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, sigma2 = sigma2, REML = REML)
  return(object)
}

#' @rdname create_designs
#' @export
designRCBD <- function(treatments, label, blocks, formula, beta = NULL, means = NULL, vcomp, sigma2, template = FALSE, REML = TRUE) {
  design_df <- df.rcbd(treatments = treatments, blocks = blocks, label = label)

  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  }
  if (!is.list(label)) label <- list(label)
  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  if (missing(formula))
    formula <- stats::reformulate(c(paste(names(label), collapse = "*"), "(1|block)"))

  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2, REML = REML)
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
                      template=FALSE,
                      REML = TRUE) {
  design_df <- df.lsd(treatments = treatments, squares = squares, reuse = reuse, label = label)

  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  }
  if (!is.list(label)) label <- list(label)
  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  if (missing(formula))
    formula <- stats::reformulate(c(paste(names(label), collapse = "*"), "(1|row)", "(1|col)"))

  # The model formula depends on block indexing. If the block indices are identical across all blocks,
  # then the following formulas apply:
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

  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2, REML = REML)
  return(object)
}

#' @rdname create_designs
#' @export
designCOD <- function(treatments, label, squares = 1, formula, beta =NULL, means =NULL, vcomp, sigma2, template=FALSE, REML = TRUE) {
  design_df <- df.cod(treatments = treatments, squares = squares, label = label)

  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  }
  if (!is.list(label)) label <- list(label)
  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  if (missing(formula))
    formula <- stats::reformulate(c(paste(names(label), collapse = "*"), "(1|subject)", "(1|period)"))

  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2, REML = REML)
  return(object)
}

#' @rdname create_designs
#' @export
designSPD <- function(trt.main, trt.sub, label, replicates, formula, beta = NULL, means = NULL, vcomp, sigma2, template = FALSE, REML = TRUE) {
  design_df <- df.spd(trt.main, trt.sub, replicates, label = label)

  if (missing(label)) {
    label.main <- lapply(trt.main, seq_len)
    label.sub <- lapply(trt.sub, seq_len)
  } else {
    label.main <- label[seq_along(trt.main)]
    label.sub <- label[seq_along(trt.sub) + length(trt.main)]
  }

  if (is.null(names(label.main))) {
    if (length(trt.main) == 1) {
      names(label.main) <- "trt.main"
    } else {
      names(label.main) <- paste0("fac", LETTERS[seq_along(trt.main)], ".main")
    }
  }

  if (is.null(names(label.sub))) {
    if (length(trt.sub) == 1) {
      names(label.sub) <- "trt.sub"
    } else {
      names(label.sub) <- paste0("fac", LETTERS[seq_along(trt.sub)], ".sub")
    }
  }

  label <- c(label.main, label.sub)

  if (missing(formula))
    formula <- stats::reformulate(c(paste(names(label), collapse = "*"), "(1|mainplot)"))

  if (template) {
    return(mkdesign(formula, design_df))
  }
  object <- mkdesign(formula = formula, data = design_df, beta = beta, means = means, vcomp = vcomp, sigma2 = sigma2, REML = REML)
  return(object)
}

