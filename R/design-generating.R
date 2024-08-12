#' Create a Completely Randomized Design
#'
#' @param treatments a vector specifying the treatment structure,
#' in which the length of the vector indicates the number of treatment factors,
#' and each value represents the number of levels for each factor.
#' A maximum of two factors is allowed, and they are arranged in a factorial design.
#' For example, \code{treatments=c(2,3)} specifies a "2x3" factorial design of
#' two treatment factors with 2 and 3 levels, respectively.
#' @param label Optional. A list specifying the names of treatment factors and
#' their corresponding level labels. Each element of the list represents a factor,
#' with the name of the element being the name of the factor. The value of each
#' element should be a vector containing the labels for that factor's levels.
#' By default, the name of one treatment factor is "trt," and the names of two
#' treatment factors are "facA" and "facB," respectively.
#' The levels of each factor are coded as "1", "2", ..., if not specified in `label`.
#' For example, `list(trt = c("A1", "A2"))` assigns the name "trt" to the treatment factor,
#' with its levels labeled as "A1" and "A2." For multiple factors, such as
#' `list(A = c("A1", "A2"), B = c("B1", "B2", "B3"))`, the first factor is named
#' "A" with levels "A1" and "A2," and the second factor is named "B" with levels
#' "B1," "B2," and "B3."
#' @param replicates The number of experimental units per group.
#' In a one-factor design, 'group' refers to the levels of that factor.
#' In a two-factor design, 'group' refers to each combination of levels between the two factors.
#' @param formula a symbolic description of the model that would be used to test
#' the treatment effects in post-experimental data analysis.
#' The way of writing formula , see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' The default formula for one treatment factor and two treatment factors are
#' \code{y ~ trt } and \code{y ~ facA*facB}, respectively,
#' and when the `label` argument is supplied, the names of treatment factors are
#' replaced by the corresponding names specified in `label`.
#' The response variable in the formula is always denoted as `y`.
#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param sigma2 error variance
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designCRD <- function(treatments,
                      label,
                      replicates,
                      formula,
                      beta,
                      sigma2){
  # raw design matrix: data frame without response variable y
  design_df <- df.crd(
    treatments = treatments,
    replicates = replicates,
    label = label
  )
  # default model formula
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
  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    sigma = sqrt(sigma2)
  )
  output <- list(design = "Completely randomized design",
                 design_df = design_df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}


#' Create an Randomized Complete Block Design
#'
#' @param treatments a vector specifying the treatment structure, in which the
#' length of the vector indicates the number of treatment factors, and each value
#' represents the number of levels for each factor. For example, \code{treatments=c(2,3)}
#' specifies two treatment factors with 2 and 3 levels, respectively.
#' @param label Optional. A list specifying the names of treatment factors and
#' their corresponding level labels. Each element of the list represents a factor,
#' with the name of the element being the name of the factor. The value of each
#' element should be a vector containing the labels for that factor's levels.
#' By default, the name of one treatment factor is "trt," and the names of two
#' treatment factors are "facA" and "facB," respectively.
#' The levels of each factor are coded as "1", "2", ..., if not specified in `label`.
#' For example, `list(trt = c("A1", "A2"))` assigns the name "trt" to the treatment factor,
#' with its levels labeled as "A1" and "A2." For multiple factors, such as
#' `list(A = c("A1", "A2"), B = c("B1", "B2", "B3"))`, the first factor is named
#' "A" with levels "A1" and "A2," and the second factor is named "B" with levels
#' "B1," "B2," and "B3."
#' @param blocks the number of blocks
#' @param formula a symbolic description of the model that would be used to test
#' the treatment effects in post-experimental data analysis.
#' The way of writing formula , see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' The default formula for one treatment factor and two treatment factors are
#' \code{y ~ trt + (1|block)} and \code{y ~ facA*facB +(1|block)}, respectively,
#' and when the `label` argument is supplied, the names of treatment factors are
#' replaced by the corresponding names specified in `label`.
#' The response variable in the formula is always denoted as `y`, and block factor
#' is always denoted as `block`.
#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor,
#' provide the variance-covariance matrix.
#' If there are multiple grouping factors,
#' supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default,
#' and the method for computing the denominator degrees of freedom,
#' with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designRCBD <- function(treatments,
                       label,
                       blocks,
                       formula,
                       beta,
                       VarCov,
                       sigma2,
                       ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.rcbd(
    treatments = treatments,
    blocks = blocks, label = label
  )
  # default formula
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

  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    VarCov = VarCov,
    sigma = sqrt(sigma2),
    ...
  )
  output <- list(design = "Randomized complete block design",
                 design_df = design_df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}

#' Create a Latin Square Design
#'
#' @param treatments a vector specifying the treatment structure, in which the
#' length of the vector indicates the number of treatment factors, and each value
#' represents the number of levels for each factor. For example, \code{treatments=c(2,3)}
#' specifies two treatment factors with 2 and 3 levels, respectively.
#' @param label Optional. A list specifying the names of treatment factors and
#' their corresponding level labels. Each element of the list represents a factor,
#' with the name of the element being the name of the factor. The value of each
#' element should be a vector containing the labels for that factor's levels.
#' By default, the name of one treatment factor is "trt," and the names of two
#' treatment factors are "facA" and "facB," respectively.
#' The levels of each factor are coded as "1", "2", ..., if not specified in `label`.
#' For example, `list(trt = c("A1", "A2"))` assigns the name "trt" to the treatment factor,
#' with its levels labeled as "A1" and "A2." For multiple factors, such as
#' `list(A = c("A1", "A2"), B = c("B1", "B2", "B3"))`, the first factor is named
#' "A" with levels "A1" and "A2," and the second factor is named "B" with levels
#' "B1," "B2," and "B3."
#' @param squares the number of replicated squares
#' @param reuse a character string: "row", "col", or "both", indicating reuse of
#' rows or columns or both when replicate a Latin square
#' @param formula a symbolic description of the model that would be used to test
#' the treatment effects in post-experimental data analysis.
#' The way of writing formula , see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' The default formula for one treatment factor and two treatment factors are
#' \code{y ~ trt + (1|row) + (1|col)} and \code{y ~ facA*facB + (1|row) + (1|col)}, respectively,
#' and when the `label` argument is supplied, the names of treatment factors are
#' replaced by the corresponding names specified in `label`.
#' Column factor `(1|col)` in the default formula is replaced by `(1|square/col)` when
#' row factors are reused, and row factor `(1|row)` is replaced by `(1|square/row)` when
#' column factors are reused. When both row and column factors are resued, the random
#' effects in the model are denoted as `(1|square) + (1|square:row) + (1|square:col)`.
#' The response variable in the formula is always denoted as `y`, and row and column block factors
#' are always denoted as `row` and `col`.
#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default,
#' and the method for computing the denominator degrees of freedom,
#' with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designLSD <- function(treatments,
                      label,
                      squares = 1,
                      reuse = c("row", "col", "both"),
                      formula,
                      beta,
                      VarCov,
                      sigma2,
                      ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.lsd(treatments = treatments, squares = squares, reuse = reuse, label = label)
  # intended model formula
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
    if (squares > 1) {
      if (reuse == "row") {
        formula <- stats::update(formula, . ~ . - (1|col) + (1|square/col))
      }
      if (reuse == "col") {
        formula <- stats::update(formula, . ~ . - (1|row) + (1|square/row))
      }
      if (reuse == "both") {
        formula <- stats::update(formula, . ~ . - (1|row) - (1|col) + (1|square) + (1|square:row) + (1|square:col))
      }
    }
  }

  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    VarCov = VarCov,
    sigma = sqrt(sigma2),
    ...
  )

  output <- list(design = "Latin square design",
                 design_df = design_df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}

#' Create a Crossover Design
#'
#' A crossover design is a specific type case of LSD, where each time period serves as a block,
#' and period blocks are reused during replicating squares.
#'
#' @param treatments a vector specifying the treatment structure, in which the
#' length of the vector indicates the number of treatment factors, and each value
#' represents the number of levels for each factor. For example, \code{treatments=c(2,3)}
#' specifies two treatment factors with 2 and 3 levels, respectively.
#' @param label Optional. A list specifying the names of treatment factors and
#' their corresponding level labels. Each element of the list represents a factor,
#' with the name of the element being the name of the factor. The value of each
#' element should be a vector containing the labels for that factor's levels.
#' By default, the name of one treatment factor is "trt," and the names of two
#' treatment factors are "facA" and "facB," respectively.
#' The levels of each factor are coded as "1", "2", ..., if not specified in `label`.
#' For example, `list(trt = c("A1", "A2"))` assigns the name "trt" to the treatment factor,
#' with its levels labeled as "A1" and "A2." For multiple factors, such as
#' `list(A = c("A1", "A2"), B = c("B1", "B2", "B3"))`, the first factor is named
#' "A" with levels "A1" and "A2," and the second factor is named "B" with levels
#' "B1," "B2," and "B3."
#' @param squares the number of replicated squares
#' @param formula a symbolic description of the model that would be used to test
#' the treatment effects in post-experimental data analysis.
#' The way of writing formula , see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' The default formula for one treatment factor and two treatment factors are
#' \code{y ~ trt + (1|subject) + (1|period)} and \code{y ~ facA*facB + (1|subject) + (1|period)}, respectively,
#' and when the `label` argument is supplied, the names of treatment factors are
#' replaced by the corresponding names specified in `label`.
#' The response variable in the formula is always denoted as `y`, and row and column block factors
#' are always denoted as `subject` and `period`.
#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default, and the method for computing the denominator degrees of freedom, with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designCOD <- function(treatments,
                      label,
                      squares,
                      formula,
                      beta,
                      VarCov,
                      sigma2,
                      ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.cod(
    treatments = treatments,
    squares = squares, label = label
  )
  # intended model formula
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

  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    VarCov = VarCov,
    sigma = sqrt(sigma2),
    ...
  )
  output <- list(design = "Crossover design",
                 design_df = design_df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}

#' Create a Split Plot Design
#'
#' @param trt.main a vector specifying the treatment structure at main plot level, in which the
#' length of the vector indicates the number of main plot factors, and each value
#' represents the number of levels for each factor.
#' @param trt.sub a vector specifying the treatment structure at main plot level, in which the
#' length of the vector indicates the number of main plot factors, and each value
#' represents the number of levels for each factor.
#' @param label Optional. A list specifying the names of treatment factors and
#' their corresponding level labels. Each element of the list represents a factor,
#' with the name of the element being the name of the factor. The value of each
#' element should be a vector containing the labels for that factor's levels.
#' By default, the name of one treatment factor is "trt," and the names of two
#' treatment factors are "facA" and "facB," respectively.
#' The levels of each factor are coded as "1", "2", ..., if not specified in `label`.
#' For example, `list(trt = c("A1", "A2"))` assigns the name "trt" to the treatment factor,
#' with its levels labeled as "A1" and "A2." For multiple factors, such as
#' `list(A = c("A1", "A2"), B = c("B1", "B2", "B3"))`, the first factor is named
#' "A" with levels "A1" and "A2," and the second factor is named "B" with levels
#' "B1," "B2," and "B3."
#' @param replicates the number of main plots per treatment group of main plot factors
#' @param formula a symbolic description of the model that would be used to test
#' the treatment effects in post-experimental data analysis.
#' The way of writing formula , see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' The response variable in the formula is always denoted as `y`, and the experimental
#' unit at main plot level (aka. block factor at subplot level) is always denoted as `mainplots`.
#' The default formula contains main effects of all treatment factors at main and sub- plot levels,
#' interactions between them factors, and random effects of main plots `y ~ . + (1|mainplots)`.
#' replaced by the corresponding names specified in `label`.
#' Column factor `(1|col)` in the default formula is replaced by `(1|square/col)` when
#' row factors are reused, and row factor `(1|row)` is replaced by `(1|square/row)` when
#' column factors are reused. When both row and column factors are resued, the random
#' effects in the model are denoted as `(1|square) + (1|square:row) + (1|square:col)`.

#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default,
#' and the method for computing the denominator degrees of freedom,
#' with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designSPD <- function(trt.main,
                      trt.sub,
                      label,
                      replicates,
                      formula,
                      beta,
                      VarCov,
                      sigma2,
                      ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.spd(
    trt.main,
    trt.sub,
    replicates, label = label
  )
  # intended model formula
  if (missing(formula)) {
    formula <- y ~ (1|mainplots)
    if (length(trt.main) == 1) {
      if (length(trt.sub) == 1) {
        formula <- stats::update(formula, . ~ . + trt.main*trt.sub)
      } else {formula <- stats::update(formula, . ~ . + trt.main*facA.sub*facB.sub)}
    }
    if (length(trt.main) == 2) {
      if (length(trt.sub) == 1) {
        formula <- stats::update(formula, . ~ . + facA.main*facB.main*trt.sub)
      } else {formula <- stats::update(formula, . ~ . + facA.main*facB.main*facA.sub*facB.sub)}
    }
  }
  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    VarCov = VarCov,
    sigma = sqrt(sigma2),
    ...
  )
  output <- list(design = "Split plot design",
                 design_df = design_df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}


#' Create a Customized Design
#'
#' @param design.df a data frame containing values of all independent variables as columns
#' @param formula an object of class \code{\link{formula}}: a symbolic description
#' of the model that would be used to test effects in post-experimental data analysis.
#' For details on model specification, see \code{\link{lm}} for fixed effects and
#' \code{\link{lmer}} for random effects.
#' @param beta a vector of the expectations of model coefficients.The coefficients
#' for categorical variables are the coefficients of dummy variables created using
#' \code{\link{contr.treatment}} contrast coding.
#' @param VarCov variance-covariance components of random effects.
#' If there are multiple random effects of one grouping factor, provide the variance-covariance matrix.
#' If there are multiple grouping factors, supply the variance-covariance matrix of each grouping factor in a list.
#' @param sigma2 error variance
#' @param design.name optional, name of the design
#' @param ... other arguments passed to the \code{anova} function in \code{lmerTest}.
#' The type of ANOVA table, with Type III as the default, and the method for computing the denominator degrees of freedom, with Satterthwaite's method as the default, can be changed.
#' For more details, see \link[lmerTest]{anova.lmerModLmerTest}.
#'
#' @return an object containing design name, design data frame, model formula,
#' and a pseudo model object with the expected effect size and variance-covariance
#' components of random effects.
#' @export
designCustom <- function(design.df,
                         formula,
                         beta,
                         VarCov,
                         sigma2,
                         design.name,
                         ...) {

  if (missing(design.name)) {
    design.name = "Customized design"

  }
  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design.df,
    beta = beta,
    VarCov = VarCov,
    sigma = sqrt(sigma2),
    ...
  )
  output <- list(design = design.name,
                 design_df = design.df,
                 formula = formula,
                 pseu_model = pseu_model
  )

  # Set the class of the output list based on the class of pseu_model
  if (inherits(pseu_model, "customLmerMod")) {
    class(output) <- "lmmDesign"
  } else if (inherits(pseu_model, "lm")) {
    class(output) <- "lmDesign"
  }

  return(output)
}


