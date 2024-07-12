#' Create Latin Square Design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param squares an integer, number of replicated squares
#' @param reuse a character string: "row", "col", or "both", indicating reuse of rows or columns or both when replicate a Latin square
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @export
#' @examples
#' mydesign <- designLSD(treatments = c(2, 2),
#'                       squares = 4,
#'                       reuse = "col",
#'                       beta = c(10, 9, 8, 7),
#'                       VarCov = list(9, 4, 6),
#'                       sigma2 = 10)
designLSD <- function(treatments,
                      squares = 1,
                      reuse = c("row", "col", "both"),
                      formula,
                      beta,
                      VarCov,
                      sigma2,
                      ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.lsd(
    treatments = treatments,
    squares = squares,
    reuse = reuse
  )
  # intended model formula
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ treatment + (1|row) + (1|col)
    } else {
      formula <- y ~ facA*facB + (1|row) + (1|col)
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
  # hint: paramters needed
  # suppressWarnings(
  #   pseudo_model <- lme4::lmer(
  #     formula = formula,
  #     data = cbind(df, y = rnorm(nrow(df), sd = 2)),
  #     control = lme4::lmerControl(check.conv.singular = "ignore")
  #   )
  # )
  # beta <- lme4::fixef(pseudo_model)
  # beta[] <- NA
  #
  # VarCov <- as.data.frame(lme4::VarCorr(pseudo_model, corr = FALSE))
  # VarCov[, 4:5] <- NA
  # VarCov <- VarCov[-nrow(VarCov), -5]
  # names(VarCov)[1] <- c("grouping factor")

  # cat("For setting fixed effects and variance-covariance components for power analysis, refer to:\n")
  # cat("beta:\n")
  # print(beta)
  # cat("\nVarCov:\n")
  # print(VarCov)

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

#' #' Create Completely Randomized Design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param replications number of experimental units in each group
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param sigma2 error variance
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @export
#' @examples
#' mydesign <- designCRD(treatments = c(2, 2),
#'                       replications = 4,
#'                       beta = c(10, 9, 8, 7),
#'                       sigma2 = 10)
designCRD <- function(treatments,
                      replications,
                      formula,
                      beta,
                      sigma2,
                      ...){
  # raw design matrix: data frame without response variable y
  design_df <- df.crd(
    treatments = treatments,
    replications = replications
  )
  # intended model formula
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ treatment
    } else {
      formula <- y ~ facA*facB
    }
  }
  pseu_model <- fit.pseu.model(
    formula = formula,
    data = design_df,
    beta = beta,
    sigma = sqrt(sigma2),
    ...
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


#' Create Randomized Complete Block Design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param blocks number of blocks
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @export
#' @examples
#' mydesign = designRCB(treatments = c(2, 2),
#'                       blocks = 10,
#'                       beta = c(10, 9, 8, 7),
#'                       VarCov = 8,
#'                       sigma2 = 9)
designRCB <- function(treatments,
                      blocks,
                      formula,
                      beta,
                      VarCov,
                      sigma2,
                      ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.rcb(
    treatments = treatments,
    blocks = blocks
  )
  # intended model formula
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ treatment + (1|block)
    } else {
      formula <- y ~ facA*facB + (1|block)
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

#' Create Crossover Design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param squares number of replicated squares
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @export
#' @examples
#' mydesign <- designCrossover(treatments = c(2, 2),
#'                             squares = 4,
#'                             beta = c(10, 9, 8, 7),
#'                             VarCov = list(9, 4),
#'                             sigma2 = 10)
designCrossover <- function(treatments,
                            squares,
                            formula,
                            beta,
                            VarCov,
                            sigma2,
                            ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.crossover(
    treatments = treatments,
    squares = squares
  )
  # intended model formula
  if (missing(formula)) {
    if (length(treatments) == 1) {
      formula <- y ~ treatment + (1|row) + (1|period)
    } else {
      formula <- y ~ facA*facB + (1|row) + (1|period)
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

#' Create Split Plot Design
#'
#' @param trt.main a vector representing the number of levels of main plot factors sequentially.
#' @param trt.sub a vector representing the number of levels of subplot factors sequentially.
#' @param replications number of main plots in each main plot group
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma2 error variance
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
#' @export
#' @examples
#' mydesign <- designSplitplot(trt.main = 2,
#'                             trt.sub = 4,
#'                            replications = 6,
#'                            beta = c(10, -2, 3, 4, 5, 1, 3, 0),
#'                            VarCov = 4,
#'                            sigma2 = 6)
designSplitplot <- function(trt.main,
                            trt.sub,
                            replications,
                            formula,
                            beta,
                            VarCov,
                            sigma2,
                            ...) {
  # raw design matrix: data frame without response variable y
  design_df <- df.splitplot(
    trt.main,
    trt.sub,
    replications
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


#' Create Customized Design
#'
#' @param design.df data frame of design matrix
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model used to assess treatment effects, which should contain fixed components and sources of variation according to experimental designs, assumptions, and objectives. The details of model specification for fixed effects see \code{\link{lm}} and random effects see \code{\link{lmer}}.
#' @param beta effect size: a vector including the expectations of intercept, main effects of factors and interactions between factors. Vector values are in the same order as the coefficients of the model that would be fitted using the exact data with \code{\link{contr.treatment}}.
#' @param VarCov variance-covariance components of random effects. If there are multiple random effects of one grouping factor, supply variance-covariance components in a matrix with the same order as the random effects specified in model formula. If there are multiple grouping factors, supply variance-covariance matrix for each grouping factor in a list.
#' @param sigma2 error variance
#' @param design.name name of the design
#' @param ... additional arguments to be passed to \code{fit.pseu.model}
#'
#' @return a list object containing design name, design data frame, model formula, and a pseudo model object with the expected effect size and variance-covariance components of random effects.
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


