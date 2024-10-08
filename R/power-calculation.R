#' Power of omnibus test
#'
#' Calculate power for testing overall effects of treatment factors and their
#' interactions, i.e., statistical power of ANOVA.
#'
#' @param design a design object created using design generating functions.
#' @param alpha significance level (type I error rate), default 0.05
#' @param ... Additional arguments passed to \link[lmerTest]{anova.lmerModLmerTest}
#' for linear mixed models and to \link[car]{Anova} for linear models. The type
#' of sum of squares (SS, default is Type III) and the method for computing denominator
#' degrees of freedom (DDF, default is Satterthwaite's method) can be modified.
#' For balanced designs, types of SS and DDF do not affect results. Note that
#' these additional arguments should be consistent in the design-generating
#' function and \code{pwr.anova} for linear mixed models.
#'
#' @return a data frame with numerator degrees of freedom (NumDF), denominator
#' degrees of freedom (DenDF), non-centrality parameter, type I error rate (alpha),
#' and power.
#'
#' @seealso [designCRD()], [designRCBD()], [designLSD()], [designCOD()], [designSPD()], [designCustom()], and [pwr.contrast()]
#' @export
#' @examples
#' # generate an RCBD
#' rcbd = designRCBD(treatments = c(2, 2), blocks = 10, beta = c(10, 9, 8, 7), VarCov = 10, sigma2 = 9)
#' # power of omnibus test
#' pwr.anova(rcbd, alpha  = 0.05)
pwr.anova <- function(design, alpha = 0.05, ...) {
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
  res <- as.data.frame(do.call(stats::anova, c(list(design$pseu_model), args)))
  res$DenDF <- design$pseu_model@DenDF
  res$`non-centrality` <- res$`F value` * res$NumDF
  res$Fcrt <- stats::qf(p = alpha, df1 = res$NumDF, df2 = res$DenDF, lower.tail = FALSE)
  res$power <- 1 - stats::pf(q = res$Fcrt, df1 = res$NumDF, df2 = res$DenDF, ncp = res$`non-centrality`)
  res <- res[, c("NumDF", "DenDF", "non-centrality", "power")]
  res$alpha <- alpha
  res <- res[, c("NumDF", "DenDF", "non-centrality", "alpha", "power")]
  attr(res, "heading") <- paste("Power Analysis of", design$design)
  class(res) <- c("anova", "data.frame")
  return(res)
}

#' @exportS3Method pwr.anova lmDesign
pwr.anova.lmDesign <- function(design, alpha = 0.05, ...){
  formula <- design$formula
  # Change contrast coding to contr.sum in lm
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
  res <- as.data.frame(do.call(car::Anova, c(list(design$pseu_model), args)))
  res$NumDF <- res$Df
  res$DenDF <- design$pseu_model$df.residual
  res$`non-centrality` <- res$`F value` * res$NumDF
  res$Fcrt <- stats::qf(p = alpha, df1 = res$NumDF, df2 = res$DenDF, lower.tail = FALSE)
  res$power <- 1 - stats::pf(q = res$Fcrt, df1 = res$NumDF, df2 = res$DenDF, ncp = res$`non-centrality`)
  res$alpha <- alpha
  res <- res[, c("NumDF", "DenDF", "non-centrality", "alpha", "power")]
  res <- res[!rownames(res) %in% c("(Intercept)", "Residuals"), ]
  attr(res, "heading") <- paste("Power analysis of", design$design)
  class(res) <- c("anova", "data.frame")
  return(res)
}

#' Power of contrasts
#'
#' Calculate power for testing various contrasts. The same syntax of
#' \link[emmeans:emmeans-package]{emmeans} package is employed to specify contrast types.
#'
#' @param design a design object created using design generating functions.
#' @param specs an argument inherited from \link[emmeans]{emmeans} specifying
#' the names of the factors over which the contrasts are performed.
#' @param method an argument inherited from \link[emmeans]{contrast} specifying
#' the method of contrasts, e.g., pairwise, linear, and polynomials.
#' @param alpha significance level (type I error rate), default 0.05
#' @param ... other arguments passed to \link[emmeans]{contrast}.
#' @export
#' @return a data frame showing the power of the specific contrast
#' @examples
#' rcbd = designRCBD(treatments = c(2, 2), blocks = 10, beta = c(10, 9, 8, 7), VarCov = 10, sigma2 = 9)
#' pwr.contrast(rcbd, specs = ~ facA|facB, method = "pairwise")
pwr.contrast <- function(design, specs, method, alpha = 0.05, ...){
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
  emm.contrast <- as.data.frame(emmeans::contrast(emm, method = method))
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
