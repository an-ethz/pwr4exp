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

#' Power of omnibus tests
#'
#' Calculates the statistical power for testing the overall effects of treatment
#' factors and their interactions, i.e., power of F-test.
#' @param object a design object created in pwr4exp
#' @param sig.level significance level, default 0.05
#' @param type the type of ANOVA table requested, default Type III
#' @return a data frame with numerator degrees of freedom (`NumDF`), denominator
#' degrees of freedom (`DenDF`), type I error rate (`sig.level`), and `power`.
#' @seealso [mkdesign], [designCRD], [designRCBD], [designLSD], [designCOD], [designSPD], [pwr.summary] and [pwr.contrast]
#' @export
#' @examples
#' # generate an RCBD
#' rcbd <- designRCBD(
#'   treatments = c(2, 2),
#'   label = list(facA = c("1", "2"), facB = c("1", "2")),
#'   blocks = 12,
#'   formula = ~ facA*facB + (1|block),
#'   means = c(32, 35, 30, 37),
#'   vcomp = 4,
#'   sigma2 = 6
#' )
#' # power of omnibus test
#' pwr.anova(rcbd)
#'
#' @importFrom utils  as.roman
pwr.anova <- function(object,
                      sig.level = 0.05,
                      type = c("III", "II", "I", "3", "2", "1")) {
  fixedfr <- object$deStruct$fxTrms$fixedfr
  type <- type[1L]
  if(!is.character(type)) type <- as.character(type)
  type <- match.arg(type)
  if(type %in% c("I", "II", "III"))
    type <- as.character(as.integer(as.roman(type)))
  # Get list of contrast matrices (L) - one for each model term:
  L_list <- if(type == "1") {
    get_contrasts_type1(fixedfr)
  } else if(type == "2") {
    get_contrasts_type2_unfolded(fixedfr)
  } else if(type == "2b") {
    get_contrasts_type2(fixedfr)
  } else if(type == "3") {
    get_contrasts_type3(fixedfr)
  } else if(type == "yates") {
    get_contrasts_yates(fixedfr)
  } else if(type == "marginal") {
    get_contrasts_marginal(fixedfr)
  } else {
    stop("'type' not recognized")
  }
  # Get power for each term and collect in table:
  table <- do.call(rbind, lapply(L_list, function(L) contrastMD(object, L, sig.level)))
  # Format power table and return:
  if(length(nm <- setdiff(names(L_list), rownames(table)))) {
    tab <- array(NA_real_, dim=c(length(nm), 6L),
                 dimnames = list(nm, colnames(table)))
    table <- rbind(table, tab)
  }
  # Format 'type':
  type <- if(type == "marginal") {
    "Marginal"
  } else if (type == "yates" || type == "3b") {
    "Yates"
  } else if(grepl("b|c", type)) {
    alph <- gsub("[0-9]", "", type)
    paste0("Type ", as.roman(as.integer(gsub("b|c", "", type))), alph)
  } else paste("type", as.roman(as.integer(type)))
  attr(table, "heading") <-
    paste("Power of", type, "F-test")
  attr(table, "hypotheses") <- L_list
  class(table) <- c("anova", "data.frame")
  table
}

#' Power of contrasts
#'
#' Computes the statistical power of t-tests for comparisons among means.
#' @param object a design object created in pwr4exp
#' @param which the factor of interest. Multiple factors can be combined using
#' \code{:} or \code{*}, e.g., \code{"facA*facB"}, which represents a single factor that combines
#' the levels of both factors.
#' @param by the variable to condition on
#' @param contrast A character string specifying the contrast method, one of
#' "pairwise", "poly", or "trt.vs.ctrl". Alternatively, a numeric vector defining
#' a single contrast or a (named) list of vectors specifying multiple custom contrasts.
#' If a list is provided, each element must be a vector whose length matches the
#' number of levels of the factor in each `by` group. In multi-factor scenarios,
#' factor levels are combined and treated as a single factor.
#' @param sig.level significance level, default 0.05
#' @param p.adj whether the sig.level should be adjusted using the Bonferroni method, default FALSE
#' @param alternative one- or two-sided test. Can be abbreviated.
#' @param strict use strict interpretation in two-sided case
#' @export
#' @return
#' For each `by` condition, returns a data frame containing the contrast value (`effect`),
#' degrees of freedom (`df`), type I error rate (`sig.level`), `power`, and the test direction
#' (by `alternative`).
#' When multiple `by` conditions are present, the results are returned as a list.
#'
#' @examples
#' rcbd <- designRCBD(
#'   treatments = c(2, 2),
#'   label = list(facA = c("1", "2"), facB = c("1", "2")),
#'   blocks = 12,
#'   formula = ~ facA*facB + (1|block),
#'   means = c(32, 35, 30, 37),
#'   vcomp = 4,
#'   sigma2 = 6
#' )
#'
#' # If contrast is not specified, pairwise comparisons are conducted
#' pwr.contrast(rcbd, which = "facA") # Marginal effect of facA
#' pwr.contrast(rcbd, which = "facA", by = "facB") # Conditional effect of facA within levels of facB
#'
#' # Custom contrast vector, identical to pairwise comparison
#' pwr.contrast(rcbd, which = "facA", contrast = c(1, -1))
#' pwr.contrast(rcbd, which = "facA", by = "facB", contrast = c(1, -1))
#'
#' # A single factor combining two factors
#' pwr.contrast(
#'   rcbd,
#'   which = "facA*facB",
#'   contrast = list(
#'     A1B1vs.A2B1 = c(1, -1, 0, 0),
#'     A1B1vs.A2B2 = c(1, 0, 0, -1)
#'   )
#' )
pwr.contrast <- function(object,
                         which,
                         by = NULL,
                         contrast = c("pairwise", "poly", "trt.vs.ctrl"),
                         sig.level = 0.05,
                         p.adj = FALSE,
                         alternative = c("two.sided", "one.sided"),
                         strict = TRUE) {
  fixedfr <- object$deStruct$fxTrms$fixedfr
  alternative <- match.arg(alternative)
  # Handle contrast argument to allow both character and user-defined contrasts
  if (is.character(contrast)) {
    contrast <- match.arg(contrast)
  }

  ## FIXME: check estimability?
  L_means <- lsmeans_contrasts(fixedfr, which, by)

  if (is.null(by)) {
    if (is.character(contrast)) {
      # Generate contrast matrix for predefined types
      if (contrast == "pairwise") {
        cl <- emmeans::pairwise.emmc(rownames(L_means))
      } else if (contrast == "poly") {
        cl <- emmeans::poly.emmc(rownames(L_means), max.degree = 3)
      } else if (contrast == "trt.vs.ctrl") {
        cl <- emmeans::trt.vs.ctrl.emmc(rownames(L_means))
      }
    } else if (is.list(contrast) | is.vector(contrast)) {
      # Use user-supplied contrast
      if (is.atomic(contrast)) contrast <- list(contrast)
      if (all(lapply(contrast, length) != nrow(L_means)))
        stop("Invalid number of contrast coefficients")
      cl <- do.call(cbind, contrast)
      rownames(cl) <- rownames(L_means)
      colnames <- colnames(cl)
      if (is.null(colnames)) colnames <- character(ncol(cl))
      wic <- which(colnames == "")
      colnames[wic] <- sapply(
        contrast[wic],
        function(x) {paste(c("c(", paste(x, collapse = ","), ")"), collapse = "")}
      )
      colnames(cl) <- colnames
    } else stop("Invalid contrast specification")

    L <- crossprod(as.matrix(cl), L_means)

    if (p.adj) sig.level <- sig.level/nrow(L)

    tab <- apply(L, 1, function(l){
      contrast1D(
        object = object, L = l,
        sig.level = sig.level, alternative = alternative,
        strict = strict
      )
    })
    tab <- do.call(rbind, tab)
    return(tab)
  }

  if (!is.null(by)) {
    tab_list <- lapply(L_means, function(l_means){
      if (is.character(contrast)) {
        # Generate contrast matrix for predefined types
        if (contrast == "pairwise") {
          cl <- emmeans::pairwise.emmc(rownames(l_means))
        } else if (contrast == "poly") {
          cl <- emmeans::poly.emmc(rownames(l_means), max.degree = 3)
        } else if (contrast == "trt.vs.ctrl") {
          cl <- emmeans::trt.vs.ctrl.emmc(rownames(l_means))
        }
      } else if (is.list(contrast) | is.vector(contrast)) {
        # Use user-supplied contrast
        if (is.atomic(contrast)) contrast <- list(contrast)
        if (all(lapply(contrast, length) != nrow(l_means)))
          stop("Invalid number of coefficients in contrast vectors")
        cl <- do.call(cbind, contrast)
        rownames(cl) <- rownames(l_means)
        colnames <- colnames(cl)
        if (is.null(colnames)) colnames <- character(ncol(cl))
        wic <- which(colnames == "")
        colnames[wic] <- sapply(
          contrast[wic],
          function(x) {paste(c("c(", paste(x, collapse = ","), ")"), collapse = "")}
        )
        colnames(cl) <- colnames
      } else stop("Invalid contrast specification")

      L <- crossprod(as.matrix(cl), l_means)

      if (p.adj) sig.level <- sig.level / nrow(L)

      tab <- apply(L, 1, function(l) {
        contrast1D(
          object = object, L = l,
          sig.level = sig.level, alternative = alternative,
          strict = strict
        )
      })

      do.call(rbind, tab)
    })

    by_names <- names(tab_list)
    list_names <- gsub(paste0("^", by), paste0(by, " = "), by_names)
    names(tab_list) <- list_names
    return(tab_list)
  }
}

#' Power for model coefficients
#'
#' Computes the statistical power for testing (t-test) model coefficients.
#'
#' @param object design object
#' @param sig.level significance level, default 0.05
#' @return a data frame containing model coefficients, degrees of freedom (`df`),
#' type I error rate (`sig.level`), power, and the test direction (`alternative`).
#' @export
#' @examples
#' rcbd <- designRCBD(
#'   treatments = c(2, 2),
#'   label = list(facA = c("1", "2"), facB = c("1", "2")),
#'   blocks = 12,
#'   formula = ~ facA*facB + (1|block),
#'   means = c(32, 35, 30, 37),
#'   vcomp = 4,
#'   sigma2 = 6
#' )
#' pwr.summary(rcbd)
pwr.summary <- function(object, sig.level = 0.05) {
  X <- object$deStruct$fxTrms$X
  p <- ncol(X)
  if(p < 1)
    return(as.matrix(contrast1D(object, numeric(0L))))
  Lmat <- diag(p)
  tab <- do.call(rbind, lapply(1:p, function(i) contrast1D(object, Lmat[i, ])))
  rownames(tab) <- colnames(X)
  tab
}

#' Computes power of t-test for one-dimensional contrast matrices
#'
#' @param object design object
#' @param L contrast vector
#' @param method DF approximation method, only "Satterthwaite" available currently
#' @param sig.level significance level, default 0.05
#' @param alternative one- or two-sided test
#' @param strict whether or not use strict interpretation in two-sided case
#' @return A data frame with columns for effect size, degrees of freedom, significance level, power, and test type.
#' @keywords modified lmerTest internal functions
contrast1D <- function(object,
                      L,
                      method=c("Satterthwaite"),
                      sig.level = 0.05,
                      alternative = c("two.sided", "one.sided"),
                      strict = TRUE) {
  param <- object$deParam
  mk_ttable <- function(effect, se, df, sig.level, alternative, strict = TRUE) {
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    ncp <- abs(effect/se)

    if (strict && tside == 2) {
      tvalue <- qt(sig.level/tside, df, lower.tail = FALSE)
      power <- pt(tvalue, df, ncp = ncp, lower.tail = FALSE) +
        pt(-tvalue, df, ncp = ncp, lower.tail = TRUE)
    } else {
      tvalue <- qt(sig.level/tside, df, lower.tail = FALSE)
      power <- pt(tvalue, df, ncp = ncp, lower.tail = FALSE)
    }

    data.frame("effect"=effect, "df"=df,
               "sig.level"=sig.level,
               "power"=power, alternative=alternative,
               check.names=FALSE)
  }
  if(is.matrix(L)) L <- drop(L) # L is a column vector
  stopifnot(is.numeric(L), length(L) == length(param$beta))
  effect <- sum(L * param$beta) # contrast effect (expectation)
  method <- match.arg(method) # currently, just one option
  alternative <- match.arg(alternative)
  if(method == "Satterthwaite") {
    var_con <- sum(L * (param$vcov_beta %*% L)) # variance of contrast
    # Compute denominator DF:
    grad_var_con <-
      vapply(param$Jac_list, function(J) sum(L * (J %*% L)), numeric(1L)) # = {L' Jac L}_i
    satt_denom <- sum(grad_var_con * (param$vcov_varpar %*% grad_var_con)) # g'Ag
    df <- drop(2 * var_con^2 / satt_denom) # denominator DF
  } # TODO: other DF methods
  mk_ttable(effect, sqrt(var_con), df, sig.level, alternative)
}


#' Computes power of F-test for multi-dimensional contrast matrices
#'
#' @param object design object
#' @param L contrast matrix
#' @param sig.level significance level
#' @param eps numeric tolerance
#' @return A data frame
#' @keywords modified lmerTest internal functions
contrastMD <- function(object, L,
                      sig.level = 0.05,
                      eps=sqrt(.Machine$double.eps)) {
  param <- object$deParam
  get_Fstat_ddf <- function(nu, tol=1e-8) {
    fun <- function(nu) {
      if(any(nu <= 2)) 2 else {
        E <- sum(nu / (nu - 2))
        2 * E / (E - (length(nu))) # q = length(nu) : number of t-statistics
      }
    }
    stopifnot(length(nu) >= 1,
              # all(nu > 0), # returns 2 if any(nu < 2)
              all(sapply(nu, is.numeric)))
    if(length(nu) == 1L) return(nu)
    if(all(abs(diff(nu)) < tol)) return(mean(nu))
    if(!is.list(nu)) fun(nu) else vapply(nu, fun, numeric(1L))
  }
  mk_Ftable <- function(ndf, ddf, ncp, sig.level) {
    Fvalue <- stats::qf(sig.level, ndf, ddf, lower.tail = FALSE)
    power <- stats::pf(Fvalue, ndf, ddf, ncp, lower.tail = FALSE)

    data.frame("NumDF"=ndf, "DenDF"=ddf,
               "sig.level"=sig.level,
               "power"=power,
               check.names=FALSE)
  }
  if(!is.matrix(L)) L <- matrix(L, ncol=length(L))
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(param$beta))
  if(nrow(L) == 0L) { # May happen if there are no fixed effects
    x <- numeric(0L)
    return(mk_Ftable(x, x, x, x))
  }
  beta <- param$beta
  vcov_con <- L %*% param$vcov_beta %*% t(L) # (co)-variance of contrasts
  eig_vcov_con <- eigen(vcov_con)
  P <- eig_vcov_con$vectors
  d <- eig_vcov_con$values
  tol <- max(eps * d[1], 0)
  pos <- d > tol
  q <- sum(pos) # rank(vcov_con)
  if(q < nrow(L))
    warning("Contrast is rank deficient and test may be affected")
  # FIXME: matrix, vector, 1D case
  PtL <- crossprod(P, L)[1:q, , drop = FALSE]
  # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector
  grad_PLcov <- lapply(1:q, function(m) {
    vapply(param$Jac_list, function(J) PtL[m, , drop = FALSE] %*% J %*% t(PtL[m, , drop = FALSE]), numeric(1L))
  })
  # Compute degrees of freedom for the q t-statistics:
  nu_m <- vapply(1:q, function(m) {
    2*(d[m])^2 / sum(grad_PLcov[[m]] * (param$vcov_varpar %*% grad_PLcov[[m]])) }, numeric(1L)) # 2D_m^2 / g'Ag
  # Compute ddf for the F-value:
  ddf <- get_Fstat_ddf(nu_m, tol=1e-8)

  # FIXME; two NCP approaches should give the same value
  # get_ncp <- function(L, beta, X, Z, G, R){
  #   if (is.null(Z) | is.null(G)) {
  #     V <- R
  #   } else {
  #     V <- Z %*% G %*% Matrix::t(Z) + R
  #   }
  #   C = solve(t(X) %*% Matrix::solve(V) %*% X)
  #   ncp = drop(t(L %*% beta) %*% solve(L %*% C %*% t(L)) %*% (L %*% beta))
  #   return(ncp)
  # }
  #
  # ncp <- get_ncp(L, beta, object$fxTrms$X, Matrix::object$reTrms$Zt, object$reTrms$G, object$rTrms$R)
  t2 <- drop(PtL %*% beta)^2 / d[1:q]
  ncp <- sum(t2)
  mk_Ftable(ndf = q, ddf = ddf, ncp = ncp, sig.level = sig.level)
}
