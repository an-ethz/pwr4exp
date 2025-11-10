###############################################################################
# This file is part of the pwr4exp package.
#
# pwr4exp is licensed under the GNU General Public License (GPL); see the LICENSE file for details.
#
# pwr4exp is provided in the hope that it will be useful, but WITHOUT ANY WARRANTY,
# including the implied warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
# For more information, please refer to the GNU General Public License.
#
# This file contains functions that have been copied or modified from internal functions
# of the lme4 and lmerTest packages. The original implementations are licensed under GPL (>= 2),
# and their use here complies with the terms of the GNU GPL.
#
# Attribution is provided in the roxygen2 documentation via:
#   @keywords lme4 internal functions
#   @keywords lmerTest internal functions
###############################################################################


################################################################################
##########                Create structural matrices                 ###########
################################################################################

#' Building Design Matrices and Covariance Structures for Linear Mixed Models
#'
#' Constructs design matrices for fixed and random effects, along with
#' variance-covariance structures for random effects (G-side) and residuals (R-side).
#'
#' @param formula a model formula.
#' @param data a data frame containing the variables used in the model.
#' @param correlation a `corStruct` object created by [corClasses] functions.
#' If NULL, an identity matrix is assumed.
#' @return A list containing:
#' - `data`: Processed data frame with NA values omitted.
#' - `fxTrms`: Fixed-effects design structure, including model frame and design matrix.
#' - `reTrms`: Random-effects structure (if applicable), including grouping factors,
#' design matrices, and variance-covariance matrix.
#' - `rTrms`: Residual structure (R-side variance-covariance components).
#' - `formula`: Expanded model formula.
#' @importFrom stats reformulate model.matrix terms
#' @keywords internal
mkStruct <- function(formula, data, correlation) {
  corform <- attr(correlation, "formula")
  allvars <- unique(c(all.vars(formula), all.vars(corform)))
  data <- stats::model.frame(stats::reformulate(allvars), data)
  if (nrow(data) == 0L) stop("0 (non-NA) cases")
  attr(data, "terms") <- NULL
  data <- factorize(formula, data, char.only = TRUE)
  n <- nrow(data)

  formula <- expandDoubleVerts(formula)

  # fixed effects
  fixedform <- nobars(formula)
  fixedfr <- stats::model.frame(fixedform, data,
                         na.action = na.omit,
                         drop.unused.levels = TRUE)
  factor_vars <- names(fixedfr)[sapply(fixedfr, is.factor)]
  contrasts_list <- lapply(factor_vars, function(var) "contr.treatment")
  names(contrasts_list) <- factor_vars
  X <- stats::model.matrix(fixedform, fixedfr, contrasts_list)
  chkRank.drop.cols(X)
  checkScaleX(X)
  fxTrms <- list(fixedfr = fixedfr, X = X)

  # random effects
  if (!is.null(bars <- findbars(formula))) {
    reTrms <- mkReTrms(bars, data)
    checkNlevels(reTrms$flist, n = n)
    checkZdims(reTrms$Ztlist, n = n)
    if (anyNA(reTrms$Zt)) {
      stop("NA in Z (random-effects model matrix): ", "please use ",
           shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
    }
    checkZrank(reTrms$Zt, n = n)
  } else reTrms <- NULL

  # residuals
  rTrms <- mkRTrms(data, correlation)

  list(data = data, fxTrms = fxTrms, reTrms = reTrms, rTrms = rTrms, formula = formula)
}

#' Design Matrices and Variance Components for Random Effects
#'
#' Adapted from lme4, this function constructs the design matrix (Z),
#' variance-covariance matrix (G), etc.
#'
#' @param bars a list of parsed random-effects terms
#' @param fr a model frame in which to evaluate these terms
#' @param drop.unused.levels (logical) drop unused factor levels?
#' @param reorder.terms arrange random effects terms in decreasing order of number of groups (factor levels)?
#' @param reorder.vars arrange columns of individual random effects terms in alphabetical order?
#' @return A list with the following components:
#' - Zt: Transposed random-effects design matrix.
#' - G: Variance-covariance matrix for random effects.
#' - Gind: Index mapping variance-covariance parameters to their positions in G.
#' - G_temp: List of individual variance component matrices for each random effect.
#' - flist: List of grouping factors used in random effects.
#' - cnms: Column names of the random-effects design matrix.
#' - Ztlist: List of per-term transposed design matrices for random effects.
#' - nl: Number of levels for each grouping factor.
#' @importFrom Matrix t sparseMatrix drop0 forceSymmetric
#' @keywords internal
mkReTrms <- function(bars,
                     fr,
                     drop.unused.levels = TRUE,
                     reorder.terms = FALSE,
                     reorder.vars = FALSE) {
  barnames <- function(bars) vapply(bars, function(x) deparse1(x[[3]]), "")
  if (!length(bars))
    stop("No random effects terms specified in formula",
         call. = FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA),
            inherits(fr, "data.frame"))
  names(bars) <- barnames(bars)
  term.names <- vapply(bars, deparse1, "")
  blist <- lapply(bars, mkBlist, fr, drop.unused.levels, reorder.vars = reorder.vars)
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (reorder.terms) {
    if (any(diff(nl) > 0)) {
      ord <- rev(order(nl))
      blist <- blist[ord]
      nl <- nl[ord]
      term.names <- term.names[ord]
    }
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1))/2)
  nb <- nc * nl
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                 sum(nb), q))
  }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  G <- Matrix::t(do.call(
    Matrix::sparseMatrix,
    do.call(
      rbind,
      lapply(
        seq_along(blist),
        function(i) {
          mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
          dd <- diag(nc[i])
          ltri <- lower.tri(dd, diag = TRUE)
          ii <- row(dd)[ltri]
          jj <- col(dd)[ltri]
          data.frame(
            i = as.vector(mm[, ii]) + boff[i],
            j = as.vector(mm[, jj]) + boff[i],
            x = as.double(rep.int(seq_along(ii), rep.int(nl[i], length(ii))) + thoff[i])
          )
        }
      )
    )
  ))
  G <- Matrix::forceSymmetric(G, uplo = "U")

  idx <- c(0, cumsum(nth))[seq_along(nc)]
  G_temp <- lapply(seq_along(nc), function(i){
    mat_temp <- matrix(0, nc[[i]], nc[[i]])
    mat_temp[lower.tri(mat_temp, diag = TRUE)] <- seq_len(nth[[i]]) + idx[i]
    mat_temp[upper.tri(mat_temp)] <- t(mat_temp)[upper.tri(mat_temp)]
    rownames(mat_temp) <- colnames(mat_temp) <- cnms[[i]]
    return(mat_temp)
  })
  names(G_temp) <- names(cnms)

  ll <- list(Zt = Matrix::drop0(Zt), G_temp = G_temp, G = G, Gind = as.integer(G@x))
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll$nl <- nl
  ll
}

#' Random effects list
#'
#' Copied from lme4.
#'
#' @importFrom Matrix fac2sparse KhatriRao sparse.model.matrix
#' @importFrom methods is
#' @noRd
#' @keywords lme4 functions
mkBlist <- function(x, frloc, drop.unused.levels = TRUE, reorder.vars = FALSE) {
  replaceTerm <- function(term, target, repl) {
    inForm <- function(form, value) {
      if (any(sapply(form, identical, value)))
        return(TRUE)
      if (all(sapply(form, length) == 1))
        return(FALSE)
      return(any(vapply(form, inForm, value, FUN.VALUE = logical(1))))
    }
    if (identical(term, target))
      return(repl)
    if (!inForm(term, target))
      return(term)
    if (length(term) == 2) {
      return(
        substitute(
          OP(x),
          list(
            OP = replaceTerm(term[[1]], target, repl),
            x = replaceTerm(term[[2]], target, repl)
          )
        )
      )
    }
    return(
      substitute(
        OP(x, y),
        list(
          OP = replaceTerm(term[[1]], target, repl),
          x = replaceTerm(term[[2]], target, repl),
          y = replaceTerm(term[[3]], target, repl)
        )
      )
    )
  }
  `%i%` <- function(f1, f2, fix.order = TRUE) {
    if (!is.factor(f1) || !is.factor(f2)) stop("both inputs must be factors")
    f12 <- paste(f1, f2, sep = ":")
    ## explicitly specifying levels is faster in any case ...
    u <- which(!duplicated(f12))
    if (!fix.order) return(factor(f12, levels = f12[u]))
    ## deal with order of factor levels
    levs_rank <- length(levels(f2))*as.numeric(f1[u])+as.numeric(f2[u])
    return(factor(f12, levels = (f12[u])[order(levs_rank)]))
  }
  frloc <- factorize(x, frloc)
  ff0 <- replaceTerm(x[[3]], quote(`:`), quote(`%i%`))
  ff <- try(eval(substitute(factor(fac), list(fac = ff0)),
                 frloc), silent = TRUE)
  if (inherits(ff, "try-error")) {
    stop("couldn't evaluate grouping factor ", deparse1(x[[3]]),
         " within model frame:", "error =", c(ff), " Try adding grouping factor to data ",
         "frame explicitly if possible", call. = FALSE)
  }
  if (all(is.na(ff)))
    stop("Invalid grouping factor specification, ", deparse1(x[[3]]),
         call. = FALSE)
  if (drop.unused.levels)
    ff <- factor(ff, exclude = NA)
  nl <- length(levels(ff))
  has.sparse.contrasts <- function(x) {
    cc <- attr(x, "contrasts")
    !is.null(cc) && methods::is(cc, "sparseMatrix")
  }
  any.sparse.contrasts <- any(vapply(frloc, has.sparse.contrasts,
                                     FUN.VALUE = TRUE))
  mMatrix <- if (!any.sparse.contrasts)
    stats::model.matrix
  else Matrix::sparse.model.matrix
  mm <- mMatrix(eval(substitute(~foo, list(foo = x[[2]]))),
                frloc)
  # if (reorder.vars) {
  #   mm <- mm[colSort(colnames(mm)), ]
  # }
  sm <- Matrix::fac2sparse(ff, to = "d", drop.unused.levels = drop.unused.levels)
  sm <- Matrix::KhatriRao(sm, t(mm))
  dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

#' Residual Variance-Covariance Matrices
#'
#' @param data a data frame with grouping factors and covariates.
#' @param correlation a `corStruct` object created by [corClasses] functions.
#' If NULL, an identity matrix is assumed.
#' @return A list containing:
#' - corStruct: An initialized correlation structure.
#' - corframe: A processed data frame with indexed grouping variables and
#'   ordering for correlation structures.
#' - R: A block-diagonal residual variance-covariance structure, not yet scaled
#' by the residual variance
#' @import stats
#' @importFrom nlme corMatrix Initialize
#' @keywords internal
mkRTrms <- function(data,
                    correlation = NULL) {
  if (is.null(correlation))
    return(list(
      corStruct = NULL,
      corframe = data,
      R = Matrix::Diagonal(nrow(data))))

  corform <- attr(correlation, "formula")
  data <- data[, colnames(data) %in% all.vars(corform), drop = FALSE]
  grp.var <- strsplit(deparse1(corform[[2]]), " \\| ")[[1]][2]
  srt.var <- strsplit(deparse1(corform[[2]]), " \\| ")[[1]][1]
  grpform <- stats::reformulate(grp.var)
  corframe <- split(data, grpform, drop = TRUE)
  corframe <- lapply(corframe, function(df){
    df[, "idx_data"] <- as.integer(rownames(df))
    return(df)
  })

  if (any(class(correlation) %in% c("corSymm", "corAR1"))) {
    if (srt.var != 1) {
      corframe <- lapply(corframe, function(df){
        if (!is.integer(df[[srt.var]])) {
          df[[srt.var]] <- as.integer(df[[srt.var]])
        }
        df <- df[do.call(order, df[srt.var]), , drop = FALSE]
        return(df)
      })
    }
  }

  if (any(class(correlation) %in% c("corCAR1"))) {
    if (srt.var != 1) {
      corframe <- lapply(corframe, function(df){
        if (!is.numeric(df[[srt.var]])) {
          df[[srt.var]] <- as.numeric(df[[srt.var]])
        }
        df <- df[do.call(order, df[srt.var]), , drop = FALSE]
        return(df)
      })
    }
  }

  if (any(class(correlation) %in% c("corARMA"))) {
    if (srt.var != 1) {
      corframe <- lapply(corframe, function(df){
        if (!is.integer(df[[srt.var]])) {
          df[[srt.var]] <- as.integer(df[[srt.var]])
        }
        df <- df[do.call(order, df[srt.var]), , drop = FALSE]
        return(df)
      })
    }
  }

  if (any(class(correlation) %in% c("corExp", "corGaus", "corLin", "corRatio", "corSpher"))) {
    srt.var <- strsplit(gsub("\\s", "", srt.var), "\\+")[[1]]
    corframe <- lapply(corframe, function(df){
      fr_d <- subset(df, select = colnames(df) %in% srt.var)
      dist <- as.matrix(dist(fr_d))
      ord <- order(as.matrix(dist)[, 1])
      df <- df[ord, , drop = FALSE]
      return(df)
    })
  }

  corframe <- do.call(rbind, corframe)
  corframe[, "idx_R"] <- 1:nrow(corframe)
  corframe <- corframe[order(corframe$idx_R), ]

  corStruct <- nlme::Initialize(correlation, data = corframe)
  R_list <- nlme::corMatrix(corStruct)
  corframe <- corframe[order(corframe$idx_data), ]
  # FIXME: more robust way of matching positions in data and in R_list?
  rownames(corframe) <- corframe$idx_data
  R <- Matrix::.bdiag(R_list)
  R <- R[corframe$idx_R, corframe$idx_R]
  rTrms <- list(corStruct = corStruct, corframe = corframe, R = R)
  return(rTrms)
}

################################################################################
########         Model formula utilities copied from lme4               ########
################################################################################

#' Convert variables to factors as necessary
#'
#' @param x a formula
#' @param frloc a data frame
#' @param char.only (logical) convert only character variables to factors?
#' @return a copy of the data frame with factors converted
#' @keywords lme4 internal functions
factorize <- function(x, frloc, char.only = FALSE) {
  for (i in all.vars(x[[length(x)]])) {
    if (!is.null(curf <- frloc[[i]])) {
      if (!is.factor(curf) && (!char.only || is.character(curf))) {
        frloc[[i]] <- factor(curf)
      } else frloc[[i]] <- curf
    }
  }
  return(frloc)
}

#' Expand terms with \code{'||'} notation into separate \code{'|'} terms
#' @param term a formula
#' @return the modified term
#' @keywords lme4 internal functions
expandDoubleVerts <- function(term) {
  expandDoubleVert <- function(term) {
    frml <- formula(substitute(~x,list(x=term[[2]])))
    ## FIXME: do this without paste and deparse if possible!
    ## need term.labels not all.vars to capture interactions too:
    newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
    if(attr(terms(frml), "intercept")!=0)
      newtrms <- c("1", newtrms)

    as.formula(paste("~(",
                     paste(vapply(newtrms, function(trm)
                       paste0(trm, "|", deparse(term[[3]])), ""),
                       collapse=")+("), ")"))[[2]]
  }

  if (!is.name(term) && is.language(term)) {
    if (term[[1]] == as.name("(")) {
      term[[2]] <- expandDoubleVerts(term[[2]])
    }
    stopifnot(is.call(term))
    if (term[[1]] == as.name('||'))
      return( expandDoubleVert(term) )
    ## else :
    term[[2]] <- expandDoubleVerts(term[[2]])
    if (length(term) != 2) {
      if(length(term) == 3)
        term[[3]] <- expandDoubleVerts(term[[3]])
    }
  }
  term
}

#' Determine random-effects expressions from a formula
#'
#' @param term a mixed-model formula
#' @return pairs of expressions that were separated by vertical bars
#' @section Note: This function is called recursively on individual
#' terms in the model, which is why the argument is called \code{term} and not
#' a name like \code{form}, indicating a formula.
#' @keywords lme4 internal functions
#' @importFrom methods is
findbars <- function(term) {
  ## Recursive function applied to individual terms
  fb <- function(term) {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(fb(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(fb(term[[2]]))
    c(fb(term[[2]]), fb(term[[3]]))
  }
  ## Expand any slashes in the grouping factors returned by fb
  expandSlash <- function(bb) {
    ## Create the interaction terms for nested effects
    makeInteraction <- function(x) {
      if (length(x) < 2) return(x)
      trm1 <- makeInteraction(x[[1]])
      trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
      list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
    }
    ## Return the list of '/'-separated terms
    slashTerms <- function(x) {
      if (!("/" %in% all.names(x))) return(x)
      if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor",call.=FALSE)
      list(slashTerms(x[[2]]), slashTerms(x[[3]]))
    }

    if (!is.list(bb))
      expandSlash(list(bb))
    else
      unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
          ## lapply(unlist(...)) - unlist returns a flattened list
          lapply(unlist(makeInteraction(trms)),
                 function(trm) substitute(foo|bar, list(foo = x[[2]], bar = trm)))
        else x
      }))
  }

  modterm <- expandDoubleVerts(
    if(methods::is(term, "formula")) term[[length(term)]] else term)
  expandSlash(fb(modterm))
}

#' Substitute Bars
#'
#' Substitute the '+' function for the '|' and '||' function in a mixed-model
#' formula.
#'
#' @param term a mixed-model formula
#' @return the formula with all |  and || operators replaced by +
#' @section Note: This function is called recursively on individual
#' terms in the model, which is why the argument is called \code{term} and not
#' a name like \code{form}, indicating a formula.
#' @keywords lme4 internal functions
subbars <- function(term) {
  if (is.name(term) || !is.language(term)) return(term)
  if (length(term) == 2) {
    term[[2]] <- subbars(term[[2]])
    return(term)
  }
  stopifnot(length(term) >= 3)
  if (is.call(term) && term[[1]] == as.name('|'))
    term[[1]] <- as.name('+')
  if (is.call(term) && term[[1]] == as.name('||'))
    term[[1]] <- as.name('+')
  for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
  term
}

#' Omit terms separated by vertical bars in a formula
#'
#' Remove the random-effects terms from a mixed-effects formula,
#' thereby producing the fixed-effects formula.
#'
#' @param term the right-hand side of a mixed-model formula
#' @return the fixed-effects part of the formula
#' @section Note: This function is called recursively on individual
#' terms in the model, which is why the argument is called \code{term} and not
#' a name like \code{form}, indicating a formula.
#' @keywords lme4 internal functions
#' @importFrom methods is
nobars <- function(term) {
  e <- environment(term)
  nb <- nobars_(term)  ## call recursive version
  if (methods::is(term,"formula") && length(term)==3 && is.symbol(nb)) {
    ## called with two-sided RE-only formula:
    ##    construct response~1 formula
    nb <- stats::reformulate("1", response=deparse(nb))
  }
  ## called with one-sided RE-only formula, or RHS alone
  if (is.null(nb)) {
    nb <- if (methods::is(term,"formula")) ~1 else 1
  }
  environment(nb) <- e
  nb
}

nobars_ <- function(term) {
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if(is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

isBar <- function(term) {
  if(is.call(term)) {
    if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for(i in seq_along(term)) {
      if(isBar(term[[i]])) return(TRUE)
    }
  }
  FALSE
}

anyBars <- function(term) {
  any(c('|','||') %in% all.names(term))
}

#' Random Effects formula only
#' @noRd
#' @keywords lme4 internal functions
reOnly <- function(f, response=FALSE) {
  stats::reformulate(paste0("(", vapply(findbars(f), deparse1, ""), ")"),
              response = if(response && length(f)==3L) f[[2]],
              env = environment(f))
}

################################################################################
########    Data and model inspection utilities adapted from lme4       ########
################################################################################

#' check scale of non-dummy columns of X
#' @noRd
# TODO: really a useful function?
checkScaleX <- function(X, kind = "warning", tol = 1000) {
  cont.cols <- apply(X, 2, function(z) !all(z %in% c(0, 1)))
  col.sd <- apply(X[, cont.cols, drop = FALSE], 2L, sd)
  sdcomp <- outer(col.sd, col.sd, "/")
  logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
  logsd <- abs(log(col.sd))

  if (any(c(logcomp, logsd) > log(tol))) {
    wmsg <- "Some predictor variables are on very different scales:"
    wmsg <- paste(wmsg, "consider rescaling")
    warning(wmsg, call. = FALSE)
  }
  else wmsg <- character()
  structure(X, msgScaleX = wmsg)
}

#' ensure full rank design matrix
#' @importFrom Matrix rankMatrix
#' @noRd
chkRank.drop.cols <- function(X, kind = "stop.deficient", tol = 1e-07, method = "qr") {
  if (kind == "ignore") return(X)
  p <- ncol(X)
  if (kind == "stop.deficient") {
    if ((rX <- Matrix::rankMatrix(X, tol = tol, method = method)) < p)
      stop(
        gettextf(
          sub(
            "\n +",
            "\n",
            "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);\nthe fixed effects will be jointly unidentifiable"
          ),
          rX,
          p
        ),
        call. = FALSE
      )
  }
  else {
    qr.X <- qr(X, tol = tol, LAPACK = FALSE)
    rnkX <- qr.X$rank
    if (rnkX == p)
      return(X)
    msg <- sprintf(
      ngettext(
        p - rnkX,
        "fixed-effect model matrix is rank deficient so dropping %d column / coefficient",
        "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"
      ),
      p - rnkX
    )
    if (kind != "silent.drop.cols")
      (if (kind == "warn+drop.cols") warning
       else message)(msg, domain = NA)
    contr <- attr(X, "contrasts")
    asgn <- attr(X, "assign")
    keep <- qr.X$pivot[seq_len(rnkX)]
    dropped.names <- colnames(X[, -keep, drop = FALSE])
    X <- X[, keep, drop = FALSE]
    if (Matrix::rankMatrix(X, tol = tol, method = method) < ncol(X))
      stop(gettextf("Dropping columns failed to produce full column rank design matrix"),
           call. = FALSE)
    if (!is.null(contr))
      attr(X, "contrasts") <- contr
    if (!is.null(asgn))
      attr(X, "assign") <- asgn[keep]
    attr(X, "msgRankdrop") <- msg
    attr(X, "col.dropped") <- setNames(qr.X$pivot[(rnkX +
                                                     1L):p], dropped.names)
  }
  X
}

#' Check that grouping factors have at least 2 and </<= nobs(.) levels
#' @noRd
# TODO: necessity to keep flexibility as originally in lme4?
checkNlevels <- function(flist, n, allow.n = FALSE) {
  nlevelVec <- vapply(flist, function(x) nlevels(factor(x, exclude = NA)), 1)
  checkCtrlLevels("check.nlev.gtr.1", cc <- "stop")
  if (doCheck(cc) && any(nlevelVec < 2)) {
    wstr <- "grouping factors must have > 1 sampled level"
    switch(
      cc,
      warning = warning(wstr, call. = FALSE),
      stop = stop(wstr, call. = FALSE),
      stop(gettextf("unknown check level for '%s'", "check.nlev.gtr.1"), domain = NA)
    )
  } else wstr <- character()
  checkCtrlLevels("check.nobs.vs.nlev", cc <- "stop")
  if (doCheck(cc) && any(if (allow.n) nlevelVec > n else nlevelVec >= n)) {
    w <- if (allow.n) which(nlevelVec > n) else which(nlevelVec >= n)
    bad_facs <- names(nlevelVec)[w]
    wst2 <- gettextf("number of levels of each grouping factor must be %s number of observations",
                     if (allow.n) "<=" else "<")
    wst2 <- paste0(wst2, " (problems: ", paste(bad_facs, collapse = ", "), ")")
    switch(cc,
           warning = warning(wst2, call. = FALSE),
           stop = stop(wst2, call. = FALSE))
  } else wst2 <- character()
  checkCtrlLevels("check.nlev.gtreq.5", cc <- "ignore")
  if (doCheck(cc) && any(nlevelVec < 5)) {
    wst3 <- "grouping factors with < 5 sampled levels may give unreliable estimates"
    switch(
      cc,
      warning = warning(wst3, call. = FALSE),
      stop = stop(wst3, call. = FALSE),
      stop(gettextf("unknown check level for '%s'", "check.nlev.gtreq.5"), domain = NA)
    )
  } else wst3 <- character()
  c(wstr, wst2, wst3)
}

#' For each r.e. term, test if Z has more columns than rows to detect unidentifiability
#' @noRd
# TODO: necessity to keep flexibility as originally in lme4?
checkZdims <- function(Ztlist, n, allow.n = FALSE) {
  stopifnot(is.list(Ztlist), is.numeric(n))
  checkCtrlLevels("check.nobs.vs.nRE", cc <- "stop")
  term.names <- names(Ztlist)
  rows <- vapply(Ztlist, nrow, 1L)
  cols <- vapply(Ztlist, ncol, 1L)
  stopifnot(all(cols == n))
  if (doCheck(cc)) {
    unique(unlist(lapply(seq_along(Ztlist), function(i) {
      ww <- wmsg(cols[i], rows[i], allow.n, "number of observations",
                 "number of random effects", sprintf(" for term (%s)", term.names[i]))
      if (ww$unident) {
        switch(cc, warning = warning(ww$wstr, call. = FALSE),
               stop = stop(ww$wstr, call. = FALSE),
               stop(gettextf("unknown check level for '%s'", "check.nobs.vs.nRE"), domain = NA))
        ww$wstr
      } else character()
    })))
  } else character()
}

#' check rank of Z
#' @noRd
#' @importFrom Matrix rankMatrix
checkZrank <- function(Zt, n, nonSmall = 1e+06, allow.n = FALSE) {
  if (doCheck(cc <- "ignore")) {
    checkCtrlLevels("check.nobs.vs.rankZ", cc, smallOK = TRUE)
    d <- dim(Zt)
    doTr <- d[1L] < d[2L]
    if (!(grepl("Small", cc) && prod(d) > nonSmall)) {
      rankZ <- Matrix::rankMatrix(if (doTr) t(Zt) else Zt, method = "qr")
      ww <- wmsg(n, rankZ, allow.n, "number of observations", "rank(Z)")
      if (is.na(rankZ)) {
        cc <- "stop"
        ww <- list(
          unident = TRUE,
          wstr = sub("^.*;", "rank(Z) is NA: invalid random effect factors?", ww$wstr)
        )
      }
      if (ww$unident) {
        switch(
          cc,
          warningSmall = ,
          warning = warning(ww$wstr, call. = FALSE),
          stopSmall = ,
          stop = stop(ww$wstr, call. = FALSE),
          stop(gettextf("unknown check level for '%s'", "check.nobs.vs.rankZ"), domain = NA))
        ww$wstr
      } else character()
    } else character()
  } else character()
}

checkCtrlLevels <- function(cstr, val, smallOK = FALSE) {
  bvals <- c("message", "warning", "stop", "ignore")
  if (smallOK)
    bvals <- outer(bvals, c("", "Small"), paste0)
  if (!is.null(val) && !val %in% bvals)
    stop("invalid control level ", sQuote(val), " in ",
         cstr, ": valid options are {", paste(sapply(bvals, sQuote), collapse = ","), "}")
  invisible(NULL)
}

doCheck <- function(x) is.character(x) && !any(x == "ignore")

wmsg <- function(n, cmp.val, allow.n, msg1 = "", msg2 = "", msg3 = "") {
  if (allow.n) {
    unident <- n < cmp.val
    cmp <- "<"
    rstr <- ""
  } else {
    unident <- n <= cmp.val
    cmp <- "<="
    rstr <- " and the residual variance (or scale parameter)"
  }
  wstr <- sprintf("%s (=%d) %s %s (=%d)%s; the random-effects parameters%s are probably unidentifiable",
                  msg1, n, cmp, msg2, cmp.val, msg3, rstr)
  list(unident = unident, wstr = wstr)
}

################################################################################
########         (co-)variance parameters related utilities             ########
################################################################################

#' Compute Variance-Covariance Matrix of Random Effects with Respect to
#' Variance Components
#'
#' Used for computing numerical derivatives.
#' @noRd
#' @keywords internal
vcovb_vp <- function(vcomp = NULL, reTrms = NULL) {
  if (length(vcomp) == 0 || is.null(reTrms)) return(NULL)
  G <- reTrms$G
  G@x <- vcomp[reTrms$Gind]
  return(G)
}

#' Compute Variance-Covariance Matrix of Residuals with Respect to Residual Variance
#' and Correlations
#'
#' Used for computing numerical derivatives.
#' @noRd
#' @keywords internal
vcove_vp <- function(corr = NULL, sigma2, rTrms) {
  if (is.null(rTrms$corStruct) || length(corr) == 0) {
    R <- Matrix::Diagonal(nrow(rTrms$corframe))
    return(sigma2 * R)
  }

  # Get correlation structure info
  corStruct <- rTrms$corStruct
  data <- rTrms$corframe
  cor_class <- class(corStruct)[1]
  cor_formula <- attr(corStruct, "formula")

  # Recreate correlation structure with new parameter value
  correlation <- switch(cor_class,
                        "corAR1" = nlme::corAR1(value = corr, form = cor_formula),
                        "corARMA" = {
                          p <- attr(corStruct, "p")
                          q <- attr(corStruct, "q")
                          nlme::corARMA(value = corr, form = cor_formula, p = p, q = q)
                        },
                        "corCAR1" = nlme::corCAR1(value = corr, form = cor_formula),
                        "corCompSymm" = nlme::corCompSymm(value = corr, form = cor_formula),
                        "corSymm" = nlme::corSymm(value = corr, form = cor_formula),
                        "corExp" = {
                          metric <- attr(corStruct, "metric")
                          nugget <- attr(corStruct, "nugget")
                          nlme::corExp(value = corr, form = cor_formula,
                                       metric = if (is.null(metric)) "euclidean" else metric,
                                       nugget = if (is.null(nugget)) FALSE else nugget)
                        },
                        "corGaus" = {
                          metric <- attr(corStruct, "metric")
                          nugget <- attr(corStruct, "nugget")
                          nlme::corGaus(value = corr, form = cor_formula,
                                        metric = if (is.null(metric)) "euclidean" else metric,
                                        nugget = if (is.null(nugget)) FALSE else nugget)
                        },
                        "corLin" = {
                          metric <- attr(corStruct, "metric")
                          nugget <- attr(corStruct, "nugget")
                          nlme::corLin(value = corr, form = cor_formula,
                                       metric = if (is.null(metric)) "euclidean" else metric,
                                       nugget = if (is.null(nugget)) FALSE else nugget)
                        },
                        "corRatio" = {
                          metric <- attr(corStruct, "metric")
                          nugget <- attr(corStruct, "nugget")
                          nlme::corRatio(value = corr, form = cor_formula,
                                         metric = if (is.null(metric)) "euclidean" else metric,
                                         nugget = if (is.null(nugget)) FALSE else nugget)
                        },
                        "corSpher" = {
                          metric <- attr(corStruct, "metric")
                          nugget <- attr(corStruct, "nugget")
                          nlme::corSpher(value = corr, form = cor_formula,
                                         metric = if (is.null(metric)) "euclidean" else metric,
                                         nugget = if (is.null(nugget)) FALSE else nugget)
                        },
                        stop(sprintf("Unsupported correlation structure: %s", cor_class))
  )

  R <- mkRTrms(data, correlation)$R
  R <- sigma2 * R
  return(R)
}

#' Compute Variance-Covariance Matrix of y with Respect to Variance Parameters
#'
#' Used for computing numerical derivatives.
#' @noRd
#' @keywords internal
vcovy_vp <- function(varpar, reTrms, rTrms) {

  na <- attr(varpar, "vcomp")
  nb <- attr(varpar, "corr")
  nc <- attr(varpar, "sigma2")

  vcomp <- varpar[seq_len(na)]
  corr <- varpar[seq.int(na + 1, length.out = nb)]
  sigma2 <- varpar[seq.int(na + nb + 1, length.out = nc)]

  G <- vcovb_vp(vcomp, reTrms)
  R <- vcove_vp(corr, sigma2, rTrms)
  if (is.null(G) || is.null(reTrms))
    V <- R else V <- Matrix::t(reTrms$Zt) %*% G %*% reTrms$Zt + R
  # return(V)
  return(as.matrix(V)) # FIXME, matrix class
}

#' Compute Variance-Covariance Matrix of beta with Respect to Variance Parameters
#'
#' Used for computing numerical derivatives.
#' @noRd
#' @keywords internal
vcovbeta_vp <- function(varpar, fxTrms, reTrms, rTrms) {
  V <- vcovy_vp(varpar, reTrms, rTrms)
  V_inv <- Matrix::solve(V)
  C_inv = t(fxTrms$X) %*% V_inv %*% fxTrms$X
  C = Matrix::solve(C_inv)
  return(C)
  # return(as.matrix(C)) # FIXME, matrix class
}

#' Compute REML or ML Information Matrix
#' @importFrom numDeriv jacobian
#' @noRd
#' @keywords internal
informat <- function(varpar, fxTrms, reTrms, rTrms, REML) {
  V <- vcovy_vp(varpar = varpar, reTrms, rTrms)
  V_inv <- Matrix::solve(V)
  if (REML)
    P <- V_inv - V_inv %*% fxTrms$X %*% solve(t(fxTrms$X) %*% V_inv %*% fxTrms$X) %*% t(fxTrms$X) %*% V_inv
  n_params <- length(varpar)
  info_mat <- matrix(0, nrow = n_params, ncol = n_params)
  dV_dvarpar <- numDeriv::jacobian(func = vcovy_vp, x = varpar, reTrms = reTrms, rTrms = rTrms)
  dV_dvarpar <- lapply(1:ncol(dV_dvarpar), function(i)
    array(dV_dvarpar[, i], dim=rep(nrow(V), 2)))
  if (REML) {
    for (i in 1:n_params) {
      for (j in i:n_params) {
        trace_ij <- sum(diag(P %*% dV_dvarpar[[i]] %*% P %*% dV_dvarpar[[j]]))
        info_mat[i, j] <- 0.5 * trace_ij
        info_mat[j, i] <- info_mat[i, j]
      }
    }
    return(info_mat)
  }
  for (i in 1:n_params) {
    for (j in i:n_params) {
      trace_ij <- sum(diag(V_inv %*% dV_dvarpar[[i]] %*% V_inv %*% dV_dvarpar[[j]]))
      info_mat[i, j] <- 0.5 * trace_ij
      info_mat[j, i] <- info_mat[i, j]
    }
  }
  return(info_mat)
}

################################################################################
########      F-test Contrast matrices, adapted from lmerTest           ########
################################################################################

## Contrast coding for factors is restricted to treatment contrasts in earlier
## processing steps for constructing design matrices.

#' Compute contrast matrices for type I tests
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_type1 <- function(fixedfr) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  p <- ncol(X) # number of parameters
  if(p == 0L) return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  if(p == 1L && attr(terms, "intercept")) # intercept-only model
    return(list(matrix(numeric(0L), ncol=1L)))
  # Compute 'normalized' doolittle factorization of XtX:
  L <- if(p == 1L) matrix(1L) else t(doolittle(crossprod(X))$L)
  dimnames(L) <- list(colnames(X), colnames(X))
  # Determine which rows of L belong to which term:
  ind.list <- term2colX(terms, X)[attr(terms, "term.labels")]
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
}

#' Create the 'genuine type II contrast' for all terms that are
#' contained in other terms.
#'
#' For all terms which are not contained in other
#' terms, the simple marginal contrast is computed.
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_type2_unfolded <- function(fixedfr, which=NULL) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(fixedfr))
  } else stopifnot(is.character(which), all(which %in% term_names))

  is_contained <- containment(terms)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  do_type2 <- setdiff(term_names, do_marginal)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(fixedfr, which=do_marginal)
  if(length(do_type2))
    Llist <- c(Llist, get_contrasts_type2(fixedfr, which=do_type2))
  Llist[term_names]
}

#' Computes the type 2 contrasts - either for all terms or for those
#' included in 'which' (a character vector naming model terms).
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_type2 <- function(fixedfr, which=NULL) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  data_classes <- attr(terms, "dataClasses") # all vars
  if(is.null(asgn <- attr(X, "assign")))
    stop("design matrix 'X' should have a non-null 'assign' attribute")
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(fixedfr))
  } else stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  # Compute containment:
  is_contained <- containment(terms)
  # Compute term asignment list: map from terms to columns in X
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if(has_intercept) c("(Intercept)", term_names)[asgn + 1] else
    term_names[asgn[asgn > 0]]
  if(!length(col_terms) == ncol(X)) # should never happen.
    stop("An error happended when computing Type II contrasts")
  term2colX <- split(seq_along(col_terms), col_terms)[unique(col_terms)]
  # Compute contrast for each term - return as named list:
  lapply(which, function(term) {
    # Reorder the cols in X to [, unrelated_to_term, term, contained_in_term]
    cols_term <- unlist(term2colX[c(term, is_contained[[term]])])
    Xnew <- cbind(X[, -cols_term, drop=FALSE], X[, cols_term, drop=FALSE])
    # Compute order of terms in Xnew:
    newXcol_terms <- c(col_terms[-cols_term], col_terms[cols_term])
    # Compute Type I contrasts for the reordered X:
    Lc <- t(doolittle(crossprod(Xnew))$L)
    dimnames(Lc) <- list(colnames(Xnew), colnames(Xnew))
    # Extract rows for term and get original order of columns:
    Lc[newXcol_terms == term, colnames(X), drop=FALSE]
  })
}

#' Create contrast matrices for type III tests
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_type3 <- function(fixedfr) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  term_names <- attr(terms, "term.labels")
  # If model has at most one term return Type I contrasts:
  if(ncol(X) <= 1L || length(term_names) <= 1L)
    return(get_contrasts_type1(fixedfr))
  # Get 'complete' design matrix:
  rdX <- get_rdX(fixedfr, do.warn = TRUE) # treatment contrasts
  # cols for aliased coefs should be removed in X; not in rdX.
  # This makes ginv(X) unique!
  L <- zapsmall(t(MASS::ginv(X) %*% rdX)) # basic contrast matrix
  dimnames(L) <- list(colnames(rdX), colnames(X))
  # Orthogonalize contrasts for terms which are contained in other terms:
  map <- term2colX(terms, X)
  is_contained <- containment(terms)
  # Orthogonalize higher order terms before lower order terms:
  terms_order <- attr(terms, "order")
  orthog_order <- term_names[order(terms_order, decreasing = TRUE)]
  for(term in orthog_order) {
    # if term is contained in other terms:
    if(length(contains <- is_contained[[term]]) > 0) {
      # orthogonalize cols in L for 'term' wrt. cols that contain 'term':
      L[, map[[term]]] <-
        zapsmall(resid(lm.fit(x=L[, unlist(map[contains]), drop=FALSE],
                              y=L[, map[[term]], drop=FALSE])))
    }
  }
  # Keep rows in L corresponding to model coefficients:
  L <- L[colnames(X), , drop=FALSE]
  # Extract list of contrast matrices from L - one for each term:
  Llist <- lapply(map[term_names], function(term) t(L[, term, drop=FALSE]))
  # Keep all non-zero rows:
  lapply(Llist, function(L) L[rowSums(abs(L)) > 1e-8, , drop=FALSE])
}



#' Compute marginal contrast matrices.
#'
#' No tests of conformity with coefficients are implemented
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_marginal <- function(fixedfr, which=NULL) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(fixedfr))
  } else stopifnot(is.character(which), all(which %in% term_names))
  ## FIXME: test use of 'which' arg.

  # Compute map from terms to columns in X and contrasts matrix
  term2colX <- term2colX(terms, X)
  L <- structure(diag(ncol(X)), dimnames = list(colnames(X), colnames(X)))

  # Extract contrast for each term - return as named list:
  which <- setNames(as.list(which), which)
  lapply(which, function(term) {
    L[term2colX[[term]], , drop=FALSE]
  })
}

term2colX <- function(terms, X) {
  if (is.null(asgn <- attr(X, "assign")))
    stop("Invalid design matrix:", "design matrix 'X' should have a non-null 'assign' attribute",
         call. = FALSE)
  term_names <- attr(terms, "term.labels")
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if (has_intercept)
    c("(Intercept)", term_names)[asgn + 1]
  else term_names[asgn[asgn > 0]]
  if (!length(col_terms) == ncol(X))
    stop("An error happended when mapping terms to columns of X")
  nm <- union(unique(col_terms), term_names)
  res <- lapply(setNames(as.list(nm), nm), function(x) numeric(0L))
  map <- split(seq_along(col_terms), col_terms)
  res[names(map)] <- map
  res[nm]
}

#' Compute yates contrast matrices.
#' @noRd
#' @keywords lmerTest internal functions
get_contrasts_yates <- function(fixedfr) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  # Is this really type 4?
  term_names <- attr(terms, "term.labels")
  is_contained <- containment(terms)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  not_marginal <- setdiff(term_names, do_marginal)
  # Split not_marginal in do_yates and do_type2:
  do_yates <- need_yates(terms)
  do_type2 <- setdiff(not_marginal, do_yates)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(fixedfr, which=do_marginal)
  if(length(do_yates))
    Llist <- c(Llist, get_yates_contrast(fixedfr, which=do_yates))
  if(length(do_type2)) {
    data_classes <- attr(terms, "dataClasses")
    Llist <- c(Llist, get_contrasts_type2(fixedfr, which=do_type2))
  }
  Llist[term_names]
}

get_yates_contrast <- function(fixedfr, which=NULL) {
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) which <- term_names
  stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  var_list <- get_var_list(fixedfr)
  grid <- get_min_data(fixedfr)
  form <- formula(terms)
  coef_nm <- colnames(X)
  uX <- stats::model.matrix(form, data=grid)
  # Compute LS-means contrast:
  Llist <- lapply(which, function(term) {
    Lt <- stats::model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
    wts <- 1/colSums(Lt) # Yates' weights
    # Lt * c(Lt %*% wts)
    # L <- diag(wts) %*% t(Lt)
    L <- t(sweep(Lt, 2, wts, "*"))
    L %*% uX
  })
  # Check estimability:
  XX <- stats::model.matrix(terms, data=fixedfr)
  # Restore contrast coding here.
  nullspaceX <- nullspace(XX)
  not_estim <- sapply(Llist, function(L)
    any(!is_estimable(L, nullspace = nullspaceX)))
  if(any(not_estim))
    warning(sprintf("Yates contrast is not uniquely defined for: %s",
                    paste(names(Llist[not_estim]), collapse = ", ")),
            call. = FALSE)
  # Make contrast for joint test of contrast among LS-means:
  lapply(Llist, function(L) {
    (t(get_trts(rownames(L))) %*% L)[, coef_nm, drop=FALSE]
  })
}

need_yates <- function(terms) {
  ## Do not need yates for:
  ## - continuous variables
  ## - factors that are not contained in other factors
  ## Need yates for all other terms, i.e. terms which are:
  ##  - contained in other terms, AND
  ##  - which are not numeric/continuous
  term_names <- attr(terms, "term.labels")
  cont <- containment(terms)
  is_contained <- names(cont[sapply(cont, function(x) length(x) > 0)])
  nmt <- numeric_terms(terms)
  num_terms <- names(nmt[nmt])
  term_names[!term_names %in% num_terms &
               term_names %in% is_contained]
}

get_trts <- function(levs) {
  nlevs <- length(levs)
  ans <- t(cbind(-1, diag(nlevs - 1)))
  rownames(ans) <- levs
  colnames(ans) <- paste(levs[-1], levs[1], sep = " - ")
  ans
}

#' Determine the Containment Structure for All Terms in a Model
#' @noRd
#' @keywords lmerTest internal functions
containment <- function(terms) { # lm or merMod
  # For all terms 'T' in object compute the terms
  # Return a list:
  #   for each term 'T' a vector of terms that contain 'T'.
  data_classes <- attr(terms, "dataClasses")
  term_names <- attr(terms, "term.labels")
  factor_mat <- attr(terms, "factors")
  lapply(setNames(term_names, term_names), function(term) {
    term_names[term_contain(term, factor_mat, data_classes, term_names)]
  })
}

term_contain <- function(term, factors, dataClasses, term_names) {
  get_vars <- function(term) rownames(factors)[factors[, term] ==
                                                 1]
  contain <- function(F1, F2) {
    all(vars[[F1]] %in% vars[[F2]]) && length(setdiff(vars[[F2]],
                                                      vars[[F1]])) > 0L && setequal(numerics[[F1]], numerics[[F2]])
  }
  vars <- lapply(setNames(term_names, term_names), get_vars)
  numerics <- lapply(vars, function(varnms) varnms[which(dataClasses[varnms] ==
                                                           "numeric")])
  sapply(term_names, function(term_nm) contain(term, term_nm))
}

#' Compute rank-deficient design-matrix X using contr.treatment coding.
#' @noRd
#' @keywords lmerTest internal functions
get_rdX <- function(fixedfr, do.warn=TRUE) {
  terms <- terms(fixedfr)
  term_names <- attr(terms, "term.labels")
  # Compute rank-deficient (full) design-matrix, X:
  rdXi <- if(length(term_names)) lapply(term_names, function(trm) {
    form <- as.formula(paste0("~ 0 + ", trm))
    stats::model.matrix(form, data = fixedfr) # no contrast arg
  }) else list(stats::model.matrix(~ 1, data = fixedfr)[, -1, drop=FALSE])
  rdX <- do.call(cbind, rdXi)
  param_names <- unlist(lapply(rdXi, colnames))
  # Potentially add intercept:
  has_intercept <- attr(terms, "intercept") != 0
  if(has_intercept) {
    rdX <- cbind('(Intercept)'=rep(1, nrow(rdX)), rdX)
    param_names <- c("(Intercept)", param_names)
  }
  colnames(rdX) <- param_names
  # Warn/message if there are cells without data:
  is_zero <- which(apply(rdX, 2, function(x) all(x == 0)))
  if(do.warn && length(is_zero)) {
    txt <- sprintf("Missing cells for: %s. ",
                   paste(param_names[is_zero], collapse = ", "))
    # warning(paste(txt, "\nInterpret type III hypotheses with care."), call.=FALSE)
    message(paste(txt, "\nInterpret type III hypotheses with care."))
  }
  rdX
}

#' Doolittle Decomposition
#' @noRd
#' @keywords lmerTest internal functions
doolittle <- function(x, eps = 1e-06) {
  if (!is.matrix(x) || ncol(x) != nrow(x) || !is.numeric(x))
    stop("argument 'x' should be a numeric square matrix")
  stopifnot(ncol(x) > 1L)
  n <- nrow(x)
  L <- U <- matrix(0, nrow = n, ncol = n)
  diag(L) <- rep(1, n)
  for (i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for (j in 1:n) {
      U[i, j] <- x[i, j]
      if (im1 > 0) {
        for (k in 1:im1) {
          U[i, j] <- U[i, j] - L[i, k] * U[k, j]
        }
      }
    }
    if (ip1 <= n) {
      for (j in ip1:n) {
        L[j, i] <- x[j, i]
        if (im1 > 0) {
          for (k in 1:im1) {
            L[j, i] <- L[j, i] - L[j, k] * U[k, i]
          }
        }
        L[j, i] <- if (abs(U[i, i]) < eps)
          0
        else L[j, i]/U[i, i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list(L = L, U = U)
}


#' Estimability of Contrasts
#' @noRd
#' @keywords lmerTest internal functions
is_estimable <- function(contrast,
                         nullspace = NULL,
                         X = NULL,
                         tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(contrast))
    contrast <- matrix(contrast, ncol = length(contrast))
  N <- if (!is.null(nullspace)) {
    nullspace
  }
  else if (!is.null(X)) {
    nullspace(X)
  }
  else {
    stop("Need non-null 'nullspace' or 'X' to compute estimability")
  }
  if (ncol(contrast) != nrow(N))
    stop(sprintf("'contrast' has %i columns: expecting %i columns",
                 ncol(contrast), nrow(N)))
  res <- if (length(N) == 0)
    rep(TRUE, nrow(contrast))
  else c(abs(rowSums(contrast %*% N)) < tol)
  setNames(res, rownames(contrast))
}

#' Nullspace
#' @noRd
#' @keywords lmerTest internal functions
nullspace <- function(A, type = c("right", "left"), tol = sqrt(.Machine$double.eps)) {
  type <- match.arg(type)
  if (type == "left")
    return(nullspace(t(A), type = "right", tol = tol))
  if (length(A) == 0L)
    return(matrix(numeric(0L)))
  svdA <- svd(A, nv = ncol(A))
  tol <- 1e-08
  positive <- svdA$d > max(tol * svdA$d[1L], 0)
  rank <- sum(positive)
  set <- if (rank == 0)
    1:ncol(A)
  else -(1:rank)
  svdA$v[, set, drop = FALSE]
}

################################################################################
########                     t-test contrast matrices                   ########
################################################################################

#' Generate contrast matrix for acquiring factor means from model coefficients
#'
#' Adapted from lmerTest:::lsmeans_contrasts. For numerical variables, means identical to model coefficients.
#' @noRd
lsmeans_contrasts <- function(fixedfr, which=NULL, by = NULL) {
  if(is.null(which) | length(which) > 1)
    stop("exaclty one factor must be specified")
  which <- sub("\\*", ":", which)

  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)

  formula <- formula(terms)
  term_names <- attr(terms, "term.labels")
  factor_terms <- term_names[!numeric_terms(terms)]
  numeric_terms <- term_names[numeric_terms(terms)]

  # Contrasts
  Contr <- attr(X, "contrasts")
  var_names <- names(get_var_list(fixedfr))
  factor_mat <- attr(terms, "factors")
  vars_in_term <- var_names[factor_mat[var_names, which] == 1]

  # Get minimal 'unique rows' data:
  grid <- get_min_data(fixedfr)
  Contrasts <- .getXlevels(terms, grid)
  Contrasts[] <- "contr.treatment"

  if (!is.null(by)) { # conditional case
    gridlist <- split(grid, grid[[by]])
    names(gridlist) <- paste0(by, levels(grid[[by]]))
    Llist <- lapply(gridlist, function(grid){
      if (which %in% factor_terms) {
        uX <- stats::model.matrix(formula, data=grid, contrasts.arg=Contr)
        Lt <- stats::model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                           contrasts.arg=Contrasts[vars_in_term])
      }
      if (which %in% numeric_terms) {
        numvar_names <- names(get_num_list(fixedfr))
        numvars_in_term <- numvar_names[factor_mat[numvar_names, which] == 1]
        facvars_in_term <- setdiff(vars_in_term, numvar_names)
        grid[[numvars_in_term]] <- 1
        grid0 <- grid
        grid0[[numvars_in_term]] <- 0
        uX <- stats::model.matrix(formula, data=grid, contrasts.arg=Contr[facvars_in_term])
        uX0 <- stats::model.matrix(formula, data=grid0, contrasts.arg=Contr[facvars_in_term])
        uX <- uX - uX0
        Lt <- stats::model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                           contrasts.arg=Contrasts[facvars_in_term])
      }
      wts <- 1/colSums(Lt)
      L <- t(sweep(Lt, 2, wts, "*"))
      return(L %*% uX)
    })
    return(Llist)
  } # marginal case continues
  if (which %in% factor_terms) {
    uX <- stats::model.matrix(formula, data=grid, contrasts.arg=Contr)
    Lt <- stats::model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                       contrasts.arg=Contrasts[vars_in_term])
  }
  if (which %in% numeric_terms) {
    numvar_names <- names(get_num_list(fixedfr))
    numvars_in_term <- numvar_names[factor_mat[numvar_names, which] == 1]
    facvars_in_term <- setdiff(vars_in_term, numvar_names)
    grid[, numvars_in_term] <- 1
    grid0 <- grid
    grid0[, numvars_in_term] <- 0
    uX <- stats::model.matrix(formula, data=grid, contrasts.arg=Contr[facvars_in_term])
    uX0 <- stats::model.matrix(formula, data=grid0, contrasts.arg=Contr[facvars_in_term])
    uX <- uX - uX0
    Lt <- stats::model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                       contrasts.arg=Contrasts[facvars_in_term])
  }
  wts <- 1/colSums(Lt)
  L <- t(sweep(Lt, 2, wts, "*"))
  return(L %*% uX)
}

#' Generate contrast matrix to convert model coefficients into means for all predictors
#'
#' For numerical predictors, regression coefficients identical to means.
#' @importFrom MASS ginv
#' @noRd
means2beta_contrasts <- function(fixedfr){
  terms <- terms(fixedfr)
  X <- stats::model.matrix(terms, fixedfr)
  term_names <- attr(terms, "term.labels")
  term_names <- setNames(term_names, term_names)
  factor_terms <- term_names[!numeric_terms(terms)]
  numeric_terms <- term_names[numeric_terms(terms)]
  factor_mat <- attr(terms, "factors")
  var_names <- names(get_var_list(fixedfr))
  numvar_names <- names(get_num_list(fixedfr))
  ordn_terms <- term_names[attr(terms, "order") > 1]
  ordn_terms <-  setNames(ordn_terms, ordn_terms)
  dropvars <- lapply(ordn_terms, function(trm){
    if (trm %in% factor_terms) {
      drop_vars <- factor_mat[var_names, trm] == 1
    } else {
      if (all(var_names[factor_mat[var_names, trm] == 1] %in% numvar_names))
        dropvars <- logical(length(var_names)) else
          drop_vars <- factor_mat[var_names, trm] == 1 & var_names %in% numvar_names
    }
  })
  dropvars <- Reduce("|", dropvars)
  dropvars <- names(dropvars[dropvars])
  which <- setdiff(term_names, dropvars)
  which <- setNames(which, which)
  # Check if 2nd-order terms should be dropped
  # FIXME: add flexibility for more higher orders
  if (any(attr(terms, "order") > 2)) {
    max_order <- max(attr(terms, "order"))
    for (i in seq(3, max_order, by = 1)) {
      ordh_terms <- term_names[attr(terms, "order") == i]
      ordl_terms <- term_names[attr(terms, "order") == i-1]
      is_contained <- function(x, y) {
        x_parts <- unlist(strsplit(x, ":", fixed = TRUE))
        y_parts <- unlist(strsplit(y, ":", fixed = TRUE))
        all(x_parts %in% y_parts)
      }
      factor_matl <- outer(ordl_terms, ordh_terms, Vectorize(is_contained))
      droptrms <- lapply(ordh_terms, function(trm){
        if (trm %in% factor_terms) {
          drop_trms <- factor_matl[, trm] == 1
        } else {
          if (all(ordl_terms[factor_matl[, trm] == 1] %in% numeric_terms))
            drop_trms <- logical(length(ordl_terms)) else
              drop_trms <- factor_matl[, trm] == 1 & ordl_terms %in% numeric_terms
        }
      })
      droptrms <- Reduce("|", droptrms)
      droptrms <- names(droptrms[droptrms])
      which <- setdiff(which, droptrms)
    }
  }
  # if (any(attr(terms, "order") > 2)) {
  #   ord3_terms <- term_names[attr(terms, "order") == 3]
  #   ord2_terms <- term_names[attr(terms, "order") == 2]
  #   is_contained <- function(x, y) {
  #     x_parts <- unlist(strsplit(x, ":", fixed = TRUE))
  #     y_parts <- unlist(strsplit(y, ":", fixed = TRUE))
  #     all(x_parts %in% y_parts)
  #   }
  #   factor_mat2 <- outer(ord2_terms, ord3_terms, Vectorize(is_contained))
  #   droptrms <- lapply(ord3_terms, function(trm){
  #     if (trm %in% factor_terms) {
  #       drop_trms <- factor_mat2[, trm] == 1
  #     } else {
  #       if (all(ord2_terms[factor_mat2[, trm] == 1] %in% numeric_terms))
  #         drop_trms <- logical(length(ord2_terms)) else
  #           drop_trms <- factor_mat2[, trm] == 1 & ord2_terms %in% numeric_terms
  #     }
  #   })
  #   droptrms <- Reduce("|", droptrms)
  #   droptrms <- names(droptrms[droptrms])
  #   which <- setdiff(which, droptrms)
  # }
  Llist <- lapply(which, function(wic){
    lsmeans_contrasts(fixedfr, wic)
  })
  L <- do.call(rbind, Llist)
  if(colSums(L)[1] == 0) L <- rbind(`(Intercept)` = c(1, numeric(ncol(L)-1)), L)
  return(L)
}

#' Convert Model Coefficients to Means and Vice Versa
#'
#' @importFrom MASS ginv
#' @noRd
means2beta <- function(L, means = NULL, beta = NULL){
  if (is.null(means) && is.null(beta)) {
    stop("at least one of 'means' and 'beta' must be provided")
  }
  if (!is.null(means) && !is.null(beta)) {
    beta2 <- drop(MASS::ginv(L) %*% drop(means))
    if (!identical(beta2, drop(beta))) {
      stop("'means' and 'beta' are not aligned")
    }
  }
  if (!is.null(means) && is.null(beta)) {
    means <- drop(means)
    beta <- drop(MASS::ginv(L) %*% means)
  }
  if (is.null(means) && !is.null(beta)) {
    beta <- drop(beta)
    means <- drop(L %*% beta)
  }
  names(beta) <- colnames(L)
  names(means) <- rownames(L)
  return(list(beta = beta, means = means))
}

################################################################################
########            Term utilities adapted from lmerTest                ########
################################################################################

#' Determines for all terms (not just all variables) if the 'dataClass'
#' is numeric.
#'
#' Interactions involving one or more numerics variables are numeric.
#'
#' @noRd
numeric_terms <- function(terms) {
  all_vars <- sapply(attr(terms, "variables")[-1], deparse)
  data_classes <- attr(terms, "dataClasses")
  var_class <- data_classes[names(data_classes) %in% all_vars]
  # FIXME: should coerce character to factor?
  factor_vars <- names(var_class[var_class %in% c("factor", "ordered", "character")])
  num_vars <- setdiff(all_vars, factor_vars)
  term_names <- attr(terms, "term.labels")
  sapply(term_names, function(term) {
    vars <- unlist(strsplit(term, ":"))
    any(vars %in% num_vars)
  })
}

#' Get a minimum complete model.frame based on the variables in the model
#' @noRd
get_min_data <- function(fixedfr, FUN=mean) {
  do.call(expand.grid, get_var_list(fixedfr, FUN=FUN))
}

#' Extract a named list of variables in the model containing the levels of
#' factors and the mean value of numeric variables
#' @noRd
get_var_list <- function(fixedfr, FUN=mean){
  c(get_fac_list(fixedfr), get_num_list(fixedfr, FUN=FUN))
}

#' Extract a named list of factor levels for each factor in the model
#' @noRd
get_fac_list <- function(fixedfr) {
  terms <- terms(fixedfr)
  res <- .getXlevels(Terms=terms, m=fixedfr)
  if(is.null(res)) list() else res
}

#' Extract named list of mean/FUN values of numeric variables in model
#' @noRd
get_num_list <- function(fixedfr, FUN=mean) { # FUN=function(x) mean(x, na.rm=TRUE)) {
  terms <- terms(fixedfr)
  xvars <- sapply(attr(terms, "variables"), deparse2)[-1L]
  if((yvar <- attr(terms, "response")) > 0)
    xvars <- xvars[-yvar]
  if(!length(xvars)) return(list())
  xlev <- lapply(fixedfr[xvars], function(x) {
    if (is.numeric(x)) FUN(x) else NULL
  })
  res <- xlev[!vapply(xlev, is.null, NA)]
  if(is.null(res)) list() else res
}

deparse2 <- function (x) {
  paste(safeDeparse(x), collapse = " ")
}

safeDeparse <- function(expr,
                        width.cutoff = 500L,
                        backtick = mode(expr) %in%
                          c("call", "expression", "(", "function"),
                        control = c("keepInteger", "showAttributes", "keepNA"),
                        nlines = -1L) {
  deparse(expr = expr, width.cutoff = width.cutoff, backtick = backtick,
          control = control, nlines = nlines)
}
