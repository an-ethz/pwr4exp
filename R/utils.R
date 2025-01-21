#' @export
#' @import stats
#' @importFrom lme4 lmerControl
designStr <- function(formula,
                      data = NULL,
                      REML = TRUE,
                      subset,
                      weights,
                      na.action,
                      offset,
                      contrasts = NULL,
                      control = lmerControl(),
                      ...) {
  mf <- match.call()
  if (is.null(lme4::findbars(formula))) {
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(model.frame)
    for (i in c("weights", "offset")) {
      if (!eval(bquote(missing(x = .(i)))))
        assign(i, get(i, parent.frame()), environment(formula))
    }
    fr <- eval(mf, parent.frame(1L))
    if (nrow(fr) == 0L)
      stop("0 (non-NA) cases")
    fr <- lme4::factorize(formula, fr, char.only = TRUE)
    attr(fr, "formula") <- formula
    attr(fr, "offset") <- mf$offset
    X <- model.matrix(formula, fr, contrasts)
    X <- chkRank.drop.cols(X, kind = "message+drop.cols", tol = 1e-07)
    X <- checkScaleX(X, kind = "warning")
    desStr <- list(frame = fr, X = X, reTrms = NULL, REML = REML, formula = update(formula, NULL ~ .))
    return(desStr)
  }
  # LMM continues
  mf[[1]] <- quote(lme4::lFormula)
  m <- match(c("formula", "data", "REML", "subset", "weights", "na.action", "offset", "contrast", "control"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  desStr <- eval(mf, parent.frame())
  names(desStr)[1] <- c("frame")
  desStr$reTrms <- desStr$reTrms[c("Zt", "Lambdat", "Lind", "cnms")]
  desStr$reTrms$Lambdat@x <- as.double(desStr$reTrms$Lind)
  desStr$reTrms$Lambdat <- Matrix::forceSymmetric(desStr$reTrms$Lambdat, uplo = "U")
  names(desStr$reTrms) <- c("Zt", "G", "Gind", "cnms")
  n_terms <- lapply(desStr$reTrms$cnms, length)
  n_vp <- lapply(n_terms, function(x){x*(x+1)/2})
  idx <- c(0, cumsum(unlist(n_vp)))[seq_along(n_terms)]
  G_temp <- lapply(seq_along(n_terms), function(i){
    mat_temp <- matrix(0, n_terms[[i]], n_terms[[i]])
    mat_temp[lower.tri(mat_temp, diag = TRUE)] <- seq_len(n_vp[[i]]) + idx[i]
    mat_temp[upper.tri(mat_temp)] <- t(mat_temp)[upper.tri(mat_temp)]
    rownames(mat_temp) <- colnames(mat_temp) <- desStr$reTrms$cnms[[i]]
    return(mat_temp)
  })
  names(G_temp) <- names(desStr$reTrms$cnms)
  desStr$reTrms$G_temp <- G_temp
  desStr$formula <- update(desStr$formula, NULL ~ .)
  return(desStr)
}

#' @export
#' @import stats
corStr <- function(data,
                   sigma2,
                   corcall = NULL) {
  if (is.null(corcall))
    return(list(
      corcall = NULL,
      corframe = data,
      R = sigma2*Matrix::Diagonal(nrow(data))))
  correlation <- try(eval(corcall), silent = TRUE)
  if (inherits(correlation, "try-error")) {
    corcall <- parse(text = paste0("nlme::", deparse(corcall)))[[1]]
    correlation <- try(eval(corcall), silent = TRUE)
    if (inherits(correlation, "try-error")) stop(correlation)
  }

  corform <- corcall$form
  data <- data[, colnames(data) %in% all.vars(corform), drop = FALSE]
  grp.var <- strsplit(deparse1(corform[[2]]), " \\| ")[[1]][2]
  srt.var <- strsplit(deparse1(corform[[2]]), " \\| ")[[1]][1]
  grpform <- reformulate(grp.var)
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

  if (any(class(correlation) %in% c("corExp", "corGaus", "corLin", "corRatio,", "corSpher."))) {
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

  R_list <- nlme::corMatrix(nlme::Initialize(correlation, data = corframe))
  corframe <- corframe[order(corframe$idx_data), ]
  # FIXME: more robust way of matching positions in data and in R_list?
  rownames(corframe) <- corframe$idx_data
  # R <- Matrix::.bdiag(R_list)*sigma2
  R <- Matrix::.bdiag(R_list)
  R <- R[corframe$idx_R, corframe$idx_R]
  R <- sigma2 * R
  rTrms <- list(corcall = corcall, corframe = corframe, R = R)
  return(rTrms)
}

#' @noRd
vcovb_vp <- function(vcomp = NULL, reTrms = NULL) {
  if (is.null(vcomp) || is.null(reTrms)) return(NULL)
  G <- reTrms$G
  G@x <- vcomp[reTrms$Gind]
  return(G)
}

#' @noRd
vcove_vp <- function(rho = NULL, sigma2, rTrms) {
  corcall <- rTrms$corcall
  data <- rTrms$corframe
  if (is.null(corcall) || is.null(rho))
    return(corStr(data, sigma2, corcall)$R)
  corcall$value <- rho
  R <- corStr(data, sigma2, corcall)$R
  return(R)
}

#' @noRd
vcovy_vp <- function(varpar, desStr) {
  reTrms <- desStr$reTrms
  rTrms <- desStr$rTrms
  if (is.null(reTrms)) {
    vcomp <- NULL
    if (is.null(rTrms$corcall))
      rho <- NULL else
        rho <- varpar[seq_along(eval(rTrms$corcall$value))]
  } else {
    vcomp <- varpar[1:max(reTrms$Gind)]
    if (is.null(rTrms$corcall))
      rho <- NULL else
        rho <- varpar[max(reTrms$Gind) + (1:length(eval(rTrms$corcall$value)))]
  }
  sigma2 <- varpar[length(varpar)]
  G <- vcovb_vp(vcomp, reTrms)
  R <- vcove_vp(rho, sigma2, rTrms)
  if (is.null(G) || is.null(reTrms)) # should not happen one is null the other not
    V <- R else V <- Matrix::t(reTrms$Zt) %*% G %*% reTrms$Zt + R
  # return(V)
  return(as.matrix(V)) # FIXME, matrix class
}

#' @noRd
vcovbeta_vp <- function(varpar, desStr) {
  V <- vcovy_vp(varpar, desStr)
  V_inv <- Matrix::solve(V)
  C_inv = t(desStr$X) %*% V_inv %*% desStr$X
  C = Matrix::solve(C_inv)
  return(C)
  # return(as.matrix(C)) # FIXME, matrix class
}

#' @noRd
informat <- function(varpar, desStr) {
  V <- vcovy_vp(varpar = varpar, desStr = desStr)
  V_inv <- Matrix::solve(V)
  P <- V_inv - V_inv %*% desStr$X %*% solve(t(desStr$X) %*% V_inv %*% desStr$X) %*% t(desStr$X) %*% V_inv
  n_params <- length(varpar)
  info_mat <- matrix(0, nrow = n_params, ncol = n_params)
  dV_dvarpar <- numDeriv::jacobian(func = vcovy_vp, x = varpar, desStr = desStr)
  dV_dvarpar <- lapply(1:ncol(dV_dvarpar), function(i)
    array(dV_dvarpar[, i], dim=rep(nrow(V), 2)))
  for (i in 1:n_params) {
    for (j in i:n_params) {
      trace_ij <- sum(diag(P %*% dV_dvarpar[[i]] %*% P %*% dV_dvarpar[[j]]))
      info_mat[i, j] <- 0.5 * trace_ij
      info_mat[j, i] <- info_mat[i, j]
    }
  }
  return(info_mat)
}

#' @noRd
get_contrasts_type3 <- function(object) {
  # Computes contrasts for type III tests with reference to treatment contrast coding
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  X <- object$X # Assumes X is full (column) rank
  # If model has at most one term return Type I contrasts:
  if(ncol(X) <= 1L || length(term_names) <= 1L)
    return(get_contrasts_type1(object))
  # Get 'complete' design matrix:
  rdX <- get_rdX(object, do.warn = TRUE) # treatment contrasts
  # cols for aliased coefs should be removed in X; not in rdX.
  # This makes ginv(X) unique!
  L <- zapsmall(t(MASS::ginv(X) %*% rdX)) # basic contrast matrix
  dimnames(L) <- list(colnames(rdX), colnames(X))
  # Orthogonalize contrasts for terms which are contained in other terms:
  map <- term2colX(terms, X)
  is_contained <- containment(object)
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



#' @noRd
get_contrasts_type2_unfolded <- function(object, which=NULL) {
  # Computes the 'genuine type II contrast' for all terms that are
  # contained in other terms. For all terms which are not contained in other
  # terms, the simple marginal contrast is computed.

  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  X <- object$X

  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(object))
  } else stopifnot(is.character(which), all(which %in% term_names))

  is_contained <- containment(object)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  do_type2 <- setdiff(term_names, do_marginal)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(object, which=do_marginal)
  if(length(do_type2))
    Llist <- c(Llist, get_contrasts_type2(object, which=do_type2))
  Llist[term_names]
}


#' @noRd
get_contrasts_type1 <- function(object) {
  terms <- terms(object$fixed.frame)
  X <- object$X
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

#' @noRd
#' @import stats
get_contrasts_marginal <- function(object, which=NULL) {
  # Computes marginal contrasts.
  #
  # No tests of conformity with coefficients are implemented
  #
  # returns a list
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  X <- object$X
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(object))
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

#' @noRd
get_contrasts_yates <- function(object) {
  # Is this really type 4?
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  X <- object$X
  is_contained <- containment(object)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  not_marginal <- setdiff(term_names, do_marginal)
  # Split not_marginal in do_yates and do_type2:
  do_yates <- need_yates(object)
  do_type2 <- setdiff(not_marginal, do_yates)

  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(object, which=do_marginal)
  if(length(do_yates))
    Llist <- c(Llist, get_yates_contrast(object, which=do_yates))
  if(length(do_type2)) {
    data_classes <- attr(terms(object$frame), "dataClasses")
    Llist <- c(Llist, get_contrasts_type2(object, which=do_type2))
  }
  Llist[term_names]
}

#' @noRd
#' @import stats
get_yates_contrast <- function(object, which=NULL) {
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  X <- object$X
  if(is.null(which)) which <- term_names
  stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  var_list <- get_var_list(object)
  grid <- get_min_data(object)
  form <- update(object$formula, NULL ~ .)
  form <- lme4::nobars(form)
  coef_nm <- colnames(object$X)
  uX <- model.matrix(form, data=grid)
  # Compute LS-means contrast:
  Llist <- lapply(which, function(term) {
    Lt <- model.matrix(formula(paste0("~ 0 + ", term)), data=grid)
    wts <- 1/colSums(Lt) # Yates' weights
    # Lt * c(Lt %*% wts)
    # L <- diag(wts) %*% t(Lt)
    L <- t(sweep(Lt, 2, wts, "*"))
    L %*% uX
  })
  # Check estimability:
  XX <- model.matrix(attr(object$fixed.frame, "terms"), data=object$frame)
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

#' @noRd
#' @import stats
get_contrasts_type2 <- function(object, which=NULL) {
  # Computes the type 2 contrasts - either for all terms or for those
  # included in 'which' (a chr vector naming model terms).
  # returns a list
  X <- object$X
  terms <- attr(object$fixed.frame, "terms") # fixed only
  data_classes <- attr(attr(object$frame, "terms"), "dataClasses") # all vars
  if(is.null(asgn <- attr(X, "assign")))
    stop("design matrix 'X' should have a non-null 'assign' attribute")
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(object))
  } else stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  # Compute containment:
  is_contained <- containment(object)
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


#' @noRd
#' @import stats
containment <- function(object) { # lm or merMod
  # For all terms 'T' in object compute the terms
  # Return a list:
  #   for each term 'T' a vector of terms that contain 'T'.
  terms <- terms(object$fixed.frame)
  data_classes <- attr(terms(object$frame), "dataClasses") # all vars
  # Note: need fixed.only for merMod objects to get dataClasses
  term_names <- attr(terms, "term.labels")
  factor_mat <- attr(terms, "factors")
  lapply(setNames(term_names, term_names), function(term) {
    term_names[term_contain(term, factor_mat, data_classes, term_names)]
  })
}

#' @noRd
#' @import stats
get_rdX <- function(object, do.warn=TRUE) {
  # Compute rank-deficient design-matrix X usign contr.treatment coding.
  #
  # model: terms(model), model.frame(model), fixef(model)
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  df <- object$frame
  # Compute rank-deficient (full) design-matrix, X:
  rdXi <- if(length(term_names)) lapply(term_names, function(trm) {
    form <- as.formula(paste0("~ 0 + ", trm))
    model.matrix(form, data=df) # no contrast arg
  }) else list(model.matrix(~ 1, data=df)[, -1, drop=FALSE])
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

#' @noRd
need_yates <- function(object) {
  ## Do not need yates for:
  ## - continuous variables
  ## - factors that are not contained in other factors
  ## Need yates for all other terms, i.e. terms which are:
  ##  - contained in other terms, AND
  ##  - which are not numeric/continuous
  term_names <- attr(terms(object$fixed.frame), "term.labels")
  cont <- containment(object)
  is_contained <- names(cont[sapply(cont, function(x) length(x) > 0)])
  nmt <- numeric_terms(object)
  num_terms <- names(nmt[nmt])
  term_names[!term_names %in% num_terms &
               term_names %in% is_contained]
}

#' @noRd
#' @import stats
numeric_terms <- function(object) {
  ## Determines for all terms (not just all variables) if the 'dataClass'
  ## is numeric
  ## (interactions involving one or more numerics variables are numeric).
  Terms <- delete.response(terms(object$fixed.frame))
  # all_vars <- all.vars(attr(Terms, "variables"))
  all_vars <- sapply(attr(Terms, "variables")[-1], deparse)
  # FIXME: dataClasses in frame terms are different from fixed.frame terms
  # data_classes <- attr(attr(object$frame, "terms"), "dataClasses")
  data_classes <- attr(terms(object$fixed.frame), "dataClasses")
  var_class <- data_classes[names(data_classes) %in% all_vars]
  # FIXME: should coerce character to factor?
  factor_vars <- names(var_class[var_class %in% c("factor", "ordered")])
  num_vars <- setdiff(all_vars, factor_vars)
  term_names <- attr(attr(object$fixed.frame, "terms"), "term.labels")
  # term_names <- setNames(as.list(term_names), term_names)
  sapply(term_names, function(term) {
    vars <- unlist(strsplit(term, ":"))
    any(vars %in% num_vars)
  })
}

#' @noRd
get_min_data <- function(object, FUN=mean)
  # Get a minimum complete model.frame based on the variables in the model
  do.call(expand.grid, get_var_list(object, FUN=FUN))

#' @noRd
get_var_list <- function(object, FUN=mean)
  # Extract a named list of variables in the model containing the levels of
  # factors and the mean value of numeric variables
  c(get_fac_list(object), get_num_list(object, FUN=FUN))

#' @noRd
get_fac_list <- function(object) {
  # Extract a named list of factor levels for each factor in the model
  res <- .getXlevels(Terms=attr(object$fixed.frame, "terms"), m=object$frame)
  if(is.null(res)) list() else res
}

#' @noRd
get_num_list <- function(object, FUN=mean) { # FUN=function(x) mean(x, na.rm=TRUE)) {
  # Extract named list of mean/FUN values of numeric variables in model
  Terms <- attr(object$fixed.frame, "terms")
  mf <- object$frame
  xvars <- sapply(attr(Terms, "variables"), deparse2)[-1L]
  if((yvar <- attr(Terms, "response")) > 0)
    xvars <- xvars[-yvar]
  if(!length(xvars)) return(list())
  xlev <- lapply(mf[xvars], function(x) {
    if (is.numeric(x)) FUN(x) else NULL
  })
  res <- xlev[!vapply(xlev, is.null, NA)]
  if(is.null(res)) list() else res
}

#' @noRd
#' @import stats
lsmeans_contrasts <- function(object, which=NULL, by = NULL) {
  if(is.null(which) | length(which) > 1)
    stop("exaclty one factor must be specified")
  # right-hand-side formula, fixed-effect-only
  formula <- update(object$formula, NULL ~ .)
  formula <- lme4::nobars(formula)

  # terms types
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  factor_terms <- term_names[!numeric_terms(object)]
  numeric_terms <- term_names[numeric_terms(object)]

  # Contrasts
  Contr <- attr(object$X, "contrasts")
  var_names <- names(get_var_list(object))
  factor_mat <- attr(terms, "factors")
  vars_in_term <- var_names[factor_mat[var_names, which] == 1]

  # Get minimal 'unique rows' data:
  grid <- get_min_data(object)
  Contrasts <- .getXlevels(terms, grid)
  Contrasts[] <- "contr.treatment"

  if (!is.null(by)) { # conditional case
    gridlist <- split(grid, grid[[by]])
    names(gridlist) <- paste0(by, levels(grid[[by]]))
    Llist <- lapply(gridlist, function(grid){
      if (which %in% factor_terms) {
        uX <- model.matrix(formula, data=grid, contrasts.arg=Contr)
        Lt <- model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                           contrasts.arg=Contrasts[vars_in_term])
      }
      if (which %in% numeric_terms) {
        numvar_names <- names(get_num_list(object))
        numvars_in_term <- numvar_names[factor_mat[numvar_names, which] == 1]
        facvars_in_term <- setdiff(vars_in_term, numvar_names)
        grid[[numvars_in_term]] <- 1
        grid0 <- grid
        grid0[[numvars_in_term]] <- 0
        uX <- model.matrix(formula, data=grid, contrasts.arg=Contr[facvars_in_term])
        uX0 <- model.matrix(formula, data=grid0, contrasts.arg=Contr[facvars_in_term])
        uX <- uX - uX0
        Lt <- model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                           contrasts.arg=Contrasts[facvars_in_term])
      }
      wts <- 1/colSums(Lt)
      L <- t(sweep(Lt, 2, wts, "*"))
      return(L %*% uX)
    })
    return(Llist)
  } # marginal case continues
  if (which %in% factor_terms) {
    uX <- model.matrix(formula, data=grid, contrasts.arg=Contr)
    Lt <- model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                       contrasts.arg=Contrasts[vars_in_term])
  }
  if (which %in% numeric_terms) {
    numvar_names <- names(get_num_list(object))
    numvars_in_term <- numvar_names[factor_mat[numvar_names, which] == 1]
    facvars_in_term <- setdiff(vars_in_term, numvar_names)
    grid[, numvars_in_term] <- 1
    grid0 <- grid
    grid0[, numvars_in_term] <- 0
    uX <- model.matrix(formula, data=grid, contrasts.arg=Contr[facvars_in_term])
    uX0 <- model.matrix(formula, data=grid0, contrasts.arg=Contr[facvars_in_term])
    uX <- uX - uX0
    Lt <- model.matrix(formula(paste0("~ 0 + ", which)), data=grid,
                       contrasts.arg=Contrasts[facvars_in_term])
  }
  wts <- 1/colSums(Lt)
  L <- t(sweep(Lt, 2, wts, "*"))
  return(L %*% uX)
}


#' @noRd
#' @import stats
means2beta_contrasts <- function(object){
  terms <- terms(object$fixed.frame)
  term_names <- attr(terms, "term.labels")
  term_names <- setNames(term_names, term_names)
  factor_terms <- term_names[!numeric_terms(object)]
  numeric_terms <- term_names[numeric_terms(object)]
  factor_mat <- attr(terms, "factors")
  var_names <- names(get_var_list(object))
  numvar_names <- names(get_num_list(object))
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
    lsmeans_contrasts(object, wic)
  })
  L <- do.call(rbind, Llist)
  if(colSums(L)[1] == 0) L <- rbind(`(Intercept)` = c(1, numeric(ncol(L)-1)), L)
  return(L)
}

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

#' @noRd
#' @import stats
checkScaleX <- function(X, kind = "warning", tol = 1000) {
  kinds <- eval(formals(lme4::lmerControl)[["check.scaleX"]])
  if (!kind %in% kinds)
    stop(gettextf("unknown check-scale option: %s", kind))
  if (is.null(kind) || kind == "ignore")
    return(X)
  cont.cols <- apply(X, 2, function(z) !all(z %in% c(0, 1)))
  col.sd <- apply(X[, cont.cols, drop = FALSE], 2L, sd)
  sdcomp <- outer(col.sd, col.sd, "/")
  logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
  logsd <- abs(log(col.sd))
  if (any(c(logcomp, logsd) > log(tol))) {
    wmsg <- "Some predictor variables are on very different scales:"
    if (kind %in% c("warning", "stop")) {
      wmsg <- paste(wmsg, "consider rescaling")
      switch(kind, warning = warning(wmsg, call. = FALSE),
             stop = stop(wmsg, call. = FALSE))
    }
    else {
      wmsg <- paste(wmsg, "auto-rescaled (results NOT adjusted)")
      X[, cont.cols] <- sweep(X[, cont.cols, drop = FALSE],
                              2, col.sd, "/")
      attr(X, "scaled:scale") <- setNames(col.sd, colnames(X)[cont.cols])
      if (kind == "warn+rescale")
        warning(wmsg, call. = FALSE)
    }
  }
  else wmsg <- character()
  structure(X, msgScaleX = wmsg)
}

#' @noRd
#' @import stats
chkRank.drop.cols <- function(X, kind, tol = 1e-07, method = "qr") {
  stopifnot(is.matrix(X))
  kinds <- eval(formals(lme4::lmerControl)[["check.rankX"]])
  if (!kind %in% kinds)
    stop(gettextf("undefined option for 'kind': %s", kind))
  if (kind == "ignore")
    return(X)
  p <- ncol(X)
  if (kind == "stop.deficient") {
    if ((rX <- Matrix::rankMatrix(X, tol = tol, method = method)) <
        p)
      stop(gettextf(sub("\n +", "\n", "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);\n                   the fixed effects will be jointly unidentifiable"),
                    rX, p), call. = FALSE)
  }
  else {
    qr.X <- qr(X, tol = tol, LAPACK = FALSE)
    rnkX <- qr.X$rank
    if (rnkX == p)
      return(X)
    msg <- sprintf(ngettext(p - rnkX, "fixed-effect model matrix is rank deficient so dropping %d column / coefficient",
                            "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"),
                   p - rnkX)
    if (kind != "silent.drop.cols")
      (if (kind == "warn+drop.cols")
        warning
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

#' @noRd
deparse2 <- function (x) {
  paste(safeDeparse(x), collapse = " ")
}

#' @noRd
safeDeparse <- function(expr,
                        width.cutoff = 500L,
                        backtick = mode(expr) %in%
                          c("call", "expression", "(", "function"),
                        control = c("keepInteger", "showAttributes", "keepNA"),
                        nlines = -1L) {
  deparse(expr = expr, width.cutoff = width.cutoff, backtick = backtick,
          control = control, nlines = nlines)
}

#' @noRd
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

#' @noRd
get_trts <- function(levs) {
  nlevs <- length(levs)
  ans <- t(cbind(-1, diag(nlevs - 1)))
  rownames(ans) <- levs
  colnames(ans) <- paste(levs[-1], levs[1], sep = " - ")
  ans
}

#' @noRd
#' @import stats
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

#' @noRd
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

#' @noRd
#' @import stats
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

#' @noRd
#' @import stats
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
