#' Create a data frame of completely randomized design
#'
#' @param treatments An integer-valued vector specifying the treatment structure,
#' in which the length of the vector indicates the number of treatment factors,
#' and each value represents the number of levels for each factor. A maximum of
#' two factors is allowed, and they are arranged in a factorial design.
#' For instance, \code{treatments = n} specifies one treatment factor with n
#' levels, and \code{treatments=c(2,3)} creates a "2x3" factorial design of
#' two treatment factors with 2 and 3 levels, respectively.
#' @param label Optional. A list of character vectors specifying the names of
#' treatment factors and factor levels. Each vector in the list represents a
#' treatment factor, where the name of the vector specifies the name of the
#' factor, and the values in the vector are the labels for that factor's levels.
#' If not provided, factors and levels for one and two treatment factors are
#' labeled as \code{list(trt = c("1", "2", ...))} and
#' \code{list(facA = c("1", "2", ...), facB = c("1", "2", ...))}, respectively.
#' @param replicates The number of experimental units per treatment.
#' @return a data.frame with columns for treatment factors and replicates
df.crd <- function(treatments, label, replicates){
  if(length(treatments) == 1){
    df = expand.grid(
      trt = 1:treatments,
      replication = 1:replicates,
      stringsAsFactors = TRUE
    )
  } else {
    df = expand.grid(
      facA = 1:treatments[1],
      facB = 1:treatments[2],
      replication = 1:replicates,
      stringsAsFactors = TRUE
    )
  }
  df[] <- lapply(df, as.factor)
  if (!missing(label)) {
    names(df)[1:length(treatments)] <- names(label)
    df[[1]] <- factor(df[[1]], labels = label[[1]])
    if(length(treatments) == 2) {
      df[[2]] <- factor(df[[2]], labels = label[[2]])
    }
  }
  return(df)
}

#' Create a data frame of randomized complete block design
#'
#' @inheritParams designCRD
#' @param blocks the number of blocks
#' @return a data.frame with columns for blocks and treatment factors
df.rcbd <- function(treatments, label, blocks){
  if(length(treatments) == 1){
    df = expand.grid(
      trt = 1:treatments,
      block = 1:blocks,
      stringsAsFactors = TRUE
    )
  } else {
    df = expand.grid(
      facA = 1:treatments[1],
      facB = 1:treatments[2],
      block = 1:blocks,
      stringsAsFactors = TRUE
    )
  }
  df[] <- lapply(df, as.factor)
  if (!missing(label)) {
    names(df)[1:length(treatments)] <- names(label)
    df[[1]] <- factor(df[[1]], labels = c(label[[1]]))
    if(length(treatments) == 2) {
      df[[2]] <- factor(df[[2]], labels = c(label[[2]]))
    }
  }
  return(df)
}

#' Create a data frame for Latin square design
#'
#' @inheritParams designCRD
#' @param squares the number of replicated squares
#' @param reuse A character string specifying how to replicate squares when
#' there are multiple squares. Options are: "row" for reusing row blocks, "col"
#' for reusing column blocks, or "both" for reusing both row and column blocks
#' to replicate a single square.
#' @return a data.frame with columns for treatment factors, row and column block factors, and squares
df.lsd <- function(treatments,
                   label,
                   squares = 1,
                   reuse = c("row", "col", "both")) {
  n <- prod(treatments)
  # lapply(1:n, function(i){(i:(i+n-1)) %% n + 1})
  # do.call(rbind, lapply(1:n, function(i){(i:(i+n-1)) %% n + 1}))
  # unlist(lapply(1:n, function(i){(i:(i+n-1)) %% n + 1}))
  df <- expand.grid(row = 1:n, col = 1:n)
  df$trt <- unlist(
    lapply(1:n, function(i) {(i:(i+n-1)) %% n + 1})
  )
  if (length(treatments) > 1) {
    df <- cbind(
      df,
      expand.grid(facA = 1:treatments[1], facB = 1:treatments[2])[df$trt, ]
    )
  }
  if (squares > 1) {
    df <- df[rep(seq_len(prod(treatments)^2), squares), ]
  }
  df[, "square"] <- rep(1:squares, each = prod(treatments)^2)
  reuse = match.arg(reuse)
  if (reuse == "row") {
    df$col <- df$col +
      rep(
        seq(from = 0, by = prod(treatments), length.out = squares),
        each = prod(treatments)^2
      )
  }
  if (reuse == "col") {
    df$row <- df$row +
      rep(
        seq(from = 0, by = prod(treatments), length.out = squares),
        each = prod(treatments)^2
      )
  }
  if (reuse == "both") {
    df$col <- df$col +
      rep(
        seq(from = 0, by = prod(treatments), length.out = squares),
        each = prod(treatments)^2
      )
    df$row <- df$row +
      rep(
        seq(from = 0, by = prod(treatments), length.out = squares),
        each = prod(treatments)^2
      )
  }
  df[] <- lapply(df, as.factor)
  rownames(df) <- 1:nrow(df)
  if (!missing(label)) {
    if (length(treatments) == 1) {
      df[["trt"]] <- factor(df[["trt"]], labels = label[[1]])
      names(df) <- sub("trt", names(label), names(df))
    }
    if (length(treatments) == 2) {
      df[["facA"]] <- factor(df[["facA"]], labels = label[[1]])
      df[["facB"]] <- factor(df[["facB"]], labels = label[[2]])
      names(df) <- sub("facA", names(label)[1], names(df))
      names(df) <- sub("facB", names(label)[2], names(df))
    }
  }
  return(df)
}

#' Create a data frame for Crossover design
#'
#' @inheritParams designCOD
#' @return a data.frame with columns for treatment factors, individuals (row block factor), period (column block factor), and squares
df.cod <- function(treatments, label, squares){
  df <- df.lsd(treatments = treatments, label = label, squares = squares, reuse = "col")
  names(df)[1:2] <- c("subject", "period")
  df[] <- lapply(df, as.factor)
  return(df)
}

#' Create data frame for split-plot design
#'
#' @inheritParams designCRD
#' @param trt.main an integer-valued vector specifying the treatment structure at
#' main plot level, similar to \code{\link{df.crd}}.
#' @param trt.sub an integer-valued vector specifying the treatment structure at
#' sub plot level, similar to `trt.main`.
#' @param replicates the number of experimental units (main plots) per treatment
#' of main plot factors.
#'
#' @return a data.frame with columns for main plots, main treatments, and sub-treatments
df.spd <- function(trt.main, trt.sub, label, replicates){
  df.main <- df.crd(treatments = trt.main, replicates = replicates)
  df.main$mainplot <- 1:(prod(trt.main)*replicates)

  if (length(trt.main) > 1) {
    names(df.main)[1:2] <- c("facA.main", "facB.main")
  } else {names(df.main)[1] <- "trt.main"}
  df.sub <- df.rcbd(treatments = trt.sub, blocks = prod(trt.main)*replicates)
  if (length(trt.sub) > 1) {
    names(df.sub) <- c("facA.sub", "facB.sub", "mainplot")
  } else {names(df.sub) <- c("trt.sub", "mainplot")}
  df <- merge(df.main, df.sub)
  df <- df[, colnames(df) != "replication"]
  if (!missing(label)) {
    if (length(trt.main) == 1) {
      df[["trt.main"]] <- factor(df[["trt.main"]], labels = label[[1]])
      names(df) <- sub("trt.main", names(label)[[1]], names(df))
      if (length(trt.sub) == 1) {
        df[["trt.sub"]] <- factor(df[["trt.sub"]], labels = label[[2]])
        names(df) <- sub("trt.sub", names(label)[[2]], names(df))
      }
      if (length(trt.sub) == 2) {
        df[["facA.sub"]] <- factor(df[["facA.sub"]], labels = label[[2]])
        df[["facB.sub"]] <- factor(df[["facB.sub"]], labels = label[[3]])
        names(df) <- sub("facA.sub", names(label)[[2]], names(df))
        names(df) <- sub("facB.sub", names(label)[[3]], names(df))
      }
    }
    if (length(trt.main) == 2) {
      df[["facA.main"]] <- factor(df[["facA.main"]], labels = label[[1]])
      df[["facB.main"]] <- factor(df[["facB.main"]], labels = label[[2]])
      names(df) <- sub("facA.main", names(label)[[1]], names(df))
      names(df) <- sub("facB.main", names(label)[[2]], names(df))
      if (length(trt.sub) == 1) {
        df[["trt.sub"]] <- factor(df[["trt.sub"]], labels = label[[3]])
        names(df) <- sub("trt.sub", names(label)[[3]], names(df))
      }
      if (length(trt.sub) == 2) {
        df[["facA.sub"]] <- factor(df[["facA.sub"]], labels = label[[3]])
        df[["facB.sub"]] <- factor(df[["facB.sub"]], labels = label[[4]])
        names(df) <- sub("facA.sub", names(label)[[3]], names(df))
        names(df) <- sub("facB.sub", names(label)[[4]], names(df))
      }
    }
  }
  return(df)
}
