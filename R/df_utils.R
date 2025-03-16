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

#' Create a data frame of completely randomized design
#'
#' @param treatments An integer vector where each element represents the number of levels
#' of the corresponding treatment factor. A single integer (e.g., \code{treatments = n})
#' specifies one treatment factor with \code{n} levels. When multiple factors are provided,
#' they are arranged in a factorial treatment factor design. For example,
#' \code{treatments = c(2, 3)} creates a 2x3 factorial design with the first factor having 2 levels
#' and the second factor having 3 levels.
#' @param label Optional. A list of character vectors, each corresponding to a treatment factor.
#' The name of each vector specifies the factor's name, and its elements provide the labels for that factor's levels.
#' If no labels are provided, default labels will be used. For a single treatment factor, the default is
#' \code{list(trt = c("1", "2", ...))}, and for two treatment factors, the default is
#' \code{list(facA = c("1", "2", ...), facB = c("1", "2", ...))}.
#' For split-plot designs, the defaults are similar but include the ".main" and ".sub" suffixes for main plot and subplot factors.
#' For example:
#' \code{list(trt.main = c("1", "2", ...), trt.sub = c("1", "2", ...))}
#' \code{list(facA.main = c("1", "2", ...), facB.main = c("1", "2", ...),
#'       facA.sub = c("1", "2", ...), facB.sub = c("1", "2", ...))}
#' Label sets should be arranged so that the main plot factors come first, followed by the subplot factors.
#' @param replicates The number of experimental units per treatment.
#' @return a data.frame representing the data structure of the design
df.crd <- function(treatments, label, replicates) {
  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  } else {
    if (!is.list(label)) label <- list(label)
    if (length(label) != length(treatments)) {
      stop(sprintf("The number of label sets does not match the number of treatment factors:
                   %d label set(s) provided, %d expected.",
                   length(label), length(treatments)))
    }
    if (!all(sapply(label, length) == treatments)) {
      stop("Some label sets do not have the expected number of labels.")
    }
    if (any(has_name <- names(label) == "")) {
      stop(sprintf("Label set(s) %s also need to be named", paste(which(has_name), collapse = ", ")))
    }
  }

  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  label <- lapply(label, function(x){factor(x, levels = x, labels = x)})

  df <- do.call(expand.grid, c(label, list(replication = seq_len(replicates))))
  return(df)
}

#' Create a data frame of randomized complete block design
#'
#' @inheritParams designCRD
#' @param blocks the number of blocks
#' @return a data.frame representing the data structure of the design
df.rcbd <- function(treatments, label, blocks){
  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  } else {
    if (!is.list(label)) label <- list(label)
    if (length(label) != length(treatments)) {
      stop(sprintf("The number of label sets does not match the number of treatment factors:
                   %d label set(s) provided, %d expected.",
                   length(label), length(treatments)))
    }
    if (!all(sapply(label, length) == treatments)) {
      stop("Some label sets do not have the expected number of labels.")
    }
    if (any(has_name <- names(label) == "")) {
      stop(sprintf("Label set(s) %s also need to be named", paste(which(has_name), collapse = ", ")))
    }
  }

  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  label <- lapply(label, function(x){factor(x, levels = x, labels = x)})

  df <- do.call(expand.grid, c(label, list(block = seq_len(blocks))))

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
#' @return a data.frame representing the data structure of the design
df.lsd <- function(treatments,
                   label,
                   squares = 1,
                   reuse = c("row", "col", "both")) {
  if (missing(label)) {
    label <- lapply(treatments, seq_len)
  } else {
    if (!is.list(label)) label <- list(label)
    if (length(label) != length(treatments)) {
      stop(sprintf("The number of label sets does not match the number of treatment factors:
                   %d label set(s) provided, %d expected.",
                   length(label), length(treatments)))
    }
    if (!all(sapply(label, length) == treatments)) {
      stop("Some label sets do not have the expected number of labels.")
    }
    if (any(has_name <- names(label) == "")) {
      stop(sprintf("Label set(s) %s also need to be named", paste(which(has_name), collapse = ", ")))
    }
  }

  if (is.null(names(label))) {
    if (length(treatments) == 1) {
      names(label) <- "trt"
    } else {
      names(label) <- paste0("fac", LETTERS[seq_along(treatments)])
    }
  }

  label <- lapply(label, function(x){factor(x, levels = x, labels = x)})

  ntrts <- prod(treatments)
  trt_seq <- unlist(lapply(1:ntrts, function(i) {(i:(i+ntrts-1)) %% ntrts + 1}))

  df <- expand.grid(row = 1:ntrts, col = 1:ntrts)
  df <- cbind(df, expand.grid(label)[trt_seq, , drop = FALSE])

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

  return(df)
}

#' Create a data frame for Crossover design
#'
#' @inheritParams designCOD
#' @return a data.frame representing the data structure of the design
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
#' @return a data.frame representing the data structure of the design
df.spd <- function(trt.main, trt.sub, label, replicates){
  if (missing(label)) {
    label.main <- lapply(trt.main, seq_len)
    label.sub <- lapply(trt.sub, seq_len)
  } else {
    if (length(label) != length(trt.main) + length(trt.sub))
      stop(sprintf("The number of label sets does not match the number of treatment factors:
                   \n%d label set(s) provided, %d expected.
                   \nLabel sets should be arranged so that the main plot factors come first,
                   followed by the subplot factors.",
                   length(label), length(trt.main) + length(trt.sub)
                       ))
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

  df.main <- df.crd(treatments = trt.main, label = label.main, replicates = replicates)
  df.main$mainplot <- 1:(prod(trt.main)*replicates)

  df.sub <- df.rcbd(treatments = trt.sub, label = label.sub, blocks = prod(trt.main)*replicates)

  names(df.sub) <- sub("block", "mainplot", names(df.sub))

  df <- merge(df.main, df.sub)
  df <- df[, colnames(df) != "replication"]

  return(df)
}
