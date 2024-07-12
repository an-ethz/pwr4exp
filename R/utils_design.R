#' Create data frame for completely randomized design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param replications number of experimental units in each group
#'
#' @return a data.frame with columns for treatment factors and replications
#' @export
#'
#' @examples df.crd(5, 3)
df.crd <- function(treatments, replications){
  if(length(treatments) == 1){
    df = expand.grid(
      treatment = 1:treatments,
      replication = 1:replications,
      stringsAsFactors = TRUE
    )
  } else {
    df = expand.grid(
      facA = 1:treatments[1],
      facB = 1:treatments[2],
      replication = 1:replications,
      stringsAsFactors = TRUE
    )
  }
  df[] <- lapply(df, as.factor)
  return(df)
}

#' Create data frame for randomized complete block design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param blocks number of blocks
#' @return a data.frame with columns for blocks and treatment factors
#' @export
#'
#' @examples df.rcb(5, 2)
df.rcb <- function(treatments, blocks){
  if(length(treatments) == 1){
    df = expand.grid(
      treatment = 1:treatments,
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
  return(df)
}

#' Create data frame for Latin square design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param squares number of replicated squares
#' @param reuse a character string: "row", "col", or "both", indicating reuse of rows or columns or both when replicate a Latin square
#'
#' @return a data.frame with columns for treatment factors, row and column block factors, and squares
#' @export
#'
#' @examples df.lsd(5, 3)
df.lsd <- function(treatments,
                   squares = 1,
                   reuse = c("row", "col", "both")) {
  n <- prod(treatments)
  # lapply(1:n, function(i){(i:(i+n-1)) %% n + 1})
  # do.call(rbind, lapply(1:n, function(i){(i:(i+n-1)) %% n + 1}))
  # unlist(lapply(1:n, function(i){(i:(i+n-1)) %% n + 1}))
  df <- expand.grid(row = 1:n, col = 1:n)
  df$treatment <- unlist(
    lapply(1:n, function(i) {(i:(i+n-1)) %% n + 1})
  )
  if (length(treatments) > 1) {
    df <- cbind(
      df,
      expand.grid(facA = 1:treatments[1], facB = 1:treatments[2])[df$treatment, ]
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
  return(df)
}

#' Create data frame for Crossover design
#'
#' @param treatments a vector representing the number of levels of treatment factors sequentially. For example, \code{treatments=c(2,3)} specifies two treatment factors with 2 and 3 levels, respectively.
#' @param squares number of replicated squares
#'
#' @return a data.frame with columns for treatment factors, individuals (row block factor), period (column block factor), and squares
#' @export
#'
#' @examples df.crossover(4, 5)
df.crossover <- function(treatments, squares){
  df <- df.lsd(treatments, squares, reuse = "col")
  names(df)[2] <- "period"
  df[] <- lapply(df, as.factor)
  return(df)
}


#' Create data frame for split-plot design
#'
#' @param trt.main a vector representing the number of levels of main plot treatment factors sequentially.
#' @param trt.sub a vector representing the number of levels of sub plot treatment factors sequentially.
#' @param replications number of replications at main plots
#'
#' @return a data.frame with columns for main plots, main treatments, and sub-treatments
#' @export
#'
#' @examples df.splitplot(3, 4, 15)
df.splitplot <- function(trt.main, trt.sub, replications){
  df.main <- df.crd(treatments = trt.main, replications)
  df.main$mainplots <- 1:(prod(trt.main)*replications)

  if (length(trt.main) > 1) {
    names(df.main)[1:2] <- c("facA.main", "facB.main")
  } else {names(df.main)[1] <- "trt.main"}


  df.sub <- df.rcb(treatments = trt.sub, blocks = prod(trt.main)*replications)

  if (length(trt.sub) > 1) {
    names(df.sub) <- c("facA.sub", "facB.sub", "mainplots")
  } else {names(df.sub) <- c("trt.sub", "mainplots")}

  df <- merge(df.main, df.sub)
  df <- df[, colnames(df) != "replication"]
  return(df)
}


