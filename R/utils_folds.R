#' Create k-fold cross-validation indices
#'
#' Internal replacement for \code{caret::createFolds} so the package no longer
#' depends on \pkg{caret}. The implementation is adapted from
#' \code{caret::createFolds} (Kuhn et al., \pkg{caret}, GPL (>= 2)) and uses only
#' base R, so it consumes the random stream in exactly the same order: for a
#' given random seed it returns fold assignments identical to
#' \code{caret::createFolds(y, k, list = TRUE, returnTrain = FALSE)}.
#'
#' For numeric \code{y} the values are stratified into up to five quantile-based
#' groups before balanced fold assignment, matching caret's behaviour.
#'
#' @param y Vector of outcomes used to balance the folds. In this package
#'   \code{y} is the vector of unique curve identifiers.
#' @param k Integer. Number of folds.
#' @param list Logical. If \code{TRUE} (default) return a named list of index
#'   vectors; otherwise an integer fold-membership vector.
#' @param returnTrain Logical. If \code{TRUE} (and \code{list = TRUE}) return the
#'   training indices for each fold instead of the held-out indices.
#'
#' @return A named list of integer index vectors, or an integer vector of fold
#'   memberships when \code{list = FALSE}.
#'
#' @keywords internal
#' @noRd
.create_folds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {
  if (is.numeric(y)) {
    cuts <- floor(length(y) / k)
    if (cuts < 2) cuts <- 2
    if (cuts > 5) cuts <- 5
    breaks <- unique(stats::quantile(y, probs = seq(0, 1, length.out = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in seq_along(numInClass)) {
      min_reps <- numInClass[i] %/% k
      if (min_reps > 0) {
        spares <- numInClass[i] %% k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) {
          seqVector <- c(seqVector, sample(1:k, spares))
        }
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      } else {
        foldVector[which(y == names(numInClass)[i])] <-
          sample(1:k, size = numInClass[i])
      }
    }
  } else {
    foldVector <- seq(along.with = y)
  }
  if (list) {
    out <- split(seq(along.with = y), foldVector)
    names(out) <- paste("Fold",
                        gsub(" ", "0", format(seq(along.with = out))),
                        sep = "")
    if (returnTrain) {
      out <- lapply(out, function(data, y) y[-data], y = seq(along.with = y))
    }
  } else {
    out <- foldVector
  }
  out
}
