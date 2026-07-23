## Text summaries for the adaptive-estimator objects, in the style of
## summary.blup_fit (cat-based, returns the object invisibly).

#' Summarise an adaptive functional time series estimator
#'
#' Compact, design-aware text summaries for the objects returned by the adaptive
#' estimators of \pkg{adaptiveFTS}. Each method prints a short report and returns
#' its argument invisibly. See [adaptiveFTS_est] for the class structure and
#' [autoplot.mean_est()] for the companion plots.
#'
#' @param object An adaptive-estimator object (see [adaptiveFTS_est]).
#' @param ... Unused; for S3 compatibility.
#' @return `object`, invisibly.
#'
#' @name adaptiveFTS-summary
#' @seealso [estimate_mean()], [estimate_autocov()], [estimate_locreg()],
#'   [estimate_cov_segment()], [adaptiveFTS_est].
#'
#' @examples
#' \dontrun{
#' data("data_far")
#' summary(estimate_mean(data = data_far, t = seq(0.1, 0.9, length.out = 9)))
#' summary(estimate_autocov(data = data_far,
#'                          s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1))
#' }
NULL

#' @rdname adaptiveFTS-summary
#' @export
summary.locreg_est <- function(object, ...) {
  cat("Adaptive local regularity estimates\n")
  cat(sprintf("  Evaluation points  : %d (t in %s)\n",
              length(unique(object$t)), .fmt_range(object$t)))
  cat(sprintf("  Kernel             : %s (centred: %s)\n",
              .est_meta(object, "kernel"), .est_meta(object, "center")))
  delta <- unique(object$Delta)
  cat(sprintf("  Delta              : %s\n",
              if (length(delta) == 1L) .fmt_num(delta) else .fmt_range(object$Delta)))
  cat(sprintf("  Curves used (Nused): %s\n", .fmt_range(object$Nused)))
  cat(sprintf("  Ht                 : %s\n", .fmt_range(object$Ht)))
  cat(sprintf("  L2t                : %s\n", .fmt_range(object$Lt)))
  invisible(object)
}

#' @rdname adaptiveFTS-summary
#' @export
summary.mean_est <- function(object, ...) {
  cat("Adaptive mean function estimate\n")
  cat(sprintf("  Evaluation points  : %d (t in %s)\n",
              length(unique(object$t)), .fmt_range(object$t)))
  cat(sprintf("  Training curves    : %s\n", .est_meta(object, "N")))
  cat(sprintf("  Kernel             : %s\n", .est_meta(object, "kernel")))
  cat(sprintf("  Optimal bandwidth  : %s\n", .fmt_range(object$optbw)))
  cat(sprintf("  Curves used (PN)   : %s\n", .fmt_range(object$PN)))
  cat(sprintf("  muhat              : %s\n", .fmt_range(object$muhat)))
  invisible(object)
}

#' @rdname adaptiveFTS-summary
#' @export
summary.cov_segment_est <- function(object, ...) {
  cat("Adaptive covariance-segment estimate\n")
  cat(sprintf("  Evaluation points  : %d (t in %s)\n",
              length(unique(object$t)), .fmt_range(object$t)))
  cat(sprintf("  Training curves    : %s\n", .est_meta(object, "N")))
  cat(sprintf("  Kernel             : %s (centred: %s)\n",
              .est_meta(object, "kernel"), .est_meta(object, "center")))
  cat(sprintf("  Optimal bandwidth  : %s\n", .fmt_range(object$optbw)))
  cat(sprintf("  Curves used (PN)   : %s\n", .fmt_range(object$PN)))
  cat(sprintf("  Segment (corrected): %s\n",
              .fmt_range(object$cov_segment_hat_corrected)))
  invisible(object)
}

#' @rdname adaptiveFTS-summary
#' @export
summary.autocov_est <- function(object, ...) {
  lag <- .est_meta(object, "lag")
  title <- if (!is.null(lag) && lag == 0)
    "Adaptive covariance estimate (lag = 0)" else
    sprintf("Adaptive autocovariance estimate (lag = %s)", lag)
  n_pair <- nrow(object)
  n_diag <- sum(object$s == object$t)
  cat(title, "\n", sep = "")
  cat(sprintf("  (s, t) pairs       : %d (%d on the diagonal s = t)\n", n_pair, n_diag))
  cat(sprintf("  Training curves    : %s\n", .est_meta(object, "N")))
  cat(sprintf("  Kernel             : %s (centred: %s)\n",
              .est_meta(object, "kernel"), .est_meta(object, "center")))
  cat(sprintf("  Same bw for s, t   : %s\n", .est_meta(object, "use_same_bw")))
  cat(sprintf("  Bandwidth (s)      : %s\n", .fmt_range(object$optbw_s)))
  cat(sprintf("  Bandwidth (t)      : %s\n", .fmt_range(object$optbw_t)))
  cat(sprintf("  Curves used (PNl)  : %s\n", .fmt_range(object$PNl)))
  cat(sprintf("  autocov            : %s\n", .fmt_range(object$autocov)))
  invisible(object)
}

## Shared body for the risk-function summaries. Lists the risk-minimising
## bandwidth per point for small grids; for larger grids reports the bandwidth
## ranges and the single best point instead of flooding the console.
#' @keywords internal
.summary_risk <- function(object, title, risk_col, by_cols, bw_cols, max_show = 10L) {
  cat(title, "\n", sep = "")
  n_pts <- nrow(unique(object[, by_cols, with = FALSE]))
  n_bw <- .est_meta(object, "n_bw")
  if (is.null(n_bw)) n_bw <- length(unique(object[[bw_cols[1]]]))
  cat(sprintf("  Evaluation points  : %d\n", n_pts))
  cat(sprintf("  Bandwidth grid     : %s values\n", n_bw))
  cat(sprintf("  Training curves    : %s\n", .est_meta(object, "N")))
  cat(sprintf("  Kernel             : %s\n", .est_meta(object, "kernel")))
  best <- object[, {
    i <- which.min(get(risk_col))
    c(lapply(bw_cols, function(b) get(b)[i]), list(risk = get(risk_col)[i]))
  }, by = by_cols]
  data.table::setnames(best, c(by_cols, bw_cols, "risk"))
  fmt_pt <- function(r) paste(sprintf("%s = %s", by_cols,
    vapply(by_cols, function(cc) .fmt_num(best[[cc]][r]), character(1))), collapse = ", ")
  if (nrow(best) <= max_show) {
    cat("  Risk-minimising bandwidth per point:\n")
    for (r in seq_len(nrow(best))) {
      bw <- paste(sprintf("%s* = %s", bw_cols, vapply(bw_cols, function(cc)
        .fmt_num(best[[cc]][r]), character(1))), collapse = ", ")
      cat(sprintf("    (%s): %s (risk = %s)\n", fmt_pt(r), bw, .fmt_num(best$risk[r])))
    }
  } else {
    for (b in bw_cols)
      cat(sprintf("  Optimal %-11s: %s\n", b, .fmt_range(best[[b]])))
    imin <- which.min(best$risk)
    cat(sprintf("  Smallest risk      : %s at (%s)\n", .fmt_num(best$risk[imin]), fmt_pt(imin)))
  }
  invisible(object)
}

#' @rdname adaptiveFTS-summary
#' @export
summary.mean_risk <- function(object, ...) {
  .summary_risk(object, "Adaptive mean risk function",
                risk_col = "mean_risk", by_cols = "t", bw_cols = "h")
}

#' @rdname adaptiveFTS-summary
#' @export
summary.cov_segment_risk <- function(object, ...) {
  .summary_risk(object, "Adaptive covariance-segment risk function",
                risk_col = "cov_segment_risk", by_cols = "t", bw_cols = "h")
}

#' @rdname adaptiveFTS-summary
#' @export
summary.autocov_risk <- function(object, ...) {
  lag <- .est_meta(object, "lag")
  title <- sprintf("Adaptive (auto)covariance risk function (lag = %s)", lag)
  .summary_risk(object, title, risk_col = "autocov_risk",
                by_cols = c("s", "t"), bw_cols = c("hs", "ht"))
}

#' @rdname adaptiveFTS-summary
#' @export
summary.adaptiveFTS_est <- function(object, ...) {
  cls <- setdiff(class(object), c("adaptiveFTS_est", "data.table", "data.frame"))
  cat(sprintf("adaptiveFTS estimator object (%s)\n",
              if (length(cls)) cls[1] else "unclassified"))
  cat(sprintf("  Rows x columns     : %d x %d\n", nrow(object), ncol(object)))
  cat(sprintf("  Columns            : %s\n", paste(names(object), collapse = ", ")))
  invisible(object)
}
