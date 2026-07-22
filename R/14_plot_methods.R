## ggplot2 visualisations for the adaptive-estimator objects. `autoplot.*`
## builds and returns a ggplot object (the ggplot2 idiom); `plot.*` is a thin
## wrapper that draws it and returns it invisibly. ggplot2 is only a Suggested
## dependency, so every entry point guards with requireNamespace() and refers to
## ggplot2 functions with the `ggplot2::` prefix.

#' @keywords internal
.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required to plot adaptiveFTS objects; ",
         "please install it.", call. = FALSE)
}

#' @keywords internal
.gg_base <- function() {
  ggplot2::theme_minimal(base_size = 11)
}

# Shared colours (kept consistent with the inst/ demos).
.col_main <- "#2c3e50"
.col_accent <- "#c0392b"

#' Plot an adaptive functional time series estimator
#'
#' [ggplot2::autoplot()][ggplot2::autoplot] methods return a ggplot object for
#' the outputs of the adaptive estimators; the matching `plot()` methods draw and
#' invisibly return it. `ggplot2` is a suggested dependency and must be installed.
#' See [adaptiveFTS_est] for the class structure and [summary.mean_est()] for the
#' text summaries.
#'
#' @param object,x An adaptive-estimator object (see [adaptiveFTS_est]).
#' @param which For `locreg_est`, the regularity parameters to display; a subset
#'   of `c("Ht", "Lt")`. Default both.
#' @param ... Passed to the corresponding `autoplot` method (for `plot`) or
#'   unused.
#' @return A [ggplot2::ggplot] object (`plot` methods return it invisibly, after
#'   drawing).
#'
#' @name adaptiveFTS-autoplot
#' @seealso [adaptiveFTS-summary], [estimate_mean()], [estimate_autocov()].
#'
#' @examples
#' \dontrun{
#' data("data_far")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::autoplot(estimate_mean(data = data_far,
#'                                   t = seq(0.1, 0.9, length.out = 20)))
#' }
#' }
NULL

# ---- mean function ---------------------------------------------------------

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.mean_est <- function(object, ...) {
  .require_ggplot2()
  ggplot2::ggplot(object, ggplot2::aes(x = t, y = muhat)) +
    ggplot2::geom_line(linewidth = 0.7, colour = .col_main) +
    ggplot2::geom_point(size = 1.1, colour = .col_main) +
    ggplot2::labs(title = "Adaptive mean function",
                  x = "t", y = expression(hat(mu)(t))) +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.mean_est <- function(x, ...) {
  p <- autoplot.mean_est(x, ...); print(p); invisible(p)
}

# ---- local regularity ------------------------------------------------------

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.locreg_est <- function(object, which = c("Ht", "Lt"), ...) {
  .require_ggplot2()
  which <- match.arg(which, choices = c("Ht", "Lt"), several.ok = TRUE)
  long <- data.table::melt(object, id.vars = "t", measure.vars = which,
                           variable.name = "parameter", value.name = "value")
  ggplot2::ggplot(long, ggplot2::aes(x = t, y = value)) +
    ggplot2::geom_line(linewidth = 0.7, colour = .col_main) +
    ggplot2::geom_point(size = 1.1, colour = .col_main) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(title = "Local regularity parameters", x = "t", y = NULL) +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.locreg_est <- function(x, ...) {
  p <- autoplot.locreg_est(x, ...); print(p); invisible(p)
}

# ---- covariance segment ----------------------------------------------------

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.cov_segment_est <- function(object, ...) {
  .require_ggplot2()
  ggplot2::ggplot(object, ggplot2::aes(x = t, y = cov_segment_hat_corrected)) +
    ggplot2::geom_line(linewidth = 0.7, colour = .col_main) +
    ggplot2::geom_point(size = 1.1, colour = .col_main) +
    ggplot2::labs(title = "Adaptive covariance-segment function",
                  x = "t", y = expression(hat(Gamma)[0](t, t))) +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.cov_segment_est <- function(x, ...) {
  p <- autoplot.cov_segment_est(x, ...); print(p); invisible(p)
}

# ---- (auto)covariance surface ----------------------------------------------

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.autocov_est <- function(object, ...) {
  .require_ggplot2()
  lag <- .est_meta(object, "lag")
  lag_txt <- if (is.null(lag)) "" else sprintf(" (lag = %s)", lag)
  if (all(object$s == object$t)) {
    # Diagonal slice: a curve rather than a surface.
    ggplot2::ggplot(object, ggplot2::aes(x = t, y = autocov)) +
      ggplot2::geom_line(linewidth = 0.7, colour = .col_main) +
      ggplot2::geom_point(size = 1.1, colour = .col_main) +
      ggplot2::labs(title = paste0("Adaptive (auto)covariance diagonal", lag_txt),
                    x = "t", y = "autocovariance") +
      .gg_base()
  } else {
    ggplot2::ggplot(object, ggplot2::aes(x = s, y = t, fill = autocov)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                                    high = .col_accent, midpoint = 0) +
      ggplot2::labs(title = paste0("Adaptive (auto)covariance surface", lag_txt),
                    x = "s", y = "t", fill = "autocov") +
      .gg_base()
  }
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.autocov_est <- function(x, ...) {
  p <- autoplot.autocov_est(x, ...); print(p); invisible(p)
}

# ---- risk curves -----------------------------------------------------------

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.mean_risk <- function(object, ...) {
  .require_ggplot2()
  best <- object[, .SD[which.min(mean_risk)], by = t]
  ggplot2::ggplot(object, ggplot2::aes(x = h, y = mean_risk, colour = factor(t))) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = best, size = 1.8) +
    ggplot2::labs(title = "Adaptive mean risk function",
                  x = "bandwidth h", y = "risk", colour = "t") +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.mean_risk <- function(x, ...) {
  p <- autoplot.mean_risk(x, ...); print(p); invisible(p)
}

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.cov_segment_risk <- function(object, ...) {
  .require_ggplot2()
  best <- object[, .SD[which.min(cov_segment_risk)], by = t]
  ggplot2::ggplot(object, ggplot2::aes(x = h, y = cov_segment_risk, colour = factor(t))) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = best, size = 1.8) +
    ggplot2::labs(title = "Adaptive covariance-segment risk function",
                  x = "bandwidth h", y = "risk", colour = "t") +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.cov_segment_risk <- function(x, ...) {
  p <- autoplot.cov_segment_risk(x, ...); print(p); invisible(p)
}

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.autocov_risk <- function(object, ...) {
  .require_ggplot2()
  d <- data.table::copy(object)
  d[, grp := sprintf("(%.2f, %.2f)", s, t)]
  best <- d[, .SD[which.min(autocov_risk)], by = grp]
  ggplot2::ggplot(d, ggplot2::aes(x = hs, y = autocov_risk, colour = grp)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = best, size = 1.8) +
    ggplot2::labs(title = "Adaptive (auto)covariance risk function",
                  x = "bandwidth (s)", y = "risk", colour = "(s, t)") +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.autocov_risk <- function(x, ...) {
  p <- autoplot.autocov_risk(x, ...); print(p); invisible(p)
}
