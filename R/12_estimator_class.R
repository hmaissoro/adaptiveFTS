## S3 identity for the adaptive-estimator outputs: a class plus a metadata
## attribute for summary()/plot()/autoplot() to dispatch on.

#' Adaptive functional time series estimator objects
#'
#' The adaptive estimators of \pkg{adaptiveFTS} (for example [estimate_mean()],
#' [estimate_autocov()], [estimate_locreg()], [estimate_cov_segment()] and their
#' risk counterparts) return a [data.table::data.table] carrying an extra S3
#' class so that dedicated `summary()`, `plot()` and
#' [ggplot2::autoplot()][ggplot2::autoplot] methods are available. The object is
#' still a genuine `data.table`: every column and every data.table operation
#' behaves exactly as before.
#'
#' Each result inherits from the shared parent class `adaptiveFTS_est` and from a
#' per-estimator subclass (`mean_est`, `mean_risk`, `autocov_est`,
#' `autocov_risk`, `locreg_est`, `cov_segment_est`, `cov_segment_risk`,
#' `fts_acf`). A short `adaptive_meta` attribute stores context used by the
#' methods (design type, kernel, number of curves, lag, ...).
#'
#' @section Subsetting caveat:
#' A data.table `[` subset returns a plain `data.table` — the leading
#' `adaptiveFTS_est` class is dropped, so `summary()` / `plot()` should be called
#' on the estimator's fresh result. Re-tagging after a subset is possible but not
#' generally needed.
#'
#' @name adaptiveFTS_est
#' @seealso [estimate_mean()], [estimate_autocov()], [estimate_locreg()],
#'   [estimate_cov_segment()].
#' @keywords internal
NULL

#' Tag an estimator result with its adaptive-estimator S3 class
#'
#' @param dt A [data.table::data.table] (modified in place by reference).
#' @param subclass \code{character(1)}. The per-estimator subclass, e.g.
#'   \code{"mean_est"}.
#' @param meta \code{list}. Context stored in the \code{adaptive_meta} attribute
#'   and read back by the \code{summary}/\code{plot} methods.
#' @return \code{dt}, invisibly re-classed and carrying the \code{adaptive_meta}
#'   attribute.
#' @keywords internal
#' @importFrom data.table setattr
.as_adaptive_est <- function(dt, subclass, meta = list()) {
  data.table::setattr(
    x = dt, name = "class",
    value = c(subclass, "adaptiveFTS_est", "data.table", "data.frame"))
  data.table::setattr(x = dt, name = "adaptive_meta", value = meta)
  invisible(dt)
}

#' Retrieve the metadata attached to an adaptive-estimator object
#'
#' @param x An \code{adaptiveFTS_est} object.
#' @param field \code{character(1)} or \code{NULL}. If \code{NULL} (default), the
#'   whole metadata list; otherwise the named element (\code{NULL} when absent).
#' @return The metadata list or the requested element.
#' @keywords internal
.est_meta <- function(x, field = NULL) {
  meta <- attr(x, "adaptive_meta", exact = TRUE)
  if (is.null(meta)) meta <- list()
  if (is.null(field)) return(meta)
  meta[[field]]
}

#' Compact numeric formatter shared by the summary methods
#'
#' @param x \code{numeric}.
#' @param digits \code{integer}. Significant digits. Default \code{4}.
#' @return A formatted \code{character} scalar.
#' @keywords internal
.fmt_num <- function(x, digits = 4L) formatC(x, digits = digits, format = "g")

#' Format a numeric vector as a bracketed min-to-max range string
#'
#' @param x \code{numeric}.
#' @param digits \code{integer}. Significant digits. Default \code{3}.
#' @return A \code{character} scalar giving the min and max, or \code{"NA"}.
#' @keywords internal
.fmt_range <- function(x, digits = 3L) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf("[%s, %s]", .fmt_num(min(x), digits), .fmt_num(max(x), digits))
}
