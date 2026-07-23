#' Convert Raw Curve Observations to the Package Data Format
#'
#' Reshapes raw curve observations into the three-column \code{data.table} that
#' every estimator of the package consumes, and checks that the result meets the
#' assumptions those estimators rely on. All exported estimators call
#' \code{format_data} on their \code{data} argument, so it rarely needs to be
#' called directly.
#'
#' @details
#' Three input layouts are accepted.
#' \itemize{
#'   \item A \code{data.table} (or \code{data.frame}) in long format, with one row
#'     per observation point and at least the columns \code{idcol}, \code{tcol}
#'     and \code{ycol}. The curve index is repeated once per observation point of
#'     that curve.
#'   \item A \code{list} with one element per curve, each element a
#'     \code{data.table} (or \code{data.frame}) holding at least the columns
#'     \code{tcol} and \code{ycol}. Set \code{idcol = NULL}: the curve index is
#'     the position in the list.
#'   \item A \code{list} with one element per curve, each element a \code{list} of
#'     two vectors of equal length named \code{tcol} and \code{ycol}. Set
#'     \code{idcol = NULL}.
#' }
#'
#' Curves are renumbered \eqn{1, \ldots, N} in order of first appearance in
#' \code{data}. That order defines the order of the series: the lag-\eqn{\ell}
#' estimators pair curve \eqn{n} with curve \eqn{n + \ell}, so the curves must
#' arrive in chronological order. Within a curve, rows need neither be contiguous
#' nor sorted; the returned table is always sorted by \code{id_curve}, then by
#' \code{tobs}.
#'
#' The observation points must lie in \eqn{[0, 1]}, the domain the estimators
#' assume for the curves; rescale them beforehand if they are recorded on another
#' scale. Missing values are rejected rather than dropped, so that the number of
#' observation points per curve is the one the caller intends.
#'
#' @param data Raw curve observations, as a \code{data.table} (or
#' \code{data.frame}) in long format, or as a \code{list} with one element per
#' curve. See \code{\link{format_data}} for the accepted layouts and for the
#' \code{id_curve} / \code{tobs} / \code{X} columns they are converted to.
#' @param idcol \code{character(1)} or \code{NULL}. Name of the column holding the
#' curve index when \code{data} is a single table. Must be \code{NULL} when
#' \code{data} is a list of curves.
#' @param tcol \code{character(1)}. Name of the column (or vector) holding the
#' observation points of the curves.
#' @param ycol \code{character(1)}. Name of the column (or vector) holding the
#' values observed at those points.
#'
#' @return A \code{data.table} with three columns, sorted by \code{id_curve} then
#' \code{tobs}:
#' \itemize{
#'   \item \code{id_curve}: the curve index, renumbered \eqn{1, \ldots, N} in
#'     order of first appearance in \code{data}.
#'   \item \code{tobs}: the observation points of each curve.
#'   \item \code{X}: the values observed at each \code{tobs}.
#' }
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' data("data_far")
#'
#' # Long format: one row per observation point.
#' dt <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")
#' head(dt)
#'
#' # One list element per curve: the curve index is the position in the list.
#' curves <- split(dt, by = "id_curve", keep.by = FALSE)
#' head(format_data(data = curves, tcol = "tobs", ycol = "X"))
#'
format_data <- function(data, idcol = NULL, tcol = "tobs", ycol = "X"){
  if (! (is.character(tcol) & length(tcol) == 1) |
      ! (is.character(ycol) & length(ycol) == 1))
    stop("'tcol' and 'ycol' must each be a single column name.", call. = FALSE)

  is_table <- is.data.frame(data)
  is_curve_list <- (! is_table) & is.list(data)
  if (! (is_table | is_curve_list))
    stop("'data' must be a data.table (or data.frame), a list of data.table ",
         "(or data.frame), or a list of list.", call. = FALSE)

  if (is_table) {
    if (is.null(idcol))
      stop("'idcol' must name the curve index column when 'data' is a ",
           "data.table (or data.frame).", call. = FALSE)
    absent <- setdiff(c(idcol, tcol, ycol), names(data))
    if (length(absent))
      stop("'data' has no column named ",
           paste0("'", absent, "'", collapse = ", "), ".", call. = FALSE)
    dt <- data.table::as.data.table(data)[, .SD, .SDcols = c(idcol, tcol, ycol)]
    data.table::setnames(x = dt, new = c("id_curve", "tobs", "X"))
  } else {
    if (! is.null(idcol))
      stop("'idcol' must be NULL when 'data' is a list of curves: the curve ",
           "index is the position in the list.", call. = FALSE)
    if (! length(data))
      stop("'data' is an empty list: there is no curve to format.", call. = FALSE)
    dt <- data.table::rbindlist(lapply(seq_along(data), function(i){
      curve <- data[[i]]
      if (! is.list(curve) | ! all(c(tcol, ycol) %in% names(curve)))
        stop("Element ", i, " of 'data' must be a data.table (or data.frame) ",
             "or a list holding '", tcol, "' and '", ycol, "'.", call. = FALSE)
      if (length(curve[[tcol]]) != length(curve[[ycol]]))
        stop("In element ", i, " of 'data', '", tcol, "' and '", ycol,
             "' have different lengths.", call. = FALSE)
      data.table::data.table("id_curve" = i, "tobs" = curve[[tcol]],
                             "X" = curve[[ycol]])
    }))
  }

  if (! is.numeric(dt[["tobs"]]) | ! is.numeric(dt[["X"]]))
    stop("The observation points ('", tcol, "') and the observed values ('",
         ycol, "') must be numeric.", call. = FALSE)
  if (anyNA(dt))
    stop("'data' holds missing values; remove them before formatting.",
         call. = FALSE)
  if (! all(data.table::between(dt[["tobs"]], 0, 1)))
    stop("The observation points must lie in [0, 1], the domain the estimators ",
         "assume; rescale them before formatting.", call. = FALSE)

  # match() against the unique values numbers the curves in order of first
  # appearance whatever the type of the original index, and without assuming
  # that the rows of a curve are contiguous.
  original_id <- dt[["id_curve"]]
  dt[, id_curve := match(original_id, unique(original_id))]
  data.table::setorder(dt, id_curve, tobs)

  if (anyDuplicated(dt, by = c("id_curve", "tobs")))
    warning("Some curves carry repeated observation points.", call. = FALSE)

  return(dt[])
}
