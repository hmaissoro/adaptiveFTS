#' Convert Data to a \code{data.table} Format
#'
#' This function converts raw data into a \code{data.table} with three columns: the index of the curve, the observed points, and the observed values for each curve.
#'
#' @param data A \code{data.table} (or \code{data.frame}), a \code{list} of \code{data.table} (or \code{data.frame}), or a \code{list} of \code{list}.
#' \itemize{
#'    \item{If \code{data.table}:}{
#'        It should contain the raw curve observations in at least three columns.
#'        \itemize{
#'          \item{\code{idcol} :}{ The name of the column containing the curve index in the sample.
#'                              Each curve index is repeated according to the number of observation points.}
#'          \item{\code{tcol} :}{ The name of the column with observation points associated with each curve index.}
#'          \item{\code{ycol} :}{ The name of the column with observed values at each observation point for each curve index.}
#'        }
#'    }
#'    \item{If \code{list} of \code{data.table}:}{
#'         In this case, each element in the \code{list} represents the observation data of a curve in the form of a \code{data.table} or \code{data.frame}.
#'         Each \code{data.table} contains at least two columns.
#'         \itemize{
#'          \item{\code{tcol} :}{ The name of the column with observation points for the curve.}
#'          \item{\code{ycol} :}{ The name of the column with observed values for the curve.}
#'        }
#'    }
#'    \item{If \code{list} of \code{list}:}{
#'      In this case, \code{data} is a list where each element is the observation data of a curve, given as a \code{list} of two vectors.
#'      \itemize{
#'          \item{\code{tcol} :}{ The vector containing observation points for the curve.}
#'          \item{\code{ycol} :}{ The vector containing observed values for the curve.}
#'        }
#'    }
#' }
#' @param idcol \code{character}. If \code{data} is given as a \code{data.table} or \code{data.frame}, this is the name of the column that holds the curve index.
#' Each curve index is repeated according to the number of observation points. If \code{data} is a \code{list} of \code{data.table} (or \code{data.frame}) or a \code{list} of \code{list}, set \code{idcol = NULL}.
#' @param tcol \code{character}. The name of the column (or vector) containing the observation points for the curves.
#' @param ycol \code{character}. The name of the column with observed values for the curves.
#'
#' @return A \code{data.table} with three columns:
#'          \itemize{
#'            \item{\code{id_curve} :}{ The index of the curve.}
#'            \item{\code{tobs} :}{ The observation points for each curve \code{id_curve}.}
#'            \item{\code{X} :}{ The observed values of the curve at each \code{tobs} point.}
#'         }
#'
#' @import data.table
#' @importFrom methods is
#'
#' @export
#'
format_data <- function(data, idcol = NULL, tcol = "tobs", ycol = "X"){
  # Check if data is a data.table or data.frame
  is_dt_or_df <- methods::is(data, "data.table") | methods::is(data, "data.frame")

  # Check if data is a list of data.table (or data.frame)
  is_list_dt_or_df <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"data.table") | methods::is(element,"data.frame")
    })))

  # Check if data is a list of list
  is_list_of_list <- methods::is(data, "list") &
    all(unlist(lapply(data, function(element){
      methods::is(element,"list")
    })))

  if (! (is_dt_or_df | is_list_dt_or_df | is_list_of_list))
    stop("'data' must of class data.table (or data.frame) or a list of data.table (or data.table) or a list of list.")

  if (is_dt_or_df) {
    if (is.null(idcol))
      stop("If the class of 'data' is data.table (or data.frame), 'idcol' need to be specifyed.")
    if (! all(c(idcol, tcol, ycol) %in% colnames(data))){
      stop("The specified column name 'idcol' or 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::as.data.table(data)
      data <- data[, .SD, .SDcols = c(idcol, tcol, ycol)]
      names(data) <- c("id_curve", "tobs", "X")
      data <- data[, list(id_curve, tobs, X)]
      Mn <- data[, .N, by = id_curve][, N]
      N <- length(Mn)
      id <- unlist(lapply(1:N, function(n, Mn){
        rep(n, Mn[n])
      }, Mn = Mn))
      data[, id_curve := id]
      rm(Mn, N, id)
    }
  } else if (is_list_dt_or_df) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_colname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% colnames(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_colname) {
      stop("The specified column name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }

  } else if (is_list_of_list) {
    if (! is.null(idcol))
      stop("If 'data' is a list of data.table (or data.table) or a list of list, 'idcol' must be NULL.")
    check_vecname <- all(unlist(lapply(data, function(element, tcol, ycol){
      all(c(tcol, ycol) %in% names(element))
    }, tcol = tcol, ycol = ycol)))
    if (! check_vecname) {
      stop("The specified vector name 'tcol' or 'ycol' is incorrect.")
    } else {
      data <- data.table::rbindlist(lapply(1:length(data), function(i, tcol, ycol){
        data[[i]] <- data.table::as.data.table(data[[i]])
        data.table::setnames(x = data[[i]], old = c(tcol, ycol), new = c("tobs", "X"))
        dt <- data.table::data.table("id_curve" = i, data[[i]][, list(tobs, X)])
      }, tcol = tcol, ycol = ycol))
    }
  } else {
    NA
  }
  return(data)
}
