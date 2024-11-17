#' Curve prediction using the Best Linear Unbiased Predictor (BLUP).
#'
#' This function predict a curve using the adaptive Best Linear Unbiased Predictor proposed by \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t A numeric vector specifying the time points at which to predict the curve \code{id_curve_to_predict}.
#' @param id_curve_to_predict An integer specifying the index of the curve to be predicted. Default is \code{NULL}, which considers the last curve in \code{data}.
#' @param bw_grid A numeric vector of bandwidth grid values for selecting optimal bandwidth parameters for (auto)covariance estimation.
#' Default is \code{NULL}, which sets it in the function.
#' @param use_same_bw A logical value indicating whether the same bandwidth should be used for the arguments \code{s} and \code{t} in the (auto)covariance estimation. Default is \code{FALSE}.
#' @param center A logical value indicating whether the data should be centered before estimating the (auto)covariance. Default is \code{TRUE}.
#' @param correct_diagonal A logical value indicating whether the diagonal of the covariances should be corrected. Default is \code{TRUE}.
#' @param kernel_name A string specifying the kernel to use for estimation. Supported values are \code{"epanechnikov"}, \code{"biweight"},
#'  \code{"triweight"}, \code{"tricube"}, \code{"triangular"}, and \code{"uniform"}. Default is \code{"epanechnikov"}.
#' @return A \code{data.table} containing the predicted curve:
#' \itemize{
#'   \item{\code{t} :}{ The time points at which the curve \code{id_curve_to_predict} is predicted.}
#'   \item{\code{muhat} :}{ The estimates of the mean function.}
#'   \item{\code{prediction} :}{ The adaptive estimates the Best Linear Unbiased Predictor.}
#' }
#' @export
#' @seealso [estimate_locreg()], [estimate_mean()], [estimate_autocov()], [estimate_nw()].
#'
#' @import data.table
#' @importFrom Rdpack reprompt
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#'
predict_curve <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = seq(0.01, 0.99, len = 99),
                          id_curve_to_predict = NULL,
                          bw_grid = NULL,
                          use_same_bw = FALSE,
                          center = TRUE,
                          correct_diagonal = TRUE,
                          kernel_name = "epanechnikov"){

  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(t > 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(use_same_bw, "logical")))
    stop("'use_same_bw' must be TRUE or FALSE.")
  if (! (methods::is(center, "logical")))
    stop("'center' must be TRUE or FALSE.")
  if (! (methods::is(correct_diagonal, "logical")))
    stop("'correct_diagonal' must be TRUE or FALSE.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate the BLUP
  res_blup_one <- estimate_curve_cpp(
    data = data, t = t, id_curve = id_curve_to_predict,
    bw_grid = NULL, use_same_bw = use_same_bw, center = use_same_bw,
    correct_diagonal = use_same_bw,
    kernel_name = "epanechnikov")

  # Return the result
  dt_res <- data.table::as.data.table(res_blup_one$res_blup)
  names(dt_res) <- c("t", "muhat", "prediction")
  return(dt_res)
}
