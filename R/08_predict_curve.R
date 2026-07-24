#' Curve prediction using the Best Linear Unbiased Predictor (BLUP).
#'
#' @description
#' \strong{Deprecated.} `predict_curve()` is deprecated and will be removed in a
#' future release.
#' It reconstructs a curve from a block system conditioning on the neighbouring
#' curve and the target curve's own partial observations. For the
#' design-weighted, Tikhonov-regularised adaptive BLUP (one-step-ahead
#' prediction of the curve following a conditioning curve), use [blup_fit()]
#' with [predict.blup_fit()], or the one-call wrapper [blup()]. Note that these
#' compute a different quantity, so results are not interchangeable.
#'
#' This function predict a curve using the adaptive Best Linear Unbiased Predictor proposed by Maissoro, Patilea and
#' Vimond (2026).
#'
#' @inheritParams format_data
#' @param t A numeric vector specifying the time points at which to predict the curve \code{id_curve_to_predict}.
#' @param id_curve_to_predict An integer specifying the index of the curve to be predicted. Default is \code{NULL},
#' which considers the last curve in \code{data}.
#' @param bw_grid A numeric vector of bandwidth grid values for selecting optimal bandwidth parameters for
#' (auto)covariance estimation.
#' Default is \code{NULL}, which sets it in the function.
#' @param common_bw A logical value indicating whether a single bandwidth is used
#' for both arguments of the (auto)covariance. Default is \code{FALSE}.
#' @param center_curves A logical value indicating whether the curves are centred
#' before smoothing. Default is \code{TRUE}.
#' @param correct_diagonal A logical value indicating whether the diagonal of the covariances should be corrected.
#' Default is \code{TRUE}.
#' @param kernel_name A string specifying the kernel to use for estimation. Supported values are \code{"epanechnikov"},
#' \code{"biweight"},
#'  \code{"triweight"}, \code{"tricube"}, \code{"triangular"}, and \code{"uniform"}. Default is \code{"epanechnikov"}.
#' @return A \code{data.table} containing the predicted curve:
#' \itemize{
#'   \item{\code{t} :}{ The time points at which the curve \code{id_curve_to_predict} is predicted.}
#'   \item{\code{muhat} :}{ The estimates of the mean function.}
#'   \item{\code{prediction} :}{ The adaptive estimates the Best Linear Unbiased Predictor.}
#' }
#' @export
#' @seealso [blup_fit()], [predict.blup_fit()], [blup()], [estimate_mean()], [estimate_autocov()].
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
#' Functional Time Series. \emph{arXiv preprint} arXiv:2609.xxxxx.
#'
#'
predict_curve <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = seq(0.01, 0.99, len = 99),
                          id_curve_to_predict = NULL,
                          bw_grid = NULL,
                          common_bw = FALSE,
                          center_curves = TRUE,
                          correct_diagonal = TRUE,
                          kernel_name = "epanechnikov"){

  .Deprecated(
    new = "blup_fit",
    msg = paste(
      "'predict_curve()' is deprecated and will be removed in a future release.",
      "Use blup_fit()/predict() or blup() for the design-weighted adaptive BLUP",
      "(note: they compute a different quantity)."))

  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(t >= 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(common_bw, "logical")))
    stop("'common_bw' must be TRUE or FALSE.")
  if (! (methods::is(center_curves, "logical")))
    stop("'center_curves' must be TRUE or FALSE.")
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
    bw_grid = bw_grid, common_bw = common_bw, center = center_curves,
    correct_diagonal = correct_diagonal,
    kernel_name = kernel_name)

  # Return the result
  dt_res <- data.table::as.data.table(res_blup_one$res_blup)
  names(dt_res) <- c("t", "muhat", "prediction")
  return(dt_res)
}
