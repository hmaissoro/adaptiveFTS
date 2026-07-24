#' Estimate the Risk of the Covariance Segment Function
#'
#' Estimates the risk \eqn{R_{\Gamma_0}(t; h)} associated with the covariance segment line estimation
#' proposed by Maissoro, Patilea and Vimond (2026).
#'
#' @inheritParams format_data
#' @param t A numeric vector. Observation points where the mean function of the underlying process is estimated.
#' @param bw_grid A numeric vector. A bandwidth grid from which the best smoothing parameter is selected for each
#' \code{t}.
#' Default is \code{NULL}, in which case it is defined as an exponential grid of \eqn{N \times \lambda}.
#' @param center_curves Logical. If \code{TRUE} (default), the curves are centred
#' before smoothing.
#' @param kernel_name Character string. Specifies the kernel function for estimation; default is \code{"epanechnikov"}.
#' Supported kernels include: \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
#' \code{"triangular"}, and \code{"uniform"}.
#'
#' @return A \link[data.table]{data.table} with columns:
#' \itemize{
#'   \item \code{t}: The points at which the risk function is estimated.
#'   \item \code{h}: The candidate bandwidth.
#'   \item \code{PN}: The number of curves used to estimate the mean at \code{t}, i.e., \eqn{P_N(t;h)}.
#'   \item \code{locreg_bw}: The bandwidth used to estimate the local regularity parameters.
#'   \item \code{Ht}: The estimates of the local exponent \eqn{H_t}.
#'   \item \code{Lt2}: The estimates of the Hölder constant \eqn{L_t^2}.
#'   \item \code{bias_term}: The bias term of the risk function.
#'   \item \code{variance_term}: The variance term of the risk function.
#'   \item \code{dependence_term}: The dependence term of the risk function.
#'   \item \code{cov_segment_risk}: The estimated risk of the covariance segment function.
#' }
#'
#' @details
#' The local regularity parameters are estimated within the function using \code{estimate_locreg_cpp}.
#'
#' The dependence term includes contributions from both a term based on \eqn{\mathbb{D}(t; h_t)} derived from
#' fourth-moment tensors,
#' and an empirical autocovariance term computed using \code{estimate_empirical_XsXt_autocov_cpp}.
#'
#'
#' @seealso \link{estimate_mean}, \link{estimate_locreg}, \link{estimate_sigma},
#'          \link{estimate_nw}, \link{estimate_empirical_autocov}
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
#' Functional Time Series. \emph{arXiv preprint} arXiv:2609.xxxxx.
#'
#' @examples
#' data("data_far")
#'
#' dt_risk <- estimate_cov_segment_risk(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.04, 0.15, length.out = 5),
#'   center_curves = TRUE, kernel_name = "epanechnikov")
#'
#' # The risk-minimising bandwidth at each t.
#' dt_risk[, list(h = h[which.min(cov_segment_risk)]), by = "t"]
#'
#' @export
#'
estimate_cov_segment_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                      t = c(1/4, 1/2, 3/4),
                                      bw_grid = NULL,
                                      center_curves = TRUE,
                                      kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    bw_grid <- .default_bw_grid(data)
  }

  # Estimate risk function using C++  function
  dt_risk <- estimate_cov_segment_risk_cpp(
    data = data, t = t, bw_grid = bw_grid,
    center = center_curves, kernel_name = kernel_name)
  dt_risk <- data.table::as.data.table(dt_risk)
  data.table::setnames(x = dt_risk,
                       new = c("t", "h", "PN", "locreg_bw", "Ht", "Lt2", "bias_term",
                               "variance_term", "dependence_term", "cov_segment_risk"))
  return(.as_adaptive_est(dt_risk, "cov_segment_risk",
                          meta = list(kernel = kernel_name, N = N, center = center_curves,
                                      n_bw = length(bw_grid))))
}


#' Estimate Covariance Segment Function for Functional Data
#'
#' Estimates the covariance segment function \eqn{\Gamma_{N,0}(t,t;h_t,h_t)} for functional data
#' using the Nadaraya–Watson estimator with a specified kernel. This is part of the methodology
#' described in Maissoro, Patilea and Vimond (2026).
#'
#' @inheritParams estimate_cov_segment_risk
#' @param bw A numeric vector. Bandwidth to use at each point of \code{t},
#' recycled if a scalar. Default \code{NULL} selects it by minimising the risk
#' estimated by \link{estimate_cov_segment_risk}.
#'
#' @return A \link[data.table]{data.table} containing the following columns:
#' \itemize{
#'   \item{\code{t} :}{ The observation points at which the covariance segment function is estimated.}
#'   \item{\code{optbw} :}{ The optimal bandwidth used to estimate covariance segment function at each \code{t}.}
#'   \item{\code{Ht} :}{ Local exponent estimates for each \code{t}, corresponding to \eqn{H_t}.}
#'   \item{\code{Lt2} :}{ Estimates of the Hölder constant for each \code{t}, corresponding to \eqn{L_t^2}.}
#'   \item{\code{PN} :}{ The number of selected curves used in the estimation for each \code{t}.}
#'   \item{\code{cov_segment_hat} :}{ Uncorrected covariance segment estimate. }
#'   \item{\code{covseg_correction} :}{ Correction term based on measurement error variance. }
#'   \item{\code{cov_segment_hat_corrected} :}{ Final corrected covariance segment estimate. }
#' }
#'
#'
#' @seealso \link{estimate_cov_segment_risk}, \link{estimate_locreg}, \link{estimate_sigma},
#'          \link{estimate_nw}, \link{estimate_empirical_autocov}
#'
#' @import data.table
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
#' Functional Time Series. \emph{arXiv preprint} arXiv:2609.xxxxx.
#'
#' @examples
#' data("data_far")
#'
#' dt_cov_segment <- estimate_cov_segment(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.04, 0.15, length.out = 5),
#'   center_curves = TRUE, kernel_name = "epanechnikov")
#' dt_cov_segment
#'
#' @export
#'
estimate_cov_segment <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                 t = c(1/4, 1/2, 3/4),
                                 bw = NULL,
                                 bw_grid = NULL,
                                 center_curves = TRUE,
                                 kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    bw_grid <- .default_bw_grid(data)
  }

  # Estimate covariance segment function using C++  function
  dt_res <- estimate_cov_segment_cpp(
    data = data, t = t, bw = bw, bw_grid = bw_grid,
    center = center_curves, kernel_name = kernel_name)
  dt_res <- data.table::as.data.table(dt_res)
  data.table::setnames(x = dt_res, new = c("t", "optbw", "Ht", "Lt2", "PN", "cov_segment_hat",
                                           "covseg_correction", "cov_segment_hat_corrected"))
  dt_res[cov_segment_hat < covseg_correction, cov_segment_hat_corrected := cov_segment_hat]
  return(.as_adaptive_est(dt_res, "cov_segment_est",
                          meta = list(kernel = kernel_name, N = N, center = center_curves)))
}
