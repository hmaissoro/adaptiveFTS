#' Estimate the Risk of the Autocovariance Function Estimator
#'
#' Estimates the risk of the adaptive lag-\eqn{\ell} autocovariance function
#' estimator over a grid of candidate bandwidths, for \eqn{\ell = 0, 1, \ldots}
#' (\eqn{\ell = 0} being the covariance function). Minimising it over the grid at
#' each pair (\code{s}, \code{t}) is what \link{estimate_autocov} does to select
#' its bandwidths.
#'
#' @details
#' Two estimators are covered. With \code{common_bw = TRUE} a single bandwidth is
#' selected for both arguments, the one-bandwidth estimator of
#' Maissoro, Patilea and Vimond (2025); the returned \code{hs}
#' and \code{ht} then hold the same value. With \code{common_bw = FALSE}
#' (default) the risk is minimised over pairs \eqn{(h_s, h_t)}, the two-bandwidth
#' estimator of Maissoro, Patilea and Vimond (2026), which adapts
#' to the regularity of the process at \code{s} and at \code{t} separately. The
#' second is the more flexible but explores the square of the grid, so it costs
#' noticeably more.
#'
#' As for the mean, the risk splits into a bias, a variance and a dependence
#' term, each returned separately. The local regularity parameters are estimated
#' internally at \code{s} and at \code{t}.
#'
#' Left to \code{NULL}, \code{bw_grid} is a 20-point geometric grid running from
#' \eqn{4(N\widehat\lambda)^{-0.9}} to \eqn{4(N\widehat\lambda)^{-1/3}}, where
#' \eqn{N} is the number of curves and \eqn{\widehat\lambda} the average number of
#' observation points per curve.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance
#' function: the points \code{s} of the pairs (\code{s}, \code{t}). Must have the
#' same length as \code{t}.
#' @param t \code{vector (numeric)}. Second argument of the autocovariance
#' function: the points \code{t} of the pairs (\code{s}, \code{t}). Must have the
#' same length as \code{s}.
#' @param lag \code{integer (non-negative)}. Lag \eqn{\ell} of the autocovariance;
#' \code{lag = 0} gives the covariance function.
#' @param bw_grid \code{vector (numeric)}. Candidate bandwidths, from which
#' \link{estimate_autocov} picks the risk-minimising one for each pair. Default
#' \code{NULL} builds the grid from the data; see Details.
#' @param common_bw \code{logical}. If \code{TRUE}, a single bandwidth is selected
#' for both arguments of the autocovariance; if \code{FALSE} (default), one
#' bandwidth per argument. See Details.
#' @param center_curves \code{logical}. If \code{TRUE} (default), the curves are
#' centred before smoothing. It governs the moment and autocovariance
#' estimators only: the local regularity step always centres the curves.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#' @inheritParams estimate_locreg
#'
#' @return A \code{data.table} with one row per (pair, bandwidth) combination and
#' columns:
#' \itemize{
#'   \item \code{s}, \code{t}: the arguments of the autocovariance function.
#'   \item \code{hs}, \code{ht}: the candidate bandwidths for \code{s} and for
#'     \code{t}. Identical when \code{common_bw = TRUE}.
#'   \item \code{PNl}: the number of curves contributing to the estimate at
#'     (\code{s}, \code{t}), \eqn{P_{N,\ell}(s,t;h_s,h_t)}.
#'   \item \code{locreg_bw}: the bandwidth used to estimate the local regularity.
#'   \item \code{Hs}, \code{Ls2}: the estimated local exponent \eqn{H_s} and
#'     squared Hölder constant \eqn{L_s^2} at \code{s}.
#'   \item \code{Ht}, \code{Lt2}: the same quantities at \code{t}.
#'   \item \code{bias_term}, \code{variance_term}, \code{dependence_term}: the
#'     three components of the risk.
#'   \item \code{autocov_risk}: the estimated risk.
#' }
#' @export
#' @seealso [estimate_autocov()], [estimate_locreg()], [estimate_sigma()].
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#' Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
#' Functional Time Series. \emph{arXiv preprint} arXiv:2609.xxxxx.
#'
#' @examples
#' data("data_far")
#'
#' dt_autocov_risk <- estimate_autocov_risk(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1,
#'   bw_grid = seq(0.04, 0.15, length.out = 5), common_bw = TRUE,
#'   center_curves = TRUE, kernel_name = "epanechnikov")
#'
#' # The risk-minimising bandwidth for each pair.
#' dt_autocov_risk[, list(hs = hs[which.min(autocov_risk)]), by = c("s", "t")]
#'
estimate_autocov_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                  s = c(1/5, 2/5, 4/5),
                                  t = c(1/4, 1/2, 3/4),
                                  lag = 1,
                                  bw_grid = NULL,
                                  common_bw = FALSE,
                                  center_curves = TRUE,
                                  presmooth_bw = NULL, Delta = NULL,
                                  presmooth_bw_grid = NULL, presmooth_nsubset = NULL,
                                  kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  .check_locreg_args(presmooth_bw = presmooth_bw, Delta = Delta,
                     presmooth_bw_grid = presmooth_bw_grid,
                     presmooth_nsubset = presmooth_nsubset)

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

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st) ; gc()

  # Estimate risk funciton using C++ function
  mat_autocov_risk <- estimate_autocov_risk_cpp(
    data = data, s = s, t = t, lag = lag, bw_grid = bw_grid,
    common_bw = common_bw, center = center_curves, kernel_name = kernel_name,
    presmooth_bw = presmooth_bw, Delta = Delta,
    presmooth_bw_grid = presmooth_bw_grid,
    presmooth_nsubset = presmooth_nsubset)
  dt_autocov_risk <- data.table::as.data.table(mat_autocov_risk)
  data.table::setnames(
    x = dt_autocov_risk,
    new = c("s", "t", "hs", "ht", "PNl", "locreg_bw", "Hs", "Ls2", "Ht", "Lt2",
            "bias_term", "variance_term", "dependence_term", "autocov_risk"))
  return(.as_adaptive_est(dt_autocov_risk, "autocov_risk",
                          meta = list(kernel = kernel_name, N = N, lag = lag,
                                      common_bw = common_bw, center = center_curves,
                                      n_bw = length(bw_grid))))
}

#' Estimate the Covariance or Autocovariance Function
#'
#' Estimates the adaptive lag-\eqn{\ell} autocovariance function for
#' \eqn{\ell = 0, 1, \ldots} (\eqn{\ell = 0} being the covariance function), using
#' at each pair of points the bandwidths that minimise the estimated risk.
#'
#' @details
#' Unless \code{bw_s} and \code{bw_t} are supplied, the bandwidths are selected
#' pair by pair by minimising the risk of \link{estimate_autocov_risk} over
#' \code{bw_grid}. Set \code{common_bw = TRUE} for the one-bandwidth estimator of
#' Maissoro, Patilea and Vimond (2025) and \code{FALSE}
#' (default) for the two-bandwidth estimator of
#' Maissoro, Patilea and Vimond (2026).
#'
#' \code{center_curves} chooses how the mean is removed. With \code{TRUE}
#' (default) the curves are centred before smoothing, which estimates
#' \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))} in one pass. With
#' \code{FALSE} the two pieces \eqn{\mathbb{E}X_0(s)X_{\ell}(t)} and
#' \eqn{\mu(s)\mu(t)} are estimated separately, the first with a bandwidth from
#' \link{estimate_autocov_risk} and the second with a bandwidth from
#' \link{estimate_mean_risk}. Both are centred estimates; they differ in which
#' bandwidth is applied to which piece.
#'
#' At \code{lag = 0} the observation noise contributes to the estimate wherever
#' the two smoothing windows overlap, that is on and near the diagonal
#' \eqn{s = t}. \code{correct_diagonal = TRUE} (default) subtracts that
#' contribution: an estimate of \eqn{\sigma(s)\sigma(t)} weighted by the overlap
#' of the two kernel weight vectors, which is largest at \eqn{s = t} and decays
#' as the points move apart. It has no effect for \code{lag > 0}, where the noise
#' of two distinct curves is uncorrelated.
#'
#' @inheritParams estimate_autocov_risk
#' @param bw_s \code{vector (numeric)}. Bandwidth to use for \code{s} at each
#' pair. Default \code{NULL} selects it by minimising the estimated risk.
#' @param bw_t \code{vector (numeric)}. Bandwidth to use for \code{t} at each
#' pair. Default \code{NULL} selects it by minimising the estimated risk.
#' @param correct_diagonal \code{logical}. If \code{TRUE} (default), the
#' observation-noise variance is subtracted from the diagonal when
#' \code{lag = 0}. See Details.
#'
#' @return A \code{data.table} with one row per pair (\code{s}, \code{t}) and
#' columns:
#' \itemize{
#'   \item \code{s}, \code{t}: the arguments of the autocovariance function.
#'   \item \code{optbw_s}, \code{optbw_t}: the bandwidths used for \code{s} and
#'     for \code{t}. Identical when \code{common_bw = TRUE}.
#'   \item \code{Hs}, \code{Ls2}: the estimated local exponent \eqn{H_s} and
#'     squared Hölder constant \eqn{L_s^2} at \code{s}.
#'   \item \code{Ht}, \code{Lt2}: the same quantities at \code{t}.
#'   \item \code{PNs}, \code{muhat_s}: the number of curves used for the mean at
#'     \code{s} and the estimated mean there.
#'   \item \code{PNt}, \code{muhat_t}: the same quantities at \code{t}.
#'   \item \code{PNl}: the number of curves contributing to the estimate at
#'     (\code{s}, \code{t}), \eqn{P_{N,\ell}(s,t;h_s,h_t)}.
#'   \item \code{autocov}: the estimated (auto)covariance.
#' }
#' @export
#' @seealso [estimate_autocov_risk()], [estimate_facf()], [estimate_mean()].
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#' Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
#' Functional Time Series. \emph{arXiv preprint} arXiv:2609.xxxxx.
#'
#' @examples
#' data("data_far")
#' dt_small <- data_far[data_far$id_curve <= 20, ]
#' bwg <- seq(0.04, 0.15, length.out = 5)
#'
#' # Lag-1 autocovariance.
#' dt_autocov <- estimate_autocov(
#'   data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1, bw_grid = bwg,
#'   common_bw = FALSE, center_curves = TRUE, correct_diagonal = FALSE,
#'   kernel_name = "epanechnikov")
#' dt_autocov[, list(s, t, optbw_s, optbw_t, PNl, autocov)]
#'
#' # Covariance, with the noise variance removed from the diagonal.
#' dt_cov <- estimate_autocov(
#'   data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/4, 1/2), t = c(1/4, 1/2), lag = 0, bw_grid = bwg,
#'   common_bw = FALSE, center_curves = TRUE, correct_diagonal = TRUE,
#'   kernel_name = "epanechnikov")
#' dt_cov[, list(s, t, PNl, autocov)]
#'
estimate_autocov <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             s = c(1/5, 2/5, 4/5),
                             t = c(1/4, 1/2, 3/4),
                             lag = 1,
                             bw_s = NULL, bw_t = NULL,
                             bw_grid = NULL,
                             common_bw = FALSE,
                             center_curves = TRUE,
                             correct_diagonal = TRUE,
                             presmooth_bw = NULL, Delta = NULL,
                             presmooth_bw_grid = NULL, presmooth_nsubset = NULL,
                             kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")
  .check_locreg_args(presmooth_bw = presmooth_bw, Delta = Delta,
                     presmooth_bw_grid = presmooth_bw_grid,
                     presmooth_nsubset = presmooth_nsubset)

  # Estimate autocovariance using C++ function
  mat_autocov <- estimate_autocov_cpp(
    data = data, s = s, t = t, lag = lag,
    bw_s = bw_s, bw_t = bw_t, bw_grid = bw_grid,
    common_bw = common_bw, center = center_curves,
    correct_diagonal = correct_diagonal, kernel_name = kernel_name,
    presmooth_bw = presmooth_bw, Delta = Delta,
    presmooth_bw_grid = presmooth_bw_grid,
    presmooth_nsubset = presmooth_nsubset)
  dt_autocov <- data.table::as.data.table(mat_autocov)

  data.table::setnames(
    x = dt_autocov,
    new = c("s", "t", "optbw_s", "optbw_t", "Hs", "Ls2", "Ht", "Lt2",
            "PNs", "muhat_s", "PNt", "muhat_t", "PNl", "autocov"))

  return(.as_adaptive_est(dt_autocov, "autocov_est",
                          meta = list(kernel = kernel_name, N = N, lag = lag,
                                      common_bw = common_bw, center = center_curves)))
}

# Autocovariance function estimator : Rubìn et Paranaretos (2020) ----
# Following the Rubìn and Panaretos Equation (B.7), we define Spq_fun and Qpq_fun

#' Weight Sum \eqn{S_{pq}^{(\ell)}} of the Rubìn-Panaretos Autocovariance Estimator
#'
#' Computes the \eqn{S_{pq}^{(\ell)}} term of Equation (B.7) of
#' Rubìn and Panaretos (2020), the kernel weight sum over all
#' pairs of observation points of two curves \eqn{\ell} apart.
#'
#' @inheritParams format_data
#' @param s \code{numeric (scalar)}. First argument of the autocovariance
#' function.
#' @param t \code{numeric (scalar)}. Second argument of the autocovariance
#' function.
#' @param lag \code{integer (non-negative)}. Lag of the autocovariance.
#' @param p,q \code{numeric (integer)}. Exponents of the centred and scaled
#' observation points in the sum.
#' @param bw \code{numeric (positive scalar)}. Bandwidth of the estimator.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @return A \code{numeric} scalar.
#' @keywords internal
#'
.Spq_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                     bw, kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1)) & length(s) == 1))
    stop("'s' must be a numeric scalar value between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))  & length(t) == 1))
    stop("'t' must be a numeric scalar value between 0 and 1.")
  smooth_ker <- .select_kernel(kernel_name)

  # Control and format data (needed before the lag check below, which uses N)
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")
  if ((any(p < 0)| (length(p) > 1) | any(p - floor(p) > 0)) |
      any(q < 0)| (length(q) > 1) | any(q - floor(q) > 0))
    stop("'p' and 'q' must be positive integers.")
  if (! (methods::is(bw, "numeric") & all(data.table::between(bw, 0, 1))  & length(bw) == 1))
    stop("'bw' must be a numeric scalar value between 0 and 1.")

  # Extract observation points
  Tn <- data[id_curve %in% 1:(N - lag), tobs]
  Tn_plus_lag <- data[id_curve %in% (1 + lag):N, tobs]
  rm(data); gc()

  dt_tobs <- data.table::CJ(Tn, Tn_plus_lag)
  if (lag == 0) {
    ind_vec <- data.table::CJ(1:length(Tn), 1:length(Tn_plus_lag))[, which(V1 != V2)]
    dt_tobs <- dt_tobs[ind_vec]
    rm(ind_vec) ; gc()
  }
  xtk_vec <- dt_tobs[, Tn]
  xthj_vec <- dt_tobs[, Tn_plus_lag]
  rm(dt_tobs, Tn, Tn_plus_lag) ; gc()

  # Calculation of the elements to be summed up
  res <- (((xthj_vec - t) / bw ) ** p) * (((xtk_vec - s) / bw ) ** q) *
    (1 / (bw ** 2)) * smooth_ker((xthj_vec - t) / bw) * smooth_ker((xtk_vec - s) / bw )
  Spq_sum <- sum(res)
  rm(res) ; gc()
  Spq <- Spq_sum / (N - lag)
  return(Spq)
}

#' Weighted Cross-Product \eqn{Q_{pq}^{(\ell)}} of the Rubìn-Panaretos Estimator
#'
#' Computes the \eqn{Q_{pq}^{(\ell)}} term of Equation (B.7) of
#' Rubìn and Panaretos (2020), the counterpart of
#' \link{.Spq_fun} weighting the centred cross-products of the observed values.
#'
#' @inheritParams .Spq_fun
#' @param mean_rp \code{data.table}. Mean function estimated at every observation
#' point of every curve, with columns \code{id_curve}, \code{tobs} and
#' \code{muhat_RP}. Default \code{NULL} estimates it from \code{bw_mean}.
#' @param bw_mean \code{numeric (positive scalar)}. Bandwidth of the mean function
#' estimator, used only when \code{mean_rp} is \code{NULL}.
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @return A \code{numeric} scalar.
#' @keywords internal
#'
.Qpq_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                     bw, mean_rp = NULL, bw_mean = NULL,
                     kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") && all(s > 0 & s <= 1) && length(s) == 1))
    stop("'s' must be a numeric scalar value between 0 and 1.")
  if (! (methods::is(t, "numeric") && all(t > 0 & t <= 1)  && length(t) == 1))
    stop("'t' must be a numeric scalar value between 0 and 1.")

  # Control and format data (needed before the lag check below, which uses N)
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")
  if ((any(p < 0)| (length(p) > 1) | any(p - floor(p) > 0)) |
      any(q < 0)| (length(q) > 1) | any(q - floor(q) > 0))
    stop("'p' and 'q' must be positive integers.")
  if (! (methods::is(bw, "numeric") && all(bw > 0 & bw < 1)  && length(bw) == 1))
    stop("'bw' must be a numeric scalar value between 0 and 1.")
  if (is.null(mean_rp)) {
    if (is.null(bw_mean)) {
      stop("If 'mean_rp' is NULL, then 'bw_mean' can not be NULL.")
    } else if (! (methods::is(bw_mean, "numeric") && all(bw_mean > 0 & bw_mean < 1) && length(bw_mean) == 1)) {
      stop("'bw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else if (! (data.table::is.data.table(mean_rp) & all(c("id_curve", "tobs", "muhat_RP") %in% colnames(mean_rp)))) {
      stop("'mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }
  smooth_ker <- .select_kernel(kernel_name)

  # Estimate mean function is it is NULL
  if (is.null(mean_rp)) {
    mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = mean_rp[, tobs], bw = bw_mean, kernel_name = kernel_name)
    mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
    rm(dt_mean) ; gc()
  } else {
    mean_rp <- mean_rp[order(id_curve)]
  }

  # Extract observation points and observed points
  data <- data[order(id_curve)]
  Tn <- data[id_curve %in% 1:(N - lag), tobs]
  Tn_plus_lag <- data[id_curve %in% (1 + lag):N, tobs]
  Yn <- data[id_curve %in% 1:(N - lag), X]
  Yn_plus_lag <- data[id_curve %in% (1 + lag):N, X]
  rm(data); gc()

  # Extract mean function estimates
  mean_rp <- mean_rp[order(id_curve)]
  muhat_Tn <- mean_rp[id_curve %in% 1:(N - lag), muhat_RP]
  muhat_Tn_plus_lag <- mean_rp[id_curve %in% (1 + lag):N, muhat_RP]
  rm(mean_rp) ; gc()
  # repeat data
  dt_tobs <- data.table::CJ(Tn, Tn_plus_lag)
  dt_Y <- data.table::CJ(Yn, Yn_plus_lag)
  dt_mean <- data.table::CJ(muhat_Tn, muhat_Tn_plus_lag)
  if (lag == 0) {
    ind_vec <- data.table::CJ(1:length(Tn), 1:length(Tn_plus_lag))[, which(V1 != V2)]
    dt_tobs <- dt_tobs[ind_vec]
    dt_Y <- dt_Y[ind_vec]
    dt_mean <- dt_mean[ind_vec]
    rm(ind_vec) ; gc()
  }
  # Extract and clean
  xtk_vec <- dt_tobs[, Tn]
  xthj_vec <- dt_tobs[, Tn_plus_lag]
  Ytk_vec <- dt_Y[, Yn]
  Ythj_vec <- dt_Y[, Yn_plus_lag]
  muhat_tk_vec <- dt_mean[, muhat_Tn]
  muhat_xthj_vec <- dt_mean[, muhat_Tn_plus_lag]
  rm(dt_tobs, dt_Y, dt_mean, Tn, Tn_plus_lag, Yn, Yn_plus_lag, muhat_Tn, muhat_Tn_plus_lag) ; gc()

  # Calculate Q function
  Gth <- (Ythj_vec - muhat_xthj_vec) * (Ytk_vec - muhat_tk_vec)
  res <- Gth * (((xthj_vec - t) / bw ) ** p) * (((xtk_vec - s) / bw ) ** q) *
    (1 / (bw ** 2)) * smooth_ker((xthj_vec - t) / bw) * smooth_ker((xtk_vec - s) / bw )

  Qpq_sum <- sum(res)
  rm(res) ; gc()
  Qpq <- Qpq_sum / (N - lag)
  return(Qpq)
}

#' Estimate the Autocovariance Function by the Rubìn-Panaretos Method
#'
#' Estimates the lag-\eqn{\ell} autocovariance function with the local-linear
#' estimator of Rubìn and Panaretos (2020), which smooths every
#' pair of observation points of curves \eqn{\ell} apart with a single bandwidth.
#' It is provided for comparison with the adaptive estimator of
#' \link{estimate_autocov}.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance
#' function: the points \code{s} of the pairs (\code{s}, \code{t}). Must have the
#' same length as \code{t}.
#' @param t \code{vector (numeric)}. Second argument of the autocovariance
#' function: the points \code{t} of the pairs (\code{s}, \code{t}). Must have the
#' same length as \code{s}.
#' @param lag \code{integer (non-negative)}. Lag \eqn{\ell} of the autocovariance.
#' @param bw \code{numeric (positive scalar)}. Bandwidth of the estimator, common
#' to every pair. See \link{estimate_autocov_bw_rp} to select it by
#' cross-validation.
#' @param mean_rp \code{data.table}. Mean function estimated at every observation
#' point of every curve, with columns \code{id_curve}, \code{tobs} and
#' \code{muhat_RP}. Default \code{NULL} estimates it from \code{bw_mean}.
#' @param bw_mean \code{numeric (positive scalar)}. Bandwidth of the mean function
#' estimator, used only when \code{mean_rp} is \code{NULL}.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @return A \code{data.table} with one row per pair (\code{s}, \code{t}) and
#' columns:
#' \itemize{
#'   \item \code{s}, \code{t}: the arguments of the autocovariance function.
#'   \item \code{lag}: the lag \eqn{\ell}.
#'   \item \code{bw_mean}: the bandwidth used for the mean function.
#'   \item \code{bw}: the bandwidth used for the autocovariance.
#'   \item \code{autocovhat_rp}: the estimated autocovariance.
#' }
#' @export
#' @seealso [estimate_autocov_bw_rp()], [estimate_autocov()].
#'
#' @examples
#' \donttest{
#' data("data_far")
#'
#' dt_autocov_rp <- estimate_autocov_rp(
#'   data = data_far[data_far$id_curve <= 10, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1,
#'   bw = 0.1, bw_mean = 0.1, mean_rp = NULL, kernel_name = "epanechnikov")
#' dt_autocov_rp
#' }
#'
estimate_autocov_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
                                lag = 1, bw, bw_mean = NULL, mean_rp = NULL,
                                kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") && all(s > 0 & s <= 1)))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") && all(t > 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )
  if (! (methods::is(bw, "numeric") && all(bw > 0 & bw < 1) && length(bw) == 1))
    stop("'bw' must be a numeric scalar value between 0 and 1.")
  if (is.null(mean_rp)) {
    if (is.null(bw_mean)) {
      stop("If 'mean_rp' is NULL, then 'bw_mean' can not be NULL.")
    } else if (! (methods::is(bw_mean, "numeric") && all(bw_mean > 0 & bw_mean < 1) && length(bw_mean) == 1)) {
      stop("'bw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else if (! (data.table::is.data.table(mean_rp) && all(c("id_curve", "tobs", "muhat_RP") %in% colnames(mean_rp)))) {
      stop("'mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st) ; gc()

  # Estimate mean function is it is NULL
  if (is.null(mean_rp)) {
    mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = mean_rp[, tobs], bw = bw_mean, kernel_name = kernel_name)
    mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
  } else {
    mean_rp <- mean_rp[order(id_curve)]
  }

  # Calculate S_{pq} and Q_{pq}
  autocov_vec <- mapply(function(si, ti, bw, lag, bw_mean, mean_rp, data, ker){
    # Calculate S_{pq} and A_1^{(\ell)},A_2^{(\ell)}, A_3^{(\ell)}
    S00 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 0, bw = bw, kernel_name = ker)
    S01 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 1, bw = bw, kernel_name = ker)
    S02 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 2, bw = bw, kernel_name = ker)
    S10 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 1, q = 0, bw = bw, kernel_name = ker)
    S11 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 1, q = 1, bw = bw, kernel_name = ker)
    S20 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 2, q = 0, bw = bw, kernel_name = ker)

    # calculate A_1^{(\ell)},A_2^{(\ell)}, A_3^{(\ell)} and B^{(\ell)}
    A1 <- S20 * S02 - (S11 ** 2)
    A2 <- S10 * S02 - S01 * S11
    A3 <- S01 * S20 - S10 * S11
    B <- A1 * S00 - A2 * S10 - A3 * S01

    # Calculate Q_{pq}
    Q00 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 0, q = 0, bw = bw, mean_rp = mean_rp, bw_mean = bw_mean, kernel_name = ker)
    Q10 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 1, q = 0, bw = bw, mean_rp = mean_rp, bw_mean = bw_mean, kernel_name = ker)
    Q01 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 0, q = 1, bw = bw, mean_rp = mean_rp, bw_mean = bw_mean, kernel_name = ker)

    # estimate autocovariance
    R <- (A1 * Q00 - A2 * Q10 - A3 * Q01) / B

    return(R)
  }, si = s, ti = t, MoreArgs = list(bw = bw, lag = lag, data = data, bw_mean = bw_mean,
                                     mean_rp = mean_rp, ker = kernel_name))
  dt_res <- data.table::data.table("s" = s, "t" = t, "lag" = lag, "bw_mean" = bw_mean, "autocovhat_rp" = autocov_vec)
  return(dt_res)
}

#' Select the Bandwidth of the Rubìn-Panaretos Autocovariance Estimator
#'
#' Selects the bandwidth of \link{estimate_autocov_rp} by \eqn{K}-fold
#' cross-validation over the curves, as described in
#' Rubìn and Panaretos (2020). Each fold is scored by the squared
#' error between the empirical cross-products of the held-out curves and the
#' lag-0 autocovariance estimated on the others.
#'
#' @details
#' Every candidate bandwidth requires a Rubìn-Panaretos estimate at every pair of
#' observation points of the held-out curves, so the runtime grows with the fourth
#' power of the number of points per curve. Keep \code{bw_grid} short and the
#' number of curves small.
#'
#' @inheritParams format_data
#' @param n_folds \code{integer (positive)}. Number of cross-validation folds.
#' @param bw_grid \code{vector (numeric)}. Candidate bandwidths.
#' @param mean_rp \code{data.table}. Mean function estimated at every observation
#' point of every curve, with columns \code{id_curve}, \code{tobs} and
#' \code{muhat_RP}. Default \code{NULL} estimates it from \code{bw_mean}.
#' @param bw_mean \code{numeric (positive scalar)}. Bandwidth of the mean function
#' estimator, used only when \code{mean_rp} is \code{NULL}.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#'
#' @return A \code{data.table} with one row per candidate bandwidth and columns:
#' \itemize{
#'   \item \code{bw}: the candidate bandwidth.
#'   \item \code{cv_error}: the cross-validation error at \code{bw}. The
#'     bandwidth minimising it is the one to pass to \link{estimate_autocov_rp}.
#' }
#' @export
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @seealso [estimate_autocov_rp()], [estimate_mean_bw_rp()].
#'
estimate_autocov_bw_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                   n_folds = 10, bw_grid = seq(0.001, 0.15, len = 45),
                                   bw_mean = NULL, mean_rp = NULL,
                                   kernel_name = "epanechnikov"){
  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(n_folds < 0)| (length(n_folds) > 1) | any(n_folds - floor(n_folds) > 0) | any(N <= n_folds))
    stop("'n_folds' must be a positive integer lower than the number of curves.")
  if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
    stop("'bw_grid' must be a vector of positive values between 0 and 1.")
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )
  if (is.null(mean_rp)) {
    if (is.null(bw_mean)) {
      stop("If 'mean_rp' is NULL, then 'bw_mean' can not be NULL.")
    } else if (! (methods::is(bw_mean, "numeric") & all(bw_mean > 0 & bw_mean <= 1)  & length(bw_mean) == 1)) {
      stop("'bw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else {
    if (! (data.table::is.data.table(mean_rp) & all(c("id_curve", "tobs", "muhat_RP") %in% colnames(mean_rp))))
      stop("'mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }

  # Estimate mean function is it is NULL
  if (is.null(mean_rp)) {
    mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = mean_rp[, tobs], bw = bw_mean, kernel_name = kernel_name)
    mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
    rm(dt_mean) ; gc()
  } else {
    mean_rp <- mean_rp[order(id_curve)]
  }

  # Create n_folds folds
  fold <- .create_folds(y = unique(data[, id_curve]), k = n_folds, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(lapply(bw_grid, function(BR0, data, mean_rp, fold, kernel_name){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, mean_rp, BR0, kernel_name){
        # split train - test
        dt_test <- data[id_curve %in% unlist(f)]
        dt_test <- dt_test[order(tobs)]
        dt_train <- data[id_curve %in% setdiff(unlist(fold), unlist(f))]
        dt_train <- dt_train[order(tobs)]

        # Extract Tn
        Tn <- dt_test[order(tobs), tobs]
        Tn_grid <- expand.grid(xti = Tn, xtj = Tn)
        xti <- Tn_grid$xti
        xtj <- Tn_grid$xtj

        # Extract Yn
        Yn <- dt_test[order(tobs), X]
        Yn_grid <- expand.grid(Yti = Yn, Ytj = Yn)
        Yti <- Yn_grid$Yti
        Ytj <- Yn_grid$Ytj

        # Extract mean function
        dt_mean_test <- mean_rp[id_curve %in% unlist(f)]
        muhat <- dt_mean_test[order(tobs), muhat_RP]
        muhat_grid <- expand.grid(muhat_ti = muhat, muhat_tj = muhat)
        muhat_ti <- muhat_grid$muhat_ti
        muhat_tj <- muhat_grid$muhat_tj

        rm(Tn, Yn, Tn_grid, Yn_grid, muhat_grid, muhat, dt_mean_test) ; gc()

        # Estimation of mean on fold\f and test on f
        dt_autocov <- estimate_autocov_rp(
          data = dt_train, idcol = "id_curve", tcol = "tobs", ycol = "X",
          s = xti, t = xtj, lag = 0, bw = BR0, bw_mean = bw_mean,
          mean_rp = mean_rp, kernel_name = kernel_name)

        # Calculate the error
        Sqerror <- ((Yti - muhat_ti) * (Ytj - muhat_tj) - dt_autocov[, autocovhat_rp]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, mean_rp = mean_rp, BR0 = BR0, kernel_name = kernel_name, simplify = TRUE),
      error = function(e){
        message("Error in estimating the autocovariance function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold, na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("bw" = BR0, "cv_error" = cv_err)
    return(dt_res)

  }, data = data, mean_rp = mean_rp, fold = fold, kernel_name = kernel_name))

  return(dt_bw)
}

