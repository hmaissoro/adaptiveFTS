#' Estimate the Local Regularity Parameters
#'
#' Estimates the local regularity parameters \eqn{H_t} and \eqn{L_t^2} of the
#' underlying process at each point of \code{t}, following Section 3 of
#' Maissoro, Patilea and Vimond (2025). \eqn{H_t} is the local
#' Hölder exponent and \eqn{L_t^2} the squared Hölder constant; both drive the
#' bandwidths of every adaptive estimator of the package.
#'
#' @details
#' Each curve is presmoothed and evaluated at three points spread over a
#' neighbourhood of length \code{Delta} around \code{t}, and the regularity is
#' read off the ratio of the mean squared increments between those points.
#' \code{Delta} drives the bias-variance trade-off of that comparison: too small
#' and the three points carry the same information, too large and the local
#' regularity is averaged away. Left to \code{NULL} it is set from the average
#' number of observation points per curve, \eqn{\widehat\lambda}, as
#' \eqn{\min\{\exp(-(\log\widehat\lambda)^{1/3}),\, 0.2\}}. Near the boundaries of
#' \eqn{[0, 1]} the neighbourhood is shifted inwards rather than truncated.
#'
#' Curves whose presmoothed values fall outside the 2.5\% and 97.5\% quantiles at
#' any of the three points are discarded, so \code{Nused} is smaller than the
#' number of curves. The exponent is clamped to \eqn{[0.1, 1]}: a boundary value
#' usually means the neighbourhood or the presmoothing bandwidth is unsuited to
#' the data rather than a genuinely extreme regularity.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Points of \eqn{[0, 1]} at which the local
#' regularity parameters are estimated.
#' @param Delta \code{numeric (positive)}. Length of the neighbourhood around each
#' point of \code{t} used to estimate the local regularity. Default \code{NULL}
#' sets it from the data; see Details.
#' @param presmooth_bw \code{numeric (positive vector or scalar)}. Bandwidth of
#' the Nadaraya-Watson estimator used to presmooth each curve before the
#' regularity is estimated. A scalar applies the same bandwidth to every curve; a
#' vector must hold one bandwidth per curve, in the order the curves appear in
#' \code{data}. Default \code{NULL} selects a single bandwidth by cross-validation
#' over every curve, as \link{get_nw_optimal_bw} does.
#' @param presmooth_bw_grid \code{vector (numeric)}. Candidate bandwidths of the
#' cross-validation that selects \code{presmooth_bw} when the latter is \code{NULL}.
#' Default \code{NULL} uses the default grid of \link{get_nw_optimal_bw}. Ignored
#' when \code{presmooth_bw} is supplied.
#' @param presmooth_nsubset \code{integer (positive)}. Number of curves used by
#' that cross-validation. Default \code{NULL} uses \code{min(70, floor(N / 2))}
#' curves, where \eqn{N} is the number of curves. Lower it to speed up the
#' selection on large samples. Ignored when \code{presmooth_bw} is supplied.
#' @param kernel_name \code{string}. Kernel of the presmoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#' @param center \code{logical}. If \code{TRUE} (default), the curves are centred
#' before the regularity is estimated.
#'
#' @return A \code{data.table} with one row per point of \code{t} and columns:
#' \itemize{
#'   \item \code{t}: the point at which the local regularity is estimated.
#'   \item \code{locreg_bw}: the presmoothing bandwidth used.
#'   \item \code{Delta}: the length of the neighbourhood used around \code{t}.
#'   \item \code{Nused}: the number of curves that contributed a non-degenerate
#'     estimate at \code{t}.
#'   \item \code{Ht}: the estimated local exponent \eqn{H_t}.
#'   \item \code{Lt2}: the estimated squared Hölder constant \eqn{L_t^2}.
#' }
#'
#' @export
#'
#' @import data.table
#' @importFrom methods is
#'
#' @seealso [get_nw_optimal_bw()], [estimate_mean()], [estimate_autocov()].
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#' @examples
#' data("data_far")
#'
#' dt_locreg <- estimate_locreg(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = seq(0.2, 0.8, length.out = 5), Delta = NULL, presmooth_bw = NULL,
#'   kernel_name = "epanechnikov", center = TRUE)
#' dt_locreg
#'
estimate_locreg <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                            t = 1/2, Delta = NULL, presmooth_bw = NULL,
                            presmooth_bw_grid = NULL, presmooth_nsubset = NULL,
                            kernel_name = "epanechnikov", center = TRUE){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(center, "logical"))
    stop("'center' must be a TRUE or FALSE.")
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

  # Estimate local regularity using C++ function
  mat_reg <- estimate_locreg_cpp(data = data, t = t, Delta = Delta, h = presmooth_bw,
                                    kernel_name = kernel_name, center = center,
                                    presmooth_bw_grid = presmooth_bw_grid,
                                    presmooth_nsubset = presmooth_nsubset)
  dt_reg <- data.table::as.data.table(mat_reg)
  data.table::setnames(x = dt_reg, new = c("t", "locreg_bw", "Delta", "Nused", "Ht", "Lt2"))

  return(.as_adaptive_est(dt_reg, "locreg_est",
                          meta = list(kernel = kernel_name, center = center)))
}
