#' Estimate the the standard deviation of the observation error
#'
#' This function estimates the the standard deviation of the observation error using the estimator proposed by Maissoro, Patilea and Vimond (2025).
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the standard deviation of the error.
#'
#' @return A data.table with two columns: \code{t} and \code{sig} corresponding to the estimated standard deviation.
#' @export
#'
#' @import data.table
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#'
#' @examples
#' # Load data
#' data("data_far")
#'
#' # Estimate the standar-deviation of the error term
#' estimate_sigma(data = data_far, t = c(1/4, 1/2, 3/4))
#'
#'
#'
#'
estimate_sigma <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = c(1/4, 1/2, 3/4)) {
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  # Estimate $\sigma$ using C++ function
  mat_sig <- estimate_sigma_cpp(data = data, t = t)
  dt_sig <- data.table::as.data.table(mat_sig)
  data.table::setnames(x = dt_sig, new = c("t", "sig"))

  return(dt_sig)
}

#' Estimate Empirical Autocovariance Function
#'
#' This function estimates the empirical autocovariance function used in the empirical study section
#' of the papers Maissoro, Patilea and Vimond (2025) and Maissoro, Patilea and Vimond (2026).
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the empirical autocovariance function.
#' @param lag \code{vector (integer)}. Lag of the autocovariance.
#' @param presmooth_bw \code{numeric (positive vector or scalar)}. Bandwidth used
#' to presmooth each curve before the estimation. A scalar applies the same
#' bandwidth to every curve; a vector must hold one bandwidth per curve, in the
#' order the curves appear in \code{data}. Default \code{NULL} selects it by
#' cross-validation, see \link{get_nw_optimal_bw}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
#' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with three columns: \code{t}, \code{lag}, and \code{autocov} corresponding to the estimated autocovariance.
#' @export
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
#' @seealso [get_nw_optimal_bw()].
#'
#' @examples
#' # Load data
#' data("data_far")
#'
#' # Estimate empirical autocovariance with a specified bandwidth
#' dt_empirical_autocov <- estimate_empirical_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), lag = c(1, 2), presmooth_bw = 0.1,
#'   kernel_name = "epanechnikov")
#' dt_empirical_autocov
#'
#' # Estimate empirical autocovariance with Cross-Validation bandwidth selection
#' dt_empirical_autocov_cv <- estimate_empirical_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), lag = c(1, 2), presmooth_bw = NULL,
#'   kernel_name = "epanechnikov")
#' dt_empirical_autocov_cv
#'
#'
estimate_empirical_autocov <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                       t = c(1/4, 1/2, 3/4), lag = c(0, 1, 2), presmooth_bw = NULL,
                                       kernel_name = "epanechnikov"){
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(N <= lag))
    stop("'lag' must be lower than the number of curves.")
  if (! all(methods::is(t, "numeric") & data.table::between(t, 0, 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  presmooth_bw <- .resolve_presmooth_bw(presmooth_bw, data, N, kernel_name)

  # Estimation using C++ function
  mat_emp_autocov <- estimate_empirical_autocov_cpp(data = data, t = t, h = presmooth_bw, lag = lag, kernel_name = kernel_name)
  dt_emp_autocov <- data.table::as.data.table(mat_emp_autocov)
  data.table::setnames(x = dt_emp_autocov, new = c("t", "lag", "autocov"))

  return(dt_emp_autocov)
}

#' Estimate empirical \eqn{p}-th order moment of \eqn{X(t)}.
#'
#' This function estimates the \eqn{p}-th order moment of \eqn{X(t)}, used in the empirical study section
#' of the papers Maissoro, Patilea and Vimond (2025) and Maissoro, Patilea and Vimond (2026).
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which the \eqn{p}-th order moment of \eqn{X(t)} is estimated.
#' Each element should be a value between 0 and 1.
#' @param mom_order \code{numeric (positive scalar)}. The order of the moment to be computed (e.g., 1 for mean, 2 for variance).
#' @param presmooth_bw \code{numeric (positive vector or scalar)}. Bandwidth used
#' to presmooth each curve before the estimation. A scalar applies the same
#' bandwidth to every curve; a vector must hold one bandwidth per curve, in the
#' order the curves appear in \code{data}. Default \code{NULL} selects it by
#' cross-validation, see \link{get_nw_optimal_bw}.
#' @param center \code{logical}. If \code{TRUE}, then the \eqn{p}-th order moment of the centered \eqn{X(t)} is estimated.
#' Default is \code{TRUE}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
#' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with three columns: \code{t}, \code{mom_order}, and \code{mom_estimate} corresponding to
#' the estimated \eqn{p}-th order moment of \eqn{X(t)} at each time point specified in \code{t}.
#'
#' @export
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
#' @seealso [get_nw_optimal_bw()].
#'
#' @examples
#' # Load example data
#' data("data_far")  # Replace with actual data containing observed curves
#'
#' # Define parameters
#' observation_points <- c(0.25, 0.5, 0.75)  # Points at which to estimate moments
#' moment_order <- 2                         # Example: 2nd order moment (variance)
#' bandwidth <- 0.1                          # Smoothing parameter; can be NULL for CV-estimated
#'
#' # Estimate the 2nd order moment (variance) at specified observation points
#' moment_estimates <- estimate_empirical_mom(
#'   data = data_far,
#'   t = observation_points,
#'   mom_order = moment_order,
#'   presmooth_bw = NULL,
#'   center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # View the result
#' print(moment_estimates)
#'
#'
estimate_empirical_mom <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                   t = c(1/4, 1/2, 3/4), mom_order = 1, presmooth_bw = NULL,
                                   center = TRUE, kernel_name = "epanechnikov") {
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (! all(methods::is(t, "numeric") & data.table::between(t, 0, 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  presmooth_bw <- .resolve_presmooth_bw(presmooth_bw, data, N, kernel_name)

  # Estimation using C++ function
  mat_mom <- estimate_empirical_mom_cpp(
    data = data, t = t, h = presmooth_bw, mom_order = mom_order,
    center = center, kernel_name = kernel_name)
  dt_mom <- data.table::as.data.table(mat_mom)
  data.table::setnames(x = dt_mom, new = c("t", "mom_order", "mom_estimate"))

  return(dt_mom)
}


#' Estimate Empirical \eqn{X_0(s)X_{\ell}(t)} Autocovariance Function for \eqn{\ell} = 0, 1, ...
#'
#' This function estimates the empirical \eqn{X_0(s)X_{\ell}(t)} autocovariance function for \eqn{\ell} = 0, 1, ...,
#' used in the empirical study of the papers Maissoro, Patilea and Vimond (2025) and Maissoro, Patilea and Vimond (2026).
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{s} in the pair (\code{s}, \code{t}).
#' Must be of the same length as \code{t}.
#' @param t \code{vector (numeric)}. Second argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{t} in the pair (\code{s}, \code{t}).
#' Must be of the same length as \code{s}.
#' @param cross_lag \code{integer (positive integer)}. The lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.
#' @param autocov_lag \code{vector (integer)}. Lags at which the autocovariance of
#' the scalar series \eqn{n \mapsto X_n(s)X_{n+\ell}(t)} is estimated, \eqn{\ell}
#' being \code{cross_lag}. If \code{NULL}, only \eqn{\mathbb{E}X_0(s)X_{\ell}(t)}
#' is returned.
#' @param presmooth_bw \code{numeric (positive vector or scalar)}. Bandwidth used
#' to presmooth each curve before the estimation. A scalar applies the same
#' bandwidth to every curve; a vector must hold one bandwidth per curve, in the
#' order the curves appear in \code{data}. Default \code{NULL} selects it by
#' cross-validation, see \link{get_nw_optimal_bw}.
#' @param center \code{logical}. If \code{TRUE}, the estimated autocovariance is centered: \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))}. Defaults to \code{FALSE}, providing \eqn{\mathbb{E}X_0(s)X_{\ell}(t)}.
#' @param kernel_name \code{string}. Kernel function for estimation; defaults to "epanechnikov". Supported kernels are: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item{s :}{ First argument in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{t :}{ Second argument in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{cross_lag :}{ Lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{lag :}{ The lags at which the autocovariance of \eqn{X_0(s)X_{\ell}(t)} is estimated; \code{NA} if \code{autocov_lag = NULL}.}
#'   \item{EXsXt_cross_lag :}{ Mean of \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{XsXt_autocov :}{ Autocovariance estimates of \eqn{X_0(s)X_{\ell}(t)} for each \code{autocov_lag}; \code{NA} if \code{autocov_lag = NULL}.}
#' }
#'
#' @export
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
#' @seealso [get_nw_optimal_bw()].
#'
#' @examples
#' # Load data
#' data("data_far")
#'
#' # Example 1: Estimate autocovariance without centering
#' dt_empirical_cov <- estimate_empirical_XsXt_autocov(
#'   data = data_far,
#'   idcol = "id_curve",
#'   tcol = "tobs",
#'   ycol = "X",
#'   s = c(1/5, 2/5, 4/5),
#'   t = c(1/4, 1/2, 3/4),
#'   cross_lag = 1,
#'   autocov_lag = c(0, 1, 2),
#'   presmooth_bw = 0.1,
#'   center = FALSE,
#'   kernel_name = "epanechnikov"
#' )
#' dt_empirical_cov
#'
#' # Example 2: Estimate autocovariance with centering
#' dt_empirical_cov_centered <- estimate_empirical_XsXt_autocov(
#'   data = data_far,
#'   idcol = "id_curve",
#'   tcol = "tobs",
#'   ycol = "X",
#'   s = c(1/5, 2/5, 4/5),
#'   t = c(1/4, 1/2, 3/4),
#'   cross_lag = 1,
#'   autocov_lag = c(0, 1, 2),
#'   presmooth_bw = 0.1,
#'   center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'dt_empirical_cov_centered
#'
#'
estimate_empirical_XsXt_autocov <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                            s = c(1/5, 2/5, 4/5),
                                            t = c(1/4, 1/2, 3/4),
                                            cross_lag = 1,
                                            autocov_lag = c(0, 1, 2), presmooth_bw = NULL,
                                            center = FALSE,
                                            kernel_name = "epanechnikov"){
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (any(N <= autocov_lag))
    stop("'autocov_lag' must be lower than the number of curves.")
  if (! all(methods::is(t, "numeric") & data.table::between(t, 0, 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (any(cross_lag < 0)| (length(cross_lag) > 1) | any(cross_lag - floor(cross_lag) > 0) | any(N <= cross_lag))
    stop("'cross_lag' must be a positive integer lower than the number of curves.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st) ; gc()

  presmooth_bw <- .resolve_presmooth_bw(presmooth_bw, data, N, kernel_name)

  # Estimation using C++ function
  mat_XsXt_autocov <- estimate_empirical_XsXt_autocov_cpp(
    data = data, t = t, s = s, lag = autocov_lag, cross_lag = cross_lag,
    h = presmooth_bw, center = center, kernel_name = kernel_name)
  dt_XsXt_autocov <- data.table::as.data.table(mat_XsXt_autocov)
  data.table::setnames(x = dt_XsXt_autocov, new = c("s", "t", "cross_lag", "lag", "EXsXt_cross_lag", "XsXt_autocov"))

  return(dt_XsXt_autocov)
}

