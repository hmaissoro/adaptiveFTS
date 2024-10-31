#' Estimate the the standard deviation of the observation error
#'
#' This function estimates the the standard deviation of the observation error using the estimator proposed by \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the standard deviation of the error.
#'
#' @return A data.table with two columns: \code{t} and \code{sig} corresponding to the estimated standard deviation.
#' @export
#'
#' @import data.table
#' @import Rdpack
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [estimate_sigma_cpp()]
#'
#' @examples
#' \dontrun{
#' # Load data
#' data("data_far")
#'
#' # Estimate the standar-deviation of the error term
#' estimate_sigma_cpp(data = data_far, t = c(1/4, 1/2, 3/4))
#'
#' }
#'
#'
estimate_sigma <- function(data, idcol = NULL, tcol = "tobs", ycol = "X", t = c(1/4, 1/2, 3/4)) {
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
#' of the papers \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the empirical autocovariance function.
#' @param lag \code{vector (integer)}. Lag of the autocovariance.
#' @param h \code{numeric (positive vector or scalar)}. The smoothing bandwidth parameter.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
#' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with three columns: \code{t}, \code{lag}, and \code{autocov} corresponding to the estimated autocovariance.
#' @export
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [get_nw_optimal_bw()], [estimate_empirical_autocov_cpp()].
#'
#' @examples
#' \dontrun{
#' # Load data
#' data("data_far")
#'
#' # Estimate empirical autocovariance with a specified bandwidth
#' dt_empirical_autocov <- estimate_empirical_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), lag = c(1, 2), h = 0.1,
#'   kernel_name = "epanechnikov")
#' dt_empirical_autocov
#'
#' # Estimate empirical autocovariance with Cross-Validation bandwidth selection
#' dt_empirical_autocov_cv <- estimate_empirical_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), lag = c(1, 2), h = NULL,
#'   kernel_name = "epanechnikov")
#' dt_empirical_autocov_cv
#' }
#'
estimate_empirical_autocov <- function(data, idcol = NULL, tcol = "tobs", ycol = "X",
                                       t = c(1/4, 1/2, 3/4), lag = c(0, 1, 2), h = NULL,
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

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, kernel_name = kernel_name)
    } else {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, kernel_name = kernel_name)
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Estimation using C++ function
  mat_emp_autocov <- estimate_empirical_autocov_cpp(data = data, t = t, h = h, lag = lag, kernel_name = kernel_name)
  dt_emp_autocov <- data.table::as.data.table(mat_emp_autocov)
  data.table::setnames(x = dt_emp_autocov, new = c("t", "lag", "autocov"))

  return(dt_emp_autocov)
}

#' Estimate empirical \eqn{p}-th order moment of \eqn{X(t)}.
#'
#' This function estimates the \eqn{p}-th order moment of \eqn{X(t)}, used in the empirical study section
#' of the papers \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which the \eqn{p}-th order moment of \eqn{X(t)} is estimated.
#' Each element should be a value between 0 and 1.
#' @param mom_order \code{numeric (positive scalar)}. The order of the moment to be computed (e.g., 1 for mean, 2 for variance).
#' @param h \code{numeric (positive vector or scalar)}. The smoothing bandwidth parameter.
#' Default \code{h = NULL}, in which case the bandwidth will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a scalar, then all curves will be smoothed with the same bandwidth.
#' If \code{h} is a vector, its length must equal the number of curves in \code{data}, with each element corresponding
#' to a curve in the same order as in \code{data}.
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
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [get_nw_optimal_bw()], [estimate_empirical_mom_cpp()].
#'
#' @examples
#' \dontrun{
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
#'   h = NULL,
#'   center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # View the result
#' print(moment_estimates)
#' }
#'
estimate_empirical_mom <- function(data, idcol = NULL, tcol = "tobs", ycol = "X",
                                   t = c(1/4, 1/2, 3/4), mom_order = 1, h = NULL,
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

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, kernel_name = kernel_name)
    } else {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, kernel_name = kernel_name)
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Estimation using C++ function
  mat_mom <- estimate_empirical_mom_cpp(
    data = data, t = t, h = h, mom_order = mom_order,
    center = center, kernel_name = kernel_name)
  dt_mom <- data.table::as.data.table(mat_mom)
  data.table::setnames(x = dt_mom, new = c("t", "mom_order", "mom_estimate"))

  return(dt_mom)
}


#' Estimate Empirical \eqn{X_0(s)X_{\ell}(t)} Autocovariance Function for \eqn{\ell} = 0, 1, ...
#'
#' This function estimates the empirical \eqn{X_0(s)X_{\ell}(t)} autocovariance function for \eqn{\ell} = 0, 1, ...,
#' used in the empirical study of the papers \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{s} in the pair (\code{s}, \code{t}).
#' Must be of the same length as \code{t}.
#' @param t \code{vector (numeric)}. Second argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{t} in the pair (\code{s}, \code{t}).
#' Must be of the same length as \code{s}.
#' @param cross_lag \code{integer (positive integer)}. The lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.
#' @param lag \code{vector (integer)}. Lag for the autocovariance of \eqn{X_0(s)X_{\ell}(t)}.
#' If \code{lag = NULL}, only \eqn{\mathbb{E}X_0(s)X_{\ell}(t)} is returned.
#' @param h \code{numeric (positive vector or scalar)}. Smoothing bandwidth parameter.
#' Defaults to \code{NULL}, in which case \code{h} is estimated via Cross-Validation on a subset of curves.
#' If \code{h} is a scalar, all curves are smoothed with the same bandwidth; if a vector, it should match the number of curves in \code{data}.
#' @param center \code{logical}. If \code{TRUE}, the estimated autocovariance is centered: \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))}. Defaults to \code{FALSE}, providing \eqn{\mathbb{E}X_0(s)X_{\ell}(t)}.
#' @param kernel_name \code{string}. Kernel function for estimation; defaults to "epanechnikov". Supported kernels are: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item{s :}{ First argument in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{t :}{ Second argument in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{cross_lag :}{ Lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{lag :}{ Lags for autocovariance estimation of \eqn{X_0(s)X_{\ell}(t)}; contains \code{NA} if \code{lag = NULL}.}
#'   \item{EXsXt_cross_lag :}{ Mean of \eqn{X_0(s)X_{\ell}(t)}.}
#'   \item{XsXt_autocov :}{ Autocovariance estimates of \eqn{X_0(s)X_{\ell}(t)} for each \code{lag}; contains \code{NA} if \code{lag = NULL}.}
#' }
#'
#' @export
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [get_nw_optimal_bw()], [estimate_empirical_XsXt_autocov_cpp()].
#'
#' #' @examples
#' \dontrun{
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
#'   lag = c(0, 1, 2),
#'   h = 0.1,
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
#'   lag = c(0, 1, 2),
#'   h = 0.1,
#'   center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'dt_empirical_cov_centered
#' }
#'
estimate_empirical_XsXt_autocov <- function(data, idcol = NULL, tcol = "tobs", ycol = "X",
                                            s = c(1/5, 2/5, 4/5),
                                            t = c(1/4, 1/2, 3/4),
                                            cross_lag = 1,
                                            lag = c(0, 1, 2), h = NULL,
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
  if (any(N <= lag))
    stop("'lag' must be lower than the number of curves.")
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

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      sample_curves <- sample(x = 1:N, size = 30)
    } else {
      sample_curves <- 1:N
    }
    if (N > 50) {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, kernel_name = kernel_name)
    } else {
      h <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, kernel_name = kernel_name)
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Estimation using C++ function
  mat_XsXt_autocov <- estimate_empirical_XsXt_autocov_cpp(
    data = data, t = t, s = s, lag = lag, cross_lag = cross_lag,
    h = h, center = center, kernel_name = kernel_name)
  dt_XsXt_autocov <- data.table::as.data.table(mat_XsXt_autocov)
  data.table::setnames(x = dt_XsXt_autocov, new = c("s", "t", "cross_lag", "lag", "EXsXt_cross_lag", "XsXt_autocov"))

  return(dt_XsXt_autocov)
}

