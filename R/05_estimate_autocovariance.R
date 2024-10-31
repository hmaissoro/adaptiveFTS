#' Estimate the Risk of the Covariance or Autocovariance Function
#'
#' This function estimates the risk function of the adaptive lag-\eqn{\ell} autocovariance function estimator, where \eqn{\ell} = 0, 1, ...,
#' using one bandwidth parameter as proposed in \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}
#' or two bandwidth parameters as proposed in \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. The first argument of the autocovariance function, corresponding to observation points \code{s} in the pair (\code{s}, \code{t}). Must be of the same length as \code{t}.
#' @param t \code{vector (numeric)}. The second argument of the autocovariance function, corresponding to observation points \code{t} in the pair (\code{s}, \code{t}). Must be of the same length as \code{s}.
#' @param lag \code{integer (positive integer)}. The lag of the autocovariance.
#' @param bw_grid \code{vector (numeric)}. Bandwidth grid for selecting the optimal smoothing parameter for each pair (\code{s}, \code{t}). Defaults to \code{NULL}, which generates an exponential grid of \eqn{N \lambda}.
#' @param use_same_bw \code{logical}. Indicates whether the same bandwidth should be used for both \code{s} and \code{t}. Defaults to \code{FALSE}.
#' @param center \code{logical}. If \code{TRUE}, centers the data before estimation. Default is \code{TRUE}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov". Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} containing the following columns:
#' \itemize{
#'   \item{\code{s} :}{ The first argument of the autocovariance function.}
#'   \item{\code{t} :}{ The second argument of the autocovariance function.}
#'   \item{\code{hs} :}{ Candidate bandwidth for \code{s}. If \code{use_same_bw = TRUE}, \code{hs} and \code{ht} will contain the same values.}
#'   \item{\code{ht} :}{ Candidate bandwidth for \code{t}.}
#'   \item{\code{PNl} :}{ Number of curves used in the autocovariance estimation at (\code{s}, \code{t}). Corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
#'   \item{\code{locreg_bw} :}{ Bandwidth used for estimating local regularity parameters.}
#'   \item{\code{Hs} :}{ Local exponent estimates for \code{s}, denoted as \eqn{H_s}.}
#'   \item{\code{Ls} :}{ Estimates of the Hölder constant for \code{s}, corresponding to \eqn{L_s^2}.}
#'   \item{\code{Ht} :}{ Local exponent estimates for \code{t}, denoted as \eqn{H_t}.}
#'   \item{\code{Lt} :}{ Estimates of the Hölder constant for \code{t}, corresponding to \eqn{L_t^2}.}
#'   \item{\code{bias_term} :}{ Bias term of the risk function.}
#'   \item{\code{variance_term} :}{ Variance term of the risk function.}
#'   \item{\code{dependence_term} :}{ Dependence term of the risk function.}
#'   \item{\code{autocov_risk} :}{ Estimated risk of the covariance/autocovariance function.}
#' }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_XsXt_autocov()].
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Load data
#' data("data_far")
#'
#' # Estimate risk function with same bandwidth for s and t
#' dt_autocov_risk <- estimate_autocov_risk(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
#'   bw_grid = NULL, use_same_bw = TRUE, center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # Visualize mean risk function for different (s, t) pairs
#' dt_dcast <- data.table::dcast(data = dt_autocov_risk, formula = hs ~ s + t, value.var = "autocov_risk")
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.2, 0.25)" = `0.2_0.25`)],
#'                       main = "lag = 1 - (s, t) = (0.2, 0.25)",
#'                       xlab = "h", ylab = "Risk Function"),
#'     dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.4, 0.5)" = `0.4_0.5`)],
#'                       main = "lag = 1 - (s, t) = (0.4, 0.5)",
#'                       xlab = "h", ylab = "Risk Function"),
#'     dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.8, 0.75)" = `0.8_0.75`)],
#'                       main = "lag = 1 - (s, t) = (0.8, 0.75)",
#'                       xlab = "h", ylab = "Risk Function")
#'   ),
#'   nrow = 3
#' )
#'
#' # Estimate risk function with separate bandwidths for s and t
#' dt_autocov_risk_2bw <- estimate_autocov_risk(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
#'   bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # Display selected columns of the results
#' dt_autocov_risk_2bw[, .(s, t, lag, hs, ht, autocov_risk)]
#' }
#'
estimate_autocov_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                  s = c(1/5, 2/5, 4/5),
                                  t = c(1/4, 1/2, 3/4),
                                  lag = 1,
                                  bw_grid = NULL,
                                  use_same_bw = FALSE,
                                  center = TRUE,
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

  if (! is.null(bw_grid)) {
    if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
      stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    K <- 20
    b0 <- 4 * (N * lambdahat) ** (- 0.9)
    bK <- 4 * (N * lambdahat) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a, lambdahat) ; gc()
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
    use_same_bw = use_same_bw, center = center, kernel_name = kernel_name)
  dt_autocov_risk <- data.table::as.data.table(mat_autocov_risk)
  data.table::setnames(
    x = dt_autocov_risk,
    new = c("s", "t", "hs", "ht", "PNl", "locreg_bw", "Hs", "Ls", "Ht", "Lt",
            "bias_term", "variance_term", "dependence_term", "autocov_risk"))
  return(dt_autocov_risk)
}

#' Estimate the Covariance or Autocovariance Function
#'
#' This function estimates the adaptive lag-\eqn{\ell} autocovariance function, where \eqn{\ell} = 0, 1, ...,
#' using one bandwidth parameter as proposed in \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}
#' or two bandwidth parameters as proposed in \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. The first argument of the autocovariance function, corresponding to observation points \code{s} in the pair (\code{s}, \code{t}). Must be of the same length as \code{t}.
#' @param t \code{vector (numeric)}. The second argument of the autocovariance function, corresponding to observation points \code{t} in the pair (\code{s}, \code{t}). Must be of the same length as \code{s}.
#' @param lag \code{integer (positive integer)}. The lag of the autocovariance.
#' @param optbw_s \code{vector (numeric)}. The optimal bandwidths for \code{s}. Default is \code{NULL}.
#' @param optbw_t \code{vector (numeric)}. The optimal bandwidths for \code{t}. Default is \code{NULL}.
#' @param bw_grid \code{vector (numeric)}. Bandwidth grid for selecting the optimal smoothing parameter for each pair (\code{s}, \code{t}). Defaults to \code{NULL}, which generates an exponential grid of \eqn{N \lambda}.
#' @param use_same_bw \code{logical}. Indicates whether the same bandwidth should be used for both \code{s} and \code{t}. Defaults to \code{FALSE}.
#' @param center \code{logical (TRUE or FALSE)}. Default \code{center = TRUE} and so the curves are centred when the autocovariance is estimated: \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))}.
#' Otherwise, the two parts \eqn{\mathbb{E}X_0(s)X_{\ell}(t)} and \eqn{\mu(s)\mu(t)} will be estimated separately.
#' The first part with a bandwidth obtained with \link{estimate_autocov_risk} and the second part with a bandwidth obtained with \link{estimate_mean_risk}.
#' @param correct_diagonal \code{logical (TRUE or FALSE)}. Indicates whether the diagonal of the covariance should be corrected when \code{lag=0}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov". Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{s : The first argument of the autocovariance function.}
#'            \item{t : The second argument of the autocovariance function.}
#'            \item{optbw_s : The optimal bandwidth for the first argument of the autocovariance function. If \code{use_same_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
#'            \item{optbw_t : The optimal bandwidth for the second argument of the autocovariance function.}
#'            \item{Hs : The estimates of the local exponent for each \code{s}. It corresponds to \eqn{H_s}.}
#'            \item{Ls : The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
#'            \item{Ht : The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}.}
#'            \item{Lt : The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
#'            \item{PNs : The number of curves used to estimate the mean at \code{s}. It corresponds to \eqn{P_N(s;h)}.}
#'            \item{muhat_s : The estimates of the mean at \code{s}. It corresponds to \eqn{\widehat{\mu}_N(s;h)}.}
#'            \item{PNt : The number of curves used to estimate the mean at \code{t}. It corresponds to \eqn{P_N(t;h)}.}
#'            \item{muhat_t : The estimates of the mean at \code{t}. It corresponds to \eqn{\widehat{\mu}_N(t;h)}.}
#'            \item{PNl : The number of curves used to estimate the autocovariance at \code{(s,t)}. It corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
#'            \item{autocov : The estimates of the covariance/autocovariance.}
#'         }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_autocov_risk()].
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' #' # Load data
#' data("data_far")
#'
#' # Estimate adaptive lag-1 autocovariance
#' dt_autocov <- estimate_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
#'   optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
#'   use_same_bw = FALSE, center = TRUE, correct_diagonal = FALSE,
#'   kernel_name = "epanechnikov")
#' dt_autocov[, .(s, t, lag, PNl, autocov)]
#'
#' # Estimate adaptive covariance
#'
#' dt_cov <- estimate_autocov(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 0,
#'   optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
#'   use_same_bw = FALSE, center = TRUE, correct_diagonal = TRUE,
#'   kernel_name = "epanechnikov")
#' dt_cov[, .(s, t, lag, PNl, autocov)]
#' }
#'
estimate_autocov <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             s = c(1/5, 2/5, 4/5),
                             t = c(1/4, 1/2, 3/4),
                             lag = 1,
                             optbw_s = NULL, optbw_t = NULL,
                             bw_grid = NULL,
                             use_same_bw = FALSE,
                             center = TRUE,
                             correct_diagonal = TRUE,
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

  # Estimate autocovariance using C++ function
  mat_autocov <- estimate_autocov_cpp(
    data = data, s = s, t = t, lag = lag,
    optbw_s = optbw_s, optbw_t = optbw_t, bw_grid = bw_grid,
    use_same_bw = use_same_bw, center = center,
    correct_diagonal = correct_diagonal, kernel_name = kernel_name)
  dt_autocov <- data.table::as.data.table(mat_autocov)

  data.table::setnames(
    x = dt_autocov,
    new = c("s", "t", "optbw_s", "optbw_t", "Hs", "Ls", "Ht", "Lt",
            "PNs", "muhat_s", "PNt", "muhat_t", "PNl", "autocov"))

  return(dt_autocov)
}

# Autocovariance function estimator : Rubìn et Paranaretos (2020) ----
# Following the Rubìn and Panaretos Equation (B.7), we define Spq_fun and Qpq_fun

#' \eqn{S_{pq}^{(\ell)}}, (\eqn{\ell \leq 0}) function. See \insertCite{rubin2020;textual}{adaptiveFTS} Equation (B.7)
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
#' @param p \code{numeric (integer)}. It is used as exponent.
#' @param q \code{numeric (integer)}. It is used as exponent.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
#' Default \code{smooth_ker = epanechnikov}.
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @return A \code{numeric} scalar.
#' @export
#'
.Spq_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                     h, smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1)) & length(s) == 1))
    stop("'s' must be a numeric scalar value between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))  & length(t) == 1))
    stop("'t' must be a numeric scalar value between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")
  if ((any(p < 0)| (length(p) > 1) | any(p - floor(p) > 0)) |
      any(q < 0)| (length(q) > 1) | any(q - floor(q) > 0))
    stop("'p' and 'q' must be positive integers.")
  if (! (methods::is(h, "numeric") & all(data.table::between(h, 0, 1))  & length(h) == 1))
    stop("'h' must be a numeric scalar value between 0 and 1.")

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

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
  res <- (((xthj_vec - t) / h ) ** p) * (((xtk_vec - s) / h ) ** q) *
    (1 / (h ** 2)) * smooth_ker((xthj_vec - t) / h) * smooth_ker((xtk_vec - s) / h )
  Spq_sum <- sum(res)
  rm(res) ; gc()
  Spq <- Spq_sum / (N - lag)
  return(Spq)
}

#' \eqn{Q_{pq}^{(\ell)}}, (\eqn{\ell \leq 0}) function. See \insertCite{rubin2020;textual}{adaptiveFTS} Equation (B.7)
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
#' @param p \code{numeric (integer)}. It is used as exponent.
#' @param q \code{numeric (integer)}. It is used as exponent.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param dt_mean_rp \code{data.table}. It contains the estimates of the mean function at each observation point for each curve.
#' The name of the curve identification column must be \code{id_curve}, the observation points column \code{tobs} and the mean estimates column \code{muhat_RP}.
#' Default \code{dt_mean_rp = NULL} and so it will be estimated.
#' @param optbw_mean \code{numeric (positive scalar)}. Optimal bandwidth for the mean function estimator.
#' It is \code{NULL} if \code{dt_mean_rp} is not \code{NULL}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
#' Default \code{smooth_ker = epanechnikov}.
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @return A \code{numeric} scalar.
#' @export
#'
.Qpq_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                     h, dt_mean_rp = NULL, optbw_mean = NULL, smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") && all(s > 0 & s <= 1) && length(s) == 1))
    stop("'s' must be a numeric scalar value between 0 and 1.")
  if (! (methods::is(t, "numeric") && all(t > 0 & t <= 1)  && length(t) == 1))
    stop("'t' must be a numeric scalar value between 0 and 1.")
  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")
  if ((any(p < 0)| (length(p) > 1) | any(p - floor(p) > 0)) |
      any(q < 0)| (length(q) > 1) | any(q - floor(q) > 0))
    stop("'p' and 'q' must be positive integers.")
  if (! (methods::is(h, "numeric") && all(h > 0 & h < 1)  && length(h) == 1))
    stop("'h' must be a numeric scalar value between 0 and 1.")
  if (is.null(dt_mean_rp)) {
    if (is.null(optbw_mean)) {
      stop("If 'dt_mean_rp' is NULL, then optbw_mean can not be NULL")
    } else if (! (methods::is(optbw_mean, "numeric") && all(optbw_mean > 0 & optbw_mean < 1) && length(optbw_mean) == 1)) {
      stop("'optbw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else if (! (data.table::is.data.table(dt_mean_rp) & all(c("id_curve", "tobs", "muhat_RP") %in% colnames(dt_mean_rp)))) {
      stop("'dt_mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate mean function is it is NULL
  if (is.null(dt_mean_rp)) {
    dt_mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = dt_mean_rp[, tobs], h = optbw_mean, smooth_ker = smooth_ker)
    dt_mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
    rm(dt_mean) ; gc()
  } else {
    dt_mean_rp <- dt_mean_rp[order(id_curve)]
  }

  # Extract observation points and observed points
  data <- data[order(id_curve)]
  Tn <- data[id_curve %in% 1:(N - lag), tobs]
  Tn_plus_lag <- data[id_curve %in% (1 + lag):N, tobs]
  Yn <- data[id_curve %in% 1:(N - lag), X]
  Yn_plus_lag <- data[id_curve %in% (1 + lag):N, X]
  rm(data); gc()

  # Extract mean function estimates
  dt_mean_rp <- dt_mean_rp[order(id_curve)]
  muhat_Tn <- dt_mean_rp[id_curve %in% 1:(N - lag), muhat_RP]
  muhat_Tn_plus_lag <- dt_mean_rp[id_curve %in% (1 + lag):N, muhat_RP]
  rm(dt_mean_rp) ; gc()
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
  res <- Gth * (((xthj_vec - t) / h ) ** p) * (((xtk_vec - s) / h ) ** q) *
    (1 / (h ** 2)) * smooth_ker((xthj_vec - t) / h) * smooth_ker((xtk_vec - s) / h )

  Qpq_sum <- sum(res)
  rm(res) ; gc()
  Qpq <- Qpq_sum / (N - lag)
  return(Qpq)
}

#' Estimate lag-\eqn{\ell} (\eqn{\ell \leq 0}) autocovariance function using \insertCite{rubin2020;textual}{adaptiveFTS} method
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param dt_mean_rp \code{data.table}. It contains the estimates of the mean function at each observation point for each curve.
#' The name of the curve identification column must be \code{id_curve}, the observation points column \code{tobs} and the mean estimates column \code{muhat_RP}.
#' Default \code{dt_mean_rp = NULL} and so it will be estimated.
#' @param optbw_mean \code{numeric (positive scalar)}. Optimal bandwidth for the mean function estimator.
#' It is \code{NULL} if \code{dt_mean_rp} is not \code{NULL}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
#' Default \code{smooth_ker = epanechnikov}.
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#'
#' @references
#' \insertAllCited{}
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{s :}{ The first argument of the autocovariance function.}
#'            \item{t :}{ The second argument of the autocovariance function.}
#'            \item{lag :}{ The lag of the autocovariance. It corresponds to \eqn{\ell \leq 0}.}
#'            \item{optbw_mean :}{ The optimal bandwidth for the mean function estimator.}
#'            \item{h :}{ The bandwidth used to estimate the lag-\eqn{\ell}, \eqn{\ell \leq 0} autocovariance function}
#'            \item{autocovhat_rp :}{ The estimates of the lag-\eqn{\ell} autocovariance function for each (\code{s}, \code{t}) using Rubìn and Panaretos (2020) method.}
#'         }
#' @export
#'
estimate_autocov_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
                                lag = 1, h, optbw_mean = NULL, dt_mean_rp = NULL,
                                smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") && all(s > 0 & s <= 1)))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") && all(t > 0 & t <= 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! (methods::is(h, "numeric") && all(h > 0 & h < 1) && length(h) == 1))
    stop("'h' must be a numeric scalar value between 0 and 1.")
  if (is.null(dt_mean_rp)) {
    if (is.null(optbw_mean)) {
      stop("If 'dt_mean_rp' is NULL, then optbw_mean can not be NULL")
    } else if (! (methods::is(optbw_mean, "numeric") && all(optbw_mean > 0 & optbw_mean < 1) && length(optbw_mean) == 1)) {
      stop("'optbw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else if (! (data.table::is.data.table(dt_mean_rp) && all(c("id_curve", "tobs", "muhat_RP") %in% colnames(dt_mean_rp)))) {
      stop("'dt_mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }
  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st) ; gc()

  # Estimate mean function is it is NULL
  if (is.null(dt_mean_rp)) {
    dt_mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = dt_mean_rp[, tobs], h = optbw_mean, smooth_ker = smooth_ker)
    dt_mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
  } else {
    dt_mean_rp <- dt_mean_rp[order(id_curve)]
  }

  # Calculate S_{pq} and Q_{pq}
  autocov_vec <- mapply(function(si, ti, h, lag, optbw_mean, dt_mean_rp, data, ker){
    # Calculate S_{pq} and A_1^{(\ell)},A_2^{(\ell)}, A_3^{(\ell)}
    S00 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 0, h = h, smooth_ker = ker)
    S01 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 1, h = h, smooth_ker = ker)
    S02 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 0, q = 2, h = h, smooth_ker = ker)
    S10 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 1, q = 0, h = h, smooth_ker = ker)
    S11 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 1, q = 1, h = h, smooth_ker = ker)
    S20 <- .Spq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    s = si, t = ti, lag = lag, p = 2, q = 0, h = h, smooth_ker = ker)

    # calculate A_1^{(\ell)},A_2^{(\ell)}, A_3^{(\ell)} and B^{(\ell)}
    A1 <- S20 * S02 - (S11 ** 2)
    A2 <- S10 * S02 - S01 * S11
    A3 <- S01 * S20 - S10 * S11
    B <- A1 * S00 - A2 * S10 - A3 * S01

    # Calculate Q_{pq}
    Q00 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 0, q = 0, h = h, dt_mean_rp = dt_mean_rp, optbw_mean = optbw_mean, smooth_ker = ker)
    Q10 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 1, q = 0, h = h, dt_mean_rp = dt_mean_rp, optbw_mean = optbw_mean, smooth_ker = ker)
    Q01 <- .Qpq_fun(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", s = si, t = ti,
                    lag = lag, p = 0, q = 1, h = h, dt_mean_rp = dt_mean_rp, optbw_mean = optbw_mean, smooth_ker = ker)

    # estimate autocovariance
    R <- (A1 * Q00 - A2 * Q10 - A3 * Q01) / B

    return(R)
  }, si = s, ti = t, MoreArgs = list(h = h, lag = lag, data = data, optbw_mean = optbw_mean,
                                     dt_mean_rp = dt_mean_rp, ker = smooth_ker))
  dt_res <- data.table::data.table("s" = s, "t" = t, "lag" = lag, "optbw_mean" = optbw_mean, "autocovhat_rp" = autocov_vec)
  return(dt_res)
}

#' Bandwidth estimation using cross-validation for the \insertCite{rubin2020;textual}{adaptiveFTS} autocovariance function estimator.
#'
#' @inheritParams format_data
#' @param Kfold \code{integer (positive)}. Number of fold for the cross-validation.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid.
#' @param dt_mean_rp \code{data.table}. It contains the estimates of the mean function at each observation point for each curve.
#' The name of the curve identification column must be \code{id_curve}, the observation points column \code{tobs} and the mean estimates column \code{muhat_RP}.
#' Default \code{dt_mean_rp = NULL} and so it will be estimated.
#' @param optbw_mean \code{numeric (positive scalar)}. Optimal bandwidth for the mean function estimator.
#' It is \code{NULL} if \code{dt_mean_rp} is not \code{NULL}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{cv_error :}{ The estimates of the Cross-Validation error for each \code{h}.}
#'         }
#' @export
#'
#' @import data.table
#' @import Rdpack
#' @importFrom methods is
#' @importFrom caret createFolds
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [estimate_mean_rp()], [estimate_mean_bw_rp()], [estimate_autocov_rp()]
#'
#'
estimate_autocov_bw_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                   Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
                                   optbw_mean = NULL, dt_mean_rp = NULL, smooth_ker = epanechnikov){
  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(Kfold < 0)| (length(Kfold) > 1) | any(Kfold - floor(Kfold) > 0) | any(N <= Kfold))
    stop("'Kfold' must be a positive integer lower than the number of curves.")
  if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
    stop("'bw_grid' must be a vector of positive values between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (is.null(dt_mean_rp)) {
    if (is.null(optbw_mean)) {
      stop("If 'dt_mean_rp' is NULL, then optbw_mean can not be NULL")
    } else if (! (methods::is(optbw_mean, "numeric") & all(optbw_mean > 0 & optbw_mean <= 1)  & length(optbw_mean) == 1)) {
      stop("'optbw_mean' must be a numeric scalar value between 0 and 1.")
    }
  } else {
    if (! (data.table::is.data.table(dt_mean_rp) & all(c("id_curve", "tobs", "muhat_RP") %in% colnames(dt_mean_rp))))
      stop("'dt_mean_rp' must be a data.table containing the columns : 'id_curve', 'tobs' and 'muhat_RP'.")
  }

  # Estimate mean function is it is NULL
  if (is.null(dt_mean_rp)) {
    dt_mean_rp <- data[order(tobs), list(id_curve, tobs)]
    dt_mean <- estimate_mean_rp(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = dt_mean_rp[, tobs], h = optbw_mean, smooth_ker = smooth_ker)
    dt_mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
    rm(dt_mean) ; gc()
  } else {
    dt_mean_rp <- dt_mean_rp[order(id_curve)]
  }

  # Create Kfold folds
  fold <- caret::createFolds(y = unique(data[, id_curve]), k = Kfold, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(parallel::mclapply(bw_grid, function(BR0, data, dt_mean_rp, fold, kernel_smooth){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, dt_mean_rp, BR0, kernel_smooth){
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
        dt_mean_test <- dt_mean_rp[id_curve %in% unlist(f)]
        muhat <- dt_mean_test[order(tobs), muhat_RP]
        muhat_grid <- expand.grid(muhat_ti = muhat, muhat_tj = muhat)
        muhat_ti <- Yn_grid$muhat_ti
        muhat_tj <- Yn_grid$muhat_tj

        rm(Tn, Yn, Tn_grid, Yn_grid, muhat_grid, muhat, dt_mean_test) ; gc()

        # Estimation of mean on fold\f and test on f
        dt_autocov <- estimate_autocov_rp(
          data = dt_train, idcol = "id_curve", tcol = "tobs", ycol = "X",
          s = xti, t = xtj, lag = 0, h = BR0, optbw_mean = optbw_mean,
          dt_mean_rp = dt_mean_rp, smooth_ker = kernel_smooth)

        # Calculate the error
        Sqerror <- ((Yti - muhat_ti) * (Ytj - muhat_tj) - dt_autocov[, autocovhat_rp]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, dt_mean_rp = dt_mean_rp, BR0 = BR0, kernel_smooth = kernel_smooth, simplify = TRUE),
      error = function(e){
        message("Error in estimating the autocovariance function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold, na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("h" = BR0, "cv_error" = cv_err)
    return(dt_res)

  }, mc.cores = 20, data = data, dt_mean_rp = dt_mean_rp, fold = fold, kernel_smooth = smooth_ker))

  return(dt_bw)
}

