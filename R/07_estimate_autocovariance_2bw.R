#' Estimate the risk of the lag-\eqn{\ell} (\eqn{\ell > 0}) autocovariance function with two bandwidths.
#'
#' @inheritParams format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameters (\eqn{h_s}, \eqn{h_t}) are selected for each pair (\code{s}, \code{t}).
#' It can be \code{NULL}, in which case it will be defined as an exponential grid of \eqn{N\times\lambda}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
#' Default \code{smooth_ker = epanechnikov}.
#' @param Hs \code{vector (numeric)}. The estimates of the local exponent for each \code{s}.
#' Default \code{Hs = NULL} and thus it will be estimated.
#' @param Ls \code{vector (numeric)}. The estimates of the Hölder constant for each \code{s}.
#' It corresponds to \eqn{L_s^2}. Default \code{Ls = NULL} and thus it will be estimated.
#' @param Ht \code{vector (numeric)}. The estimates of the local exponent for each \code{t}.
#' Default \code{Ht = NULL} and thus it will be estimated.
#' @param Lt \code{vector (numeric)}. The estimates of the Hölder constant for each \code{t}.
#' It corresponds to \eqn{L_t^2}. Default \code{Lt = NULL} and thus it will be estimated.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{s} and \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{s :}{ The first argument of the autocovariance function.}
#'            \item{t :}{ The second argument of the autocovariance function.}
#'            \item{hs :}{ The candidate bandwidth for the argument \code{s}.}
#'            \item{ht :}{ The candidate bandwidth for the argument \code{t}.}
#'            \item{Hs :}{ The estimates of the local exponent for each \code{s}. It corresponds to \eqn{H_s}}
#'            \item{Ls :}{ The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
#'            \item{Ht :}{ The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}}
#'            \item{Lt :}{ The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
#'            \item{bias_term :}{ The bias term of the risk function.}
#'            \item{varriance_term :}{ The variance term of the risk function.}
#'            \item{dependence_term :}{ The dependence term of the risk function.}
#'            \item{autocov_risk :}{ The estimates of the risk function of the autocovariance function.}
#'         }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_XsXt_autocov()].
#'
#' @import data.table
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Generate a sample path of FTS
#' dt_far <- simulate_far(N = 50, lambda = 70,
#'                        tdesign = "random",
#'                        Mdistribution = rpois,
#'                        tdistribution = runif,
#'                        tcommon = NULL,
#'                        hurst_fun = hurst_logistic,
#'                        L = 4,
#'                        far_kernel = get_real_data_far_kenel,
#'                        far_mean = get_real_data_mean,
#'                        int_grid = 100L,
#'                        burnin = 100L,
#'                        remove_burnin = TRUE)
#'
#' # Add noise
#' dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]
#'
#' # Estimate risk function
#' dt_autocov_risk <- estimate_autocov_risk(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   s = c(1/5, 2/5, 4/5),
#'   t = c(1/4, 1/2, 3/4),
#'   lag = 3,
#'   bw_grid = seq(0.005, 0.15, len = 45),
#'   Delta = NULL, h = NULL, smooth_ker = epanechnikov
#' )
#'
#' # Plot mean risk function
#' dt_dcast <- data.table::dcast(data = dt_autocov_risk,
#'                               formula = h ~ s + t ,
#'                               value.var = "autocov_risk")
#'
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.2, 0.25)" = `0.2_0.25`)],
#'                       main = "lag = 3 - (s, t) = (0.2, 0.25)",
#'                       xlab = "h",
#'                       ylab = "risk function"),
#'     dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.4, 0.5)" = `0.4_0.5`)],
#'                       main = "lag = 3 - (s, t) = (0.4, 0.5)",
#'                       xlab = "h",
#'                       ylab = "risk function"),
#'     dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.8, 0.75)" = `0.8_0.75`)],
#'                       main = "lag = 3 - (s, t) = (0.8, 0.75)",
#'                       xlab = "h",
#'                       ylab = "risk function")
#'   ),
#'   nrow = 3
#' )
#'
#' }
#'
estimate_autocov_risk_2bw <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                      s = c(1/5, 2/5, 4/5),
                                      t = c(1/4, 1/2, 3/4),
                                      lag = 1,
                                      bw_grid = seq(0.005, 0.15, len = 20),
                                      smooth_ker = epanechnikov,
                                      Hs = NULL, Ls = NULL,
                                      Ht = NULL, Lt = NULL,
                                      Delta = NULL, h = NULL){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  # Control on local regularity parameters
  if (((!is.null(Hs)) & length(Hs) != length(s)) | ((!is.null(Ls)) & length(Ls) != length(s)))
    stop("If 'Hs' or 'Ls' is not NULL, it must be the same length as 's'.")
  if ( ((!is.null(Ht)) & length(Ht) != length(t)) | ((!is.null(Lt)) & length(Lt) != length(t)))
    stop("If 'Ht' or 'Lt' is not NULL, it must be the same length as 't'.")
  if ((!is.null(Hs) & ! methods::is(Hs, "numeric")) |
      (!is.null(Ls) & ! methods::is(Ls, "numeric")) |
      (!is.null(Ht) & ! methods::is(Ht, "numeric")) |
      (!is.null(Lt) & ! methods::is(Lt, "numeric")))
    stop("If 'Hs', 'Ls', 'Ht' or 'Lt' is not NULL, it must be numeric.")

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

  if (any(lag < 0) | (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

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
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, smooth_ker = smooth_ker)
      h <- dt_optbw[, median(optbw)]
      rm(dt_optbw) ; gc()
    } else {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, smooth_ker = smooth_ker)
      h <- dt_optbw[, optbw]
      rm(dt_optbw) ; gc()
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Estimate local regularity parameters
  # This function controls the remaining arguments
  if (is.null(Hs) | is.null(Ls)) {
    dt_locreg_s <- estimate_locreg(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = s, Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    Hs <- dt_locreg_s[, Ht]
    Ls <- dt_locreg_s[, Lt]
    rm(dt_locreg_s) ; gc()
  }
  if (is.null(Ht) | is.null(Lt)) {
    dt_locreg_t <- estimate_locreg(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = t, Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    Ht <- dt_locreg_t[, Ht]
    Lt <- dt_locreg_t[, Lt]
    rm(dt_locreg_t) ; gc()
  }
  # Estimation of the observation error standard deviation
  ## If s and t are all equal OR not
  if (all(unique(s) %in% unique(t))) {
    dt_sigma <- estimate_sigma(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = sort(unique(s)))
    dt_sigma <- data.table::merge.data.table(
      x = data.table::merge.data.table(x = dt_sigma[, .("s" = t, "sig_error_s" = sig)], y = data.table::data.table("s" = s, "t" = t), by = "s"),
      y = data.table::merge.data.table(x = dt_sigma[, .("t" = t, "sig_error_t" = sig)], y = data.table::data.table("s" = s, "t" = t), by = "t"),
      by = c("s", "t")
    )
  } else {
    dt_sigma_s <- estimate_sigma(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = s)
    dt_sigma_t <- estimate_sigma(data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = t)
    dt_sigma <- data.table::merge.data.table(
      x = data.table::merge.data.table(x = dt_sigma_s[, .("s" = t, "sig_error_s" = sig)], y = data.table::data.table("s" = s, "t" = t), by = "s"),
      y = data.table::merge.data.table(x = dt_sigma_t[, .("t" = t, "sig_error_t" = sig)], y = data.table::data.table("s" = s, "t" = t), by = "t"),
      by = c("s", "t")
    )
    rm(dt_sigma_s, dt_sigma_t) ; gc()
  }

  # Estimation of the centered second order moment using the presmoothing bandwidth
  ## Smooth curves at s and at t
  dt_Xhat <- data.table::rbindlist(
    lapply(1:N, function(curve_index, s, t, presmooth_bw_s, presmooth_bw_t, kernel_smooth, data){
      Xhat_s <- estimate_nw(
        y = data[id_curve %in% curve_index][order(tobs), X],
        t = data[id_curve %in% curve_index][order(tobs), tobs],
        tnew = s,
        h = presmooth_bw_s[curve_index],
        smooth_ker = kernel_smooth)
      Xhat_s <- Xhat_s[, yhat]

      Xhat_t <- estimate_nw(
        y = data[id_curve %in% curve_index][order(tobs), X],
        t = data[id_curve %in% curve_index][order(tobs), tobs],
        tnew = t,
        h = presmooth_bw_t[curve_index],
        smooth_ker = kernel_smooth)
      Xhat_t <- Xhat_t[, yhat]

      dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t, "Xhat_s" = Xhat_s, "Xhat_t" = Xhat_t)
    }, s = s, t = t,
    presmooth_bw_s = h, presmooth_bw_t = h,
    kernel_smooth = smooth_ker, data = data))

  ##  Transform NaN to NA
  dt_Xhat[is.nan(Xhat_s), Xhat_s := NA]
  dt_Xhat[is.nan(Xhat_t), Xhat_t := NA]

  # Estimate mean
  dt_Xhat[, mu_s := mean(Xhat_s, na.rm = TRUE), by = "s"]
  dt_Xhat[, mu_t := mean(Xhat_t, na.rm = TRUE), by = "t"]

  dt_EX2 <- dt_Xhat[, list("EX2_centered_s" = mean((Xhat_s - mu_s) ** 2, na.rm = TRUE),
                           "EX2_centered_t" = mean((Xhat_t - mu_t) ** 2, na.rm = TRUE)),
                    by = c("s", "t")]
  rm(dt_Xhat) ; gc()

  # Estimation of the "weighted" empirical autocovariance of X_0(s)X_{\ell}(s)
  ## The bandwidth is taking as the presmoothing bandwidth

  dt_XsXt_autocov <- estimate_empirical_XsXt_autocov(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = s, t = t, cross_lag = lag,
    lag = 0:(N - lag - 1), h = h,
    smooth_ker = smooth_ker, center = TRUE)
  # If the variance is NaN, set it as 0
  dt_XsXt_autocov[is.nan(XsXt_autocov) & lag == 0, XsXt_autocov := 0]

  ## compute weights p_k(s,t;h)
  ## We assume assume h is the presmoothing bandwidth
  # dt_pi_raw <- data.table::rbindlist(lapply(1:N, function(curve_index, bw_grid, s, t, data){
  #   # Extract the observation points per curve
  #   Tn <- data[id_curve == curve_index, tobs]
  #
  #   dt_pi_hi <- data.table::rbindlist(lapply(bw_grid, function(hi, Tn){
  #     # \pi_n(s;h)
  #     pin_s <- outer(X = s, Y = Tn, FUN = function(si, Tnm, hi) as.numeric(abs(Tnm - si) <= hi), hi = hi)
  #     pin_s <- as.numeric(rowSums(pin_s) >= 1)
  #
  #     # \pi_n(t;h)
  #     pin_t <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, hi) as.numeric(abs(Tnm - ti) <= hi), hi = hi)
  #     pin_t <- as.numeric(rowSums(pin_t) >= 1)
  #
  #     # data.table return
  #     dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t, "pin_s" = pin_s, "pin_t" = pin_t, "hi" = hi)
  #   }, Tn = Tn))
  #   return(dt_pi_hi)
  # }, bw_grid = bw_grid, s = s, t = t, data = data))
  #
  # ## Calculate P_{N,\ell}(s,t,h)
  # dt_pin_s <- data.table::merge.data.table(
  #   x = data.table::data.table("id_curve" = 1:(N - lag), "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N)),
  #   y = dt_pi_raw[id_curve %in% 1:(N - lag), list(id_curve, hi, s, t, pin_s)],
  #   by = "id_curve")
  #
  # dt_pin_t <- data.table::merge.data.table(
  #   x = data.table::data.table("id_curve" = (1 + lag):N, "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N)),
  #   y = dt_pi_raw[id_curve %in% (1 + lag):N, list(id_curve, hi, s, t, pin_t)],
  #   by = "id_curve")
  # dt_PNl <- data.table::merge.data.table(x = dt_pin_s, y = dt_pin_t, by = c("id_lag", "s", "t", "hi"))
  # dt_PNl <- dt_PNl[, .("PNl" = sum(pin_s * pin_t)), by = c("s", "t", "hi")]
  # rm(dt_pin_s, dt_pin_t) ; gc()
  #
  # ## Compuet weights
  # dt_p <- data.table::rbindlist(lapply(1:(N - lag - 1), function(k, dt_pi_raw, dt_PNl){
  #   ## Extract lagged data
  #   dt_pi_s_i <- dt_pi_raw[id_curve %in% 1:(N - lag - k), .(id_curve, hi, s, t, "pin_s_i" = pin_s)]
  #   dt_pi_s_i_plus_k <- dt_pi_raw[id_curve %in% (1 + k):(N - lag), .(id_curve, hi, s, t, "pin_s_i_plus_k" = pin_s)]
  #
  #   dt_pi_t_i_plus_ell <- dt_pi_raw[id_curve %in% (1 + lag):(N - k), .(id_curve, hi, s, t, "pin_t_i_plus_ell" = pin_t)]
  #   dt_pi_t_i_plus_ell_plus_k <- dt_pi_raw[id_curve %in% (1 + k + lag):N, .(id_curve, hi, s, t, "pin_t_i_plus_ell_plus_k" = pin_t)]
  #
  #   ## Add lag id to \pi_i(s;h)
  #   dt_id_lag_s_i <- data.table::data.table(
  #     "id_curve" = sort(unique(dt_pi_s_i[, id_curve])),
  #     "id_lag" = paste0(1:(N - lag - k), "_", (1 + k):(N - lag), "_", (1 + lag):(N - k), "_", (1 + k + lag):N)
  #   )
  #   dt_pi_s_i <- data.table::merge.data.table(x = dt_id_lag_s_i, y = dt_pi_s_i, by = "id_curve")
  #   rm(dt_id_lag_s_i) ; gc()
  #
  #   ## Add lag id to \pi_{i + k}(s;h)
  #   dt_id_lag_s_i_plus_k <- data.table::data.table(
  #     "id_curve" = sort(unique(dt_pi_s_i_plus_k[, id_curve])),
  #     "id_lag" = paste0(1:(N - lag - k), "_", (1 + k):(N - lag), "_", (1 + lag):(N - k), "_", (1 + k + lag):N)
  #   )
  #   dt_pi_s_i_plus_k <- data.table::merge.data.table(x = dt_id_lag_s_i_plus_k, y = dt_pi_s_i_plus_k, by = "id_curve")
  #   rm(dt_id_lag_s_i_plus_k) ; gc()
  #
  #   ## Add lag id to \pi_{i + \ell}(t;h)
  #   dt_id_lag_t_i_plus_ell <- data.table::data.table(
  #     "id_curve" = sort(unique(dt_pi_t_i_plus_ell[, id_curve])),
  #     "id_lag" = paste0(1:(N - lag - k), "_", (1 + k):(N - lag), "_", (1 + lag):(N - k), "_", (1 + k + lag):N)
  #   )
  #   dt_pi_t_i_plus_ell <- data.table::merge.data.table(x = dt_id_lag_t_i_plus_ell, y = dt_pi_t_i_plus_ell, by = "id_curve")
  #   rm(dt_id_lag_t_i_plus_ell) ; gc()
  #
  #   ## Add lag id to \pi_{i + \ell + k}(t;h)
  #   dt_id_lag_t_i_plus_ell_plus_k <- data.table::data.table(
  #     "id_curve" = sort(unique(dt_pi_t_i_plus_ell_plus_k[, id_curve])),
  #     "id_lag" = paste0(1:(N - lag - k), "_", (1 + k):(N - lag), "_", (1 + lag):(N - k), "_", (1 + k + lag):N)
  #   )
  #   dt_pi_t_i_plus_ell_plus_k <- data.table::merge.data.table(x = dt_id_lag_t_i_plus_ell_plus_k, y = dt_pi_t_i_plus_ell_plus_k, by = "id_curve")
  #   rm(dt_id_lag_t_i_plus_ell_plus_k) ; gc()
  #
  #   ## All four weight dt
  #   dt_pi <- data.table::merge.data.table(
  #     x = data.table::merge.data.table(x = dt_pi_s_i, y = dt_pi_s_i_plus_k, by = c("id_lag", "s", "t", "hi")),
  #     y = data.table::merge.data.table(x = dt_pi_t_i_plus_ell, y = dt_pi_t_i_plus_ell_plus_k, by = c("id_lag", "s", "t", "hi")),
  #     by = c("id_lag", "s", "t", "hi")
  #   )
  #   rm(dt_pi_s_i, dt_pi_s_i_plus_k, dt_pi_t_i_plus_ell, dt_pi_t_i_plus_ell_plus_k) ; gc()
  #
  #   ## Calculate the numerator of p_k
  #   dt_pk <- dt_pi[, .("pk" = sum(pin_s_i * pin_s_i_plus_k * pin_t_i_plus_ell * pin_t_i_plus_ell_plus_k)), by = c("s", "t", "hi")]
  #
  #   dt_pk <- data.table::merge.data.table(x = dt_pk, y = dt_PNl, by = c("s", "t", "hi"))
  #   dt_pk[, pk := pk / PNl, by = c("s", "t", "hi")]
  #   dt_pk[PNl == 0, pk := 0]
  #   dt_pk[, lag := k]
  # }, dt_pi_raw = dt_pi_raw, dt_PNl = dt_PNl))
  #
  #   dt_long_run_var <- data.table::merge.data.table(x = dt_XsXt_autocov[lag > 0], y = dt_p, by = c("s", "t", "lag"))
  #   dt_long_run_var <- dt_long_run_var[! is.nan(XsXt_autocov), list("long_run_var" = sum(2 * XsXt_autocov * pk)), by = c("s", "t", "hi")]
  #   dt_dependence_coef <- data.table::merge.data.table(
  #     x = dt_long_run_var,
  #     y = dt_XsXt_autocov[lag == 0][, .(s, t, "XsXt_var" = XsXt_autocov)],
  #     by = c("s", "t")
  #   )
  #   dt_dependence_coef[, dependence_coef := XsXt_var + long_run_var]
  dt_long_run_var <- dt_XsXt_autocov[(! is.nan(XsXt_autocov)) & (lag != 0), list("long_run_var" = sum(2 * abs(XsXt_autocov))), by = c("s", "t")]
  dt_dependence_coef <- data.table::merge.data.table(
    x = dt_long_run_var,
    y = dt_XsXt_autocov[lag == 0][, .(s, t, "XsXt_var" = XsXt_autocov)],
    by = c("s", "t")
  )
  dt_dependence_coef[, dependence_coef := XsXt_var + long_run_var]

  # Take the bandwidth grid
  bwg <- expand.grid("bw_grid_s" = bw_grid, "bw_grid_t" = bw_grid)
  bw_grid_s <- bwg$bw_grid_s
  bw_grid_t <- bwg$bw_grid_t
  rm(bwg)

  # Estimate the risk function
  dt_autocov_risk <- data.table::rbindlist(
    lapply(1:length(bw_grid_s), function(idx_bw, s, t, lag, Hs, Ls, Ht, Lt, bw_grid_s, bw_grid_t,
                                         kernel_smooth, dt_sigma, N, data, dt_EX2, dt_dependence_coef){
      hs <- bw_grid_s[idx_bw]
      ht <- bw_grid_t[idx_bw]
      # compute quantities to estimate estimate the risk
      dt_risk_raw <- data.table::rbindlist(lapply(1:N, function(curve_index, hs, ht, s, t, Hs, Ls, Ht, Lt, kernel_smooth, data){
        # Extract the observation points per curve
        # and compute the weight of the estimator
        ## By convention NaN = 0
        Tn <- data[id_curve == curve_index, tobs]

        ## For the argument s
        ker_s <- outer(X = s, Y = Tn, FUN = function(si, Tnm, Ker, h) Ker((Tnm - si) / h), Ker = kernel_smooth, h = hs)
        wker_s <- apply(ker_s, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
        wker_s <- t(wker_s)

        ### Take the maximum and c_n(s,h_s), the sum of the weight
        wmax_s <- apply(wker_s, 1, max)
        cn_s <- apply(wker_s, 1, function(r) sum(abs(r)))

        Tn_s_2H_s <- outer(
          X = 1:length(s), Y = Tn,
          FUN = function(sid, Tnm, s, Hs, h) abs((Tnm - s[sid]) / h) ** (2 * Hs[sid]),
          s = s, Hs = Hs, h = hs
        )

        bn2H_s <- diag(abs(wker_s) %*% t(Tn_s_2H_s))

        ### \pi_n(s;h)
        pin_s <- outer(X = s, Y = Tn, FUN = function(si, Tnm, h) as.numeric(abs(Tnm - si) <= h), h = hs)
        pin_s <- as.numeric(rowSums(pin_s) >= 1)

        ## For argument t
        ### By convention NaN = 0
        ker_t <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, Ker, h) Ker((Tnm - ti) / h), Ker = kernel_smooth, h = ht)
        wker_t <- apply(ker_t, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
        wker_t <- t(wker_t)

        ### Take the maximum and c_n(t,h), the sum of the weight
        wmax_t <- apply(wker_t, 1, max)
        cn_t <- apply(wker_t, 1, function(r) sum(abs(r)))

        Tn_t_2H_t <- outer(
          X = 1:length(t), Y = Tn,
          FUN = function(tid, Tnm, t, Ht, h) abs((Tnm - t[tid]) / h) ** (2 * Ht[tid]),
          t = t, Ht = Ht, h = ht
        )

        bn2H_t <- diag(abs(wker_t) %*% t(Tn_t_2H_t))

        # \pi_n(t;h)
        pin_t <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, h) as.numeric(abs(Tnm - ti) <= h), h = ht)
        pin_t <- as.numeric(rowSums(pin_t) >= 1)

        # data.table return
        dt_res <- data.table::data.table(
          "id_curve" = curve_index, "s" = s, "t" = t,
          "wmax_s" = wmax_s, "wmax_t" = wmax_t,
          "cn_s" = cn_s, "cn_t" = cn_t,
          "bn2H_s" = bn2H_s, "bn2H_t" = bn2H_t,
          "pin_s" = pin_s, "pin_t" = pin_t
        )
        return(dt_res)
      }, hs = hs, ht = ht, s = s, t = t, Hs = Hs, Ht = Ht, kernel_smooth = kernel_smooth, data = data))

      # Split and merge by lag
      ## The argument s is associated to the curves n = 1,..., N-lag
      dt_risk_s <- dt_risk_raw[, list(id_curve, s, t, wmax_s, cn_s, bn2H_s, pin_s)]
      dt_risk_s <- dt_risk_s[id_curve %in% 1:(N - lag)]
      dt_id_lag_s <- data.table::data.table(
        "id_curve" = sort(unique(dt_risk_s[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_risk_s <- data.table::merge.data.table(x = dt_id_lag_s, y = dt_risk_s, by = "id_curve")

      ## The argument t is associated to the curves n = 1 + lag,..., N
      dt_risk_t <- dt_risk_raw[, list(id_curve, s, t, wmax_t, cn_t, bn2H_t, pin_t)]
      dt_risk_t <- dt_risk_t[id_curve %in% (1 + lag):N]
      dt_id_lag_t <- data.table::data.table(
        "id_curve" = sort(unique(dt_risk_t[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_risk_t <- data.table::merge.data.table(x = dt_id_lag_t, y = dt_risk_t, by = "id_curve")

      ## Merge
      dt_risk_merge <- data.table::merge.data.table(
        x = dt_risk_s,
        y = dt_risk_t,
        by = c("id_lag", "s", "t"))
      dt_risk_merge <- dt_risk_merge[order(id_curve.x)]

      ## Clean unnecessary object
      rm(dt_id_lag_s, dt_id_lag_t, dt_risk_s, dt_risk_t) ; gc()

      # Compute P_N(s, t; h_s, h_t)
      dt_risk_merge[, PNl := sum(pin_s * pin_t), by = c("s", "t")]

      # compute \mathbb B(s, h_s, 2H_s|t,h, 0) and \mathbb B(t, h_t, 2H_t|s, h_s, \ell)
      dt_risk_merge[, Bs := sum(pin_s * pin_t * cn_s * bn2H_s) / PNl, by = c("s", "t")]
      dt_risk_merge[, Bt := sum(pin_s * pin_t * cn_t * bn2H_t) / PNl, by = c("s", "t")]

      # Compute \mathbb V_{\gamma, 1}(s, t; h) and \mathbb V_{\gamma, 2}(s, t; h)
      dt_risk_merge[, Vgamma1 := sum(pin_s * pin_t * cn_s * wmax_s) / PNl ** 2, by = c("s", "t")]
      dt_risk_merge[, Vgamma2 := sum(pin_s * pin_t * cn_t * wmax_t) / PNl ** 2, by = c("s", "t")]
      dt_risk_merge[, Vgamma := sum(pin_s * pin_t * cn_s * cn_t * wmax_s * wmax_t) / PNl ** 2, by = c("s", "t")]

      # Extract the elements by (s,t)
      dt_risk <- unique(dt_risk_merge[, list(s, t, Bs, Bt, Vgamma1, Vgamma2, Vgamma, PNl)])

      # Add second order moment and error standard deviation
      dt_risk <- data.table::merge.data.table(x = dt_risk, y = dt_EX2, by = c("s", "t"))
      dt_risk <- data.table::merge.data.table(x = dt_risk, y = dt_sigma, by = c("s", "t"))
      dt_risk <- dt_risk[order(s,t)]

      ## Biais term
      bias_term <- 4 * dt_risk[, EX2_centered_t] * Ls * (hs ** (2 * Hs)) * dt_risk[, Bs] +
        4 * dt_risk[, EX2_centered_s] * Lt * (ht ** (2 * Ht)) * dt_risk[, Bt]

      ## Variance term
      varriance_term <- 4 * (dt_risk[, sig_error_s] ** 2) * dt_risk[, EX2_centered_t] * dt_risk[, Vgamma1] +
        4 * (dt_risk[, sig_error_t] ** 2) * dt_risk[, EX2_centered_s] * dt_risk[, Vgamma2] +
        4 * (dt_risk[, sig_error_s] ** 2) * (dt_risk[, sig_error_t] ** 2) * dt_risk[, Vgamma]

      ## Dependence term
      # dependence_coef <- dt_dependence_coef[hi == h][order(s,t), dependence_coef]
      dependence_coef <- dt_dependence_coef[order(s,t), dependence_coef]
      dependence_coef <- abs(dependence_coef)
      dependence_term <- dependence_coef  / dt_risk[order(s, t), PNl]

      ## Final risk function
      autocov_risk <- 2 * (bias_term + varriance_term + dependence_term)

      # Result to returned
      dt_res <- data.table::data.table("s" = s, "t" = t, "hs" = hs, "ht" = ht, "Hs" = Hs, "Ls" = Ls, "Ht" = Ht, "Lt" = Lt,
                                       "bias_term" = bias_term, "varriance_term" = varriance_term,
                                       "dependence_term" = dependence_term, "autocov_risk" = autocov_risk)
      return(dt_res)
    }, s = s, t = t, lag = lag, Hs = Hs, Ls = Ls, Ht = Ht, Lt = Lt,
    bw_grid_s = bw_grid_s, bw_grid_t = bw_grid_t, dt_sigma = dt_sigma,
    kernel_smooth = smooth_ker, N = N, data = data,
    dt_EX2 = dt_EX2, dt_dependence_coef = dt_dependence_coef))

  return(dt_autocov_risk)
}

#' Estimate lag-\eqn{\ell} (\eqn{\ell > 0}) autocovariance function
#'
#' @inheritParams estimate_autocov_risk
#' @param optbw \code{vector (numeric)}. The optimal bandwidth parameter for lag-\eqn{\ell} (\eqn{\ell > 0}) autocovariance function estimation for each pair (\code{s}, \code{t}).
#' Default \code{optbw = NULL} and thus it will be estimated using \link{estimate_autocov_risk()} function.
#' @param center \code{logical (TRUE or FALSE)}. Default \code{center = TRUE} and so the curves are centred when the autocovariance is estimated: \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))}.
#' Otherwise, the two parts \eqn{\mathbb{E}X_0(s)X_{\ell}(t)} and \eqn{\mu(s)\mu(t)} will be estimated separately.
#' The first part with a bandwidth obtained with \link{estimate_autocov_risk} and the second part with a bandwidth obtained with \link{estimate_mean_risk}.
#' @param mean_estimates_s \code{vector (numeric)}. Mean function estimates at each \code{s}.
#' The default is \code{mean_estimates_s = NULL} and so it is estimated if \code{center = FALSE}.
#' @param mean_estimates_t \code{vector (numeric)}. Mean function estimates at each \code{t} if \code{center = FALSE}.
#' The default is \code{mean_estimates_t = NULL} and so it is estimated.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{s :}{ The first argument of the autocovariance function.}
#'            \item{t :}{ The second argument of the autocovariance function.}
#'            \item{locreg_bw :}{ The bandwidth used to estimate the local regularity parameters.}
#'            \item{Hs :}{ The estimates of the local exponent for each \code{s}.}
#'            \item{Ls :}{ The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
#'            \item{Ht :}{ The estimates of the local exponent for each \code{t}.}
#'            \item{Lt :}{ The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
#'            \item{lag :}{ The lag of the autocovariance. It corresponds to \eqn{\ell > 0}.}
#'            \item{optbw_gamma :}{ The optimal bandwidth for \eqn{\gamma_\ell(s,t)}, \eqn{\ell > 0}.}
#'            \item{optbw_muhat_s :}{ The optimal bandwidth for \eqn{\mu(s)}.}
#'            \item{optbw_muhat_t :}{ The optimal bandwidth for \eqn{\mu(t)}.}
#'            \item{muhat_s :}{ The estimates of the mean function \eqn{\mu(s)} for each \code{s}.}
#'            \item{muhat_t :}{ The estimates of the mean function \eqn{\mu(t)} for each \code{t}.}
#'            \item{gammahat :}{ The estimates of the \eqn{\gamma_\ell(s,t)} function for each (\code{s}, \code{t}).}
#'            \item{autocovhat :}{ The estimates of the lag-\eqn{\ell} autocovariance function for each (\code{s}, \code{t}).}
#'         }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_autocov_risk()].
#' @export
#'
#' @import data.table
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Coming ...
#'
#' # Coming ...
#'
#' # Coming ...
#' }
estimate_autocov_2bw <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                 s = c(1/5, 2/5, 4/5),
                                 t = c(1/4, 1/2, 3/4),
                                 lag = 1,
                                 optbw_s = NULL,
                                 optbw_t = NULL,
                                 bw_grid = seq(0.005, 0.15, len = 45),
                                 Hs = NULL, Ls = NULL,
                                 Ht = NULL, Lt = NULL,
                                 Delta = NULL, h = NULL,
                                 center = TRUE,
                                 mean_estimates_s = NULL,
                                 mean_estimates_t = NULL,
                                 smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  # Control on local regularity parameters
  if (((!is.null(Hs)) & length(Hs) != length(s)) | ((!is.null(Ls)) & length(Ls) != length(s)))
    stop("If 'Hs' or 'Ls' is not NULL, it must be the same length as 's'.")
  if ( ((!is.null(Ht)) & length(Ht) != length(t)) | ((!is.null(Lt)) & length(Lt) != length(t)))
    stop("If 'Ht' or 'Lt' is not NULL, it must be the same length as 't'.")
  if ((!is.null(Hs) & ! methods::is(Hs, "numeric")) |
      (!is.null(Ls) & ! methods::is(Ls, "numeric")) |
      (!is.null(Ht) & ! methods::is(Ht, "numeric")) |
      (!is.null(Lt) & ! methods::is(Lt, "numeric")))
    stop("If 'Hs', 'Ls', 'Ht' or 'Lt' is not NULL, it must be numeric.")

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Control the given optimal bandwidth and bandwidth grid
  if (((!is.null(optbw_s)) & length(optbw_s) != length(s)) | (!is.null(optbw_t)) & length(optbw_t) != length(t)) {
    stop("If 'optbw_s' or 'optbw_t' is not NULL, it must be the same length as 's' and 't'.")
  } else if (is.null(optbw_s) | is.null(optbw_s)) {
    if ((! is.null(bw_grid)) ) {
      if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
        stop("If 'bw_grid' is not NULL, it must be a vector of positive values between 0 and 1.")
    } else {
      lambdahat <- mean(data[, .N, by = "id_curve"][, N])
      K <- 20
      b0 <- 4 * max((N * lambdahat) ** (- 0.9), (N * (lambdahat ** 2)) ** (- 0.9))
      bK <- 4 * max((N * lambdahat) ** (- 1 / 3), (N * (lambdahat ** 2)) ** (- 1 / 3))
      a <- exp((log(bK) - log(b0)) / K)
      bw_grid <- b0 * a ** (seq_len(K))
      rm(K, b0, bK, a) ; gc()
    }
  }

  # If Non centered version of the autocovariance function
  if (! center) {
    if (((!is.null(mean_estimates_s)) & length(mean_estimates_s) != length(s)) |
        ((!is.null(mean_estimates_t)) & length(mean_estimates_t) != length(t)))
      stop("If 'mean_estimates_s' or 'mean_estimates_t' is not NULL, it must be the same length as 's' and 't'.")
    if ((!is.null(mean_estimates_s) & ! methods::is(mean_estimates_s, "numeric")) |
        (!is.null(mean_estimates_t) & ! methods::is(mean_estimates_t, "numeric")))
      stop("If 'mean_estimates_s', 'mean_estimates_t' is not NULL, it must be numeric.")
  }

  # Sort by s and t
  # dt_st <- data.table::data.table("s" = s, "t" = t)
  # dt_st <- dt_st[order(s,t)]
  # s <- dt_st[, s]
  # t <- dt_st[, t]
  # rm(dt_st) ; gc()

  # Estimate mean function at s and at t
  if ((! center) & (is.null(mean_estimates_s) | is.null(mean_estimates_t))) {
    ## Estimate mean function at s
    dt_mean_s <- estimate_mean(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = s, optbw = optbw_s, bw_grid = bw_grid,
      Ht = Hs, Lt = Ls,
      Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    data.table::setnames(x = dt_mean_s,
                         old = c("t", "locreg_bw", "Ht", "Lt", "optbw_s", "muhat"),
                         new = c("s", "locreg_bw_s", "Hs", "Ls", "optbw_muhat_s", "muhat_s"))

    ## Estimate mean function at s
    dt_mean_t <- estimate_mean(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = t, optbw = optbw_t, bw_grid = bw_grid,
      Ht = Ht, Lt = Lt,
      Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    data.table::setnames(x = dt_mean_t,
                         old = c("t", "locreg_bw", "Ht", "Lt", "optbw_s", "muhat"),
                         new = c("t", "locreg_bw_t", "Ht", "Lt", "optbw_muhat_t", "muhat_t"))
    ## Merge
    dt_mean <- cbind(dt_mean_s[, .SD, .SDcols = c("s", "locreg_bw_s", "Hs", "Ls", "optbw_muhat_s", "muhat_s")],
                     dt_mean_t[, .SD, .SDcols = c("t", "locreg_bw_t", "Ht", "Lt", "optbw_muhat_t", "muhat_t")])
    dt_mean[, locreg_bw := (locreg_bw_s + locreg_bw_t) / 2]
    dt_mean[, c("locreg_bw_s", "locreg_bw_t") := NULL]
    data.table::setcolorder(x = dt_mean,
                            neworder = c("s", "t", "locreg_bw", "Hs", "Ls", "Ht", "Lt",
                                         "optbw_muhat_s", "optbw_muhat_t", "muhat_s", "muhat_t"))

    ## Clean
    dt_mean <- dt_mean[order(s,t)]
    rm(dt_mean_s, dt_mean_t) ; gc()

    ## Set local regularity parameters estimates
    Hs <- dt_mean[, Hs]
    Ls <- dt_mean[, Ls]
    Ht <- dt_mean[, Ht]
    Lt <- dt_mean[, Lt]
    h <- dt_mean[, unique(locreg_bw)]
  } else if (! center) {
    dt_mean <- data.table::data.table(
      "s" = s, "t" = t, "locreg_bw_s" = NA, "locreg_bw_t" = NA, "Hs" = NA, "Ls" = NA, "Ht" = NA, "Lt" = NA,
      "optbw_muhat_s" = NA, "optbw_muhat_t" = NA,
      "muhat_s" = mean_estimates_s, "muhat_t" = mean_estimates_t)
  }

  # If the optimal bandwidth is not given
  if (is.null(optbw_s) | is.null(optbw_t)) {
    # Estimate the risk function
    dt_autocov_risk <- estimate_autocov_risk_2bw(
      data = data, idcol = idcol, tcol = tcol, ycol = ycol,
      s = s, t = t, lag = lag,
      bw_grid = bw_grid, smooth_ker = smooth_ker,
      Hs = Hs, Ls = Ls,
      Ht = Ht, Lt = Lt,
      Delta = Delta, h = h)

    # Take the optimum of the risk function
    dt_autocov_risk[, c("optbw_s", "optbw_t") := .(hs[which.min(autocov_risk)], ht[which.min(autocov_risk)]), by = c("s", "t")]
    dt_autocov_optbw <- unique(dt_autocov_risk[, list(s, t, optbw_s, optbw_t)])
    dt_autocov_optbw <- dt_autocov_optbw[order(s,t)]
  } else {
    # Set optimal autocovariance bandwidth function
    dt_autocov_optbw <- data.table::data.table("s" = s, "t" = t, "optbw_s" = optbw_s, "optbw_t" = optbw_t)
  }

  # Smooth curves with optimal bandwidth parameters
  dt_Xhat <- data.table::rbindlist(lapply(1:N, function(curve_index, s, t, data, dt_autocov_optbw, smooth_ker){
    # Extract data
    Tn <- data[id_curve == curve_index, tobs]
    Yn <- data[id_curve == curve_index, X]

    # \pi_n(s,h)
    pin_s <- sapply(X = 1:length(s), function(sidx, Tn, s, optbw_autocov){
      as.numeric(abs(Tn - s[sidx]) <= optbw_autocov[sidx])
    }, Tn = Tn, s, optbw_autocov = dt_autocov_optbw[order(s, t), optbw_s])
    pin_s <- t(pin_s)
    pin_s <- as.numeric(rowSums(pin_s, na.rm = TRUE) >= 1)

    # \pi_{n + \ell}(t,h)
    pin_t <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw_autocov){
      as.numeric(abs(Tn - t[tidx]) <= optbw_autocov[tidx])
    }, Tn = Tn, t, optbw_autocov = dt_autocov_optbw[order(s, t), optbw_t])
    pin_t <- t(pin_t)
    pin_t <- as.numeric(rowSums(pin_t, na.rm = TRUE) >= 1)

    # \widehat X(s;h)
    Xhat_s <- mapply(function(s, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = s, h = h, smooth_ker = ker)$yhat
    }, s = dt_autocov_optbw[, s], h = dt_autocov_optbw[, optbw_s],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    # \widehat X(t;h)
    Xhat_t <- mapply(function(t, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
    }, t = dt_autocov_optbw[, t], h = dt_autocov_optbw[, optbw_t],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t,
                                     "pin_s" = pin_s, "pin_t" = pin_t,
                                     "Xhat_s" = Xhat_s, "Xhat_t" = Xhat_t)
    return(dt_res)
  }, data = data, s = s, t = t, dt_autocov_optbw = dt_autocov_optbw, smooth_ker = smooth_ker))

  # Estimate \gamma_\ell(s,t, hs_opt, ht_opt)
  dt_gamma <- data.table::rbindlist(lapply(1:length(s), function(sid, s, t, lag, dt_Xhat, center){
    N <- dt_Xhat[s == s[sid] & t == t[sid], .N]

    ## For the argument s
    Xhat_s_vec <- dt_Xhat[s == s[sid] & t == t[sid], Xhat_s]
    pin_s_vec <- dt_Xhat[s == s[sid] & t == t[sid], pin_s]
    Xhat_s_i <- Xhat_s_vec[1:(N - lag)]
    Xhat_s_i[is.nan(Xhat_s_i)] <- 0
    pin_s_i <- pin_s_vec[1:(N - lag)]

    ## For the argument t
    Xhat_t_vec <- dt_Xhat[s == s[sid] & t == t[sid], Xhat_t]
    pin_t_vec <- dt_Xhat[s == s[sid] & t == t[sid], pin_t]
    Xhat_t_i_plus_lag <- Xhat_t_vec[(1 + lag):N]
    Xhat_t_i_plus_lag[is.nan(Xhat_t_i_plus_lag)] <- 0
    pin_t_i_plus_lag <- pin_t_vec[(1 + lag):N]

    ## Number of used points for autocovariance
    PNl <- sum(pin_s_i *pin_t_i_plus_lag)

    ## Estimate gamma
    XsXt_vec <- pin_s_i * pin_t_i_plus_lag * Xhat_s_i * Xhat_t_i_plus_lag
    gammahat <- sum(XsXt_vec) / PNl

    if (center) {
      ## Estimate mean fun
      PNs <- sum(pin_s_vec)
      PNt <- sum(pin_t_vec)
      Xhat_s_vec[is.nan(Xhat_s_vec)] <- 0
      Xhat_t_vec[is.nan(Xhat_t_vec)] <- 0
      muhat_s <- sum(Xhat_s_vec * pin_s_vec) / PNs
      muhat_t <- sum(Xhat_t_vec * pin_t_vec) / PNt
      XsXt_vec_centered <- pin_s_i * pin_t_i_plus_lag * (Xhat_s_i - muhat_s) * (Xhat_t_i_plus_lag - muhat_t)
      autocovhat <- sum(XsXt_vec_centered) / PNl
      dt_res <- data.table::data.table("s" = s[sid], "t" = t[sid], "lag" = lag, "PNl" = PNl,
                                       "gammahat" = gammahat, "autocovhat" = autocovhat,
                                       "muhat_s" = muhat_s, "muhat_t" = muhat_t)
    } else {
      dt_res <- data.table::data.table("s" = s[sid], "t" = t[sid], "lag" = lag, "PNl" = PNl, "gammahat" = gammahat)
    }
  }, s = s, t = t, lag = lag, dt_Xhat = dt_Xhat, center = center))

  ## Centered version
  if (center) {
    # Estimate  lag-\ell autocovariance
    # dt_gamma[, autocovhat := gammahat - muhat_s * muhat_t]
    dt_autocov <- data.table::merge.data.table(
      x = dt_gamma, y = dt_autocov_optbw[, list("optbw_gamma_s" = optbw_s, "optbw_gamma_t" = optbw_t), by = c("s", "t")],
      by = c("s", "t")
    )
    dt_autocov[, c("locreg_bw_s", "locreg_bw_t", "Hs", "Ls", "Ht", "Lt", "lag", "PNl",
                   "optbw_muhat_s", "optbw_muhat_t") := .(NA, NA, NA, NA, NA, NA, lag, NA, NA, NA)]
    data.table::setcolorder(
      x = dt_autocov,
      neworder = c("s", "t", "locreg_bw_s", "locreg_bw_t", "Hs", "Ls", "Ht", "Lt", "lag", "PNl",
                   "optbw_gamma_s", "optbw_gamma_t", "optbw_muhat_s", "optbw_muhat_t",
                   "muhat_s", "muhat_t", "gammahat", "autocovhat"))
  } else {

    # Estimate  lag-\ell autocovariance
    dt_autocov <- data.table::merge.data.table(
      x = dt_mean, y = dt_gamma, by = c("s", "t")
    )
    dt_autocov <- data.table::merge.data.table(
      x = dt_autocov, y = dt_autocov_optbw[, list("optbw_gamma_s" = optbw_s, "optbw_gamma_t" = optbw_t), by = c("s", "t")],
      by = c("s", "t")
    )
    dt_autocov[, autocovhat := gammahat - muhat_s * muhat_t, by = c("s", "t")]
    data.table::setcolorder(
      x = dt_autocov,
      neworder = c("s", "t", "locreg_bw_s", "locreg_bw_t", "Hs", "Ls", "Ht", "Lt", "lag", "PNl",
                   "optbw_gamma_s", "optbw_gamma_t", "optbw_muhat_s", "optbw_muhat_t",
                   "muhat_s", "muhat_t", "gammahat", "autocovhat"))
  }
  return(dt_autocov)
}


