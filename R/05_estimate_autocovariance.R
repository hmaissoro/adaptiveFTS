#' Estimate the risk of the lag-\eqn{\ell} (\eqn{\ell > 0}) autocovariance function
#'
#' @inheritParams .format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each pair (\code{s}, \code{t}).
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
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{Hs :}{ The estimates of the local exponent for each \code{s}.}
#'            \item{Ls :}{ The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
#'            \item{Ht :}{ The estimates of the local exponent for each \code{t}.}
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
estimate_autocov_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                  s = c(1/5, 2/5, 4/5),
                                  t = c(1/4, 1/2, 3/4),
                                  lag = 1,
                                  bw_grid = seq(0.005, 0.15, len = 45),
                                  smooth_ker = epanechnikov,
                                  Hs = NULL, Ls = NULL,
                                  Ht = NULL, Lt = NULL,
                                  Delta = NULL, h = NULL){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
    stop("'bw_grid' must be a vector of positive values between 0 and 1.")

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
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st)
  gc()

  # Estimate local regularity parameters
  # This function controls the remaining arguments
  if (is.null(Hs) | is.null(Ls)) {
    dt_locreg_s <- estimate_locreg(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = s, Delta = Delta, h = h,
      smooth_ker = smooth_ker)
    Hs <- dt_locreg_s[, H]
    Ls <- dt_locreg_s[, L]
    hs <- dt_locreg_s[, unique(h)]
  } else {
    hs <- h
  }
  if (is.null(Ht) | is.null(Lt)) {
    dt_locreg_t <- estimate_locreg(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t, Delta = Delta, h = h,
    smooth_ker = smooth_ker)
    Ht <- dt_locreg_t[, H]
    Lt <- dt_locreg_t[, L]
    ht <- dt_locreg_t[, unique(h)]
  } else {
    ht <- h
  }
  rm(dt_locreg_s, dt_locreg_t)
  # Estimation of the observation error standard deviation
  dt_sigma_s <- estimate_sigma(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = s)
  dt_sigma_t <- estimate_sigma(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t)
  dt_sigma <- cbind(
    dt_sigma_s[, list("s" = t, "sig_error_s" = sig)],
    dt_sigma_t[, list("t" = t, "sig_error_t" = sig)]
  )
  rm(dt_sigma_s, dt_sigma_t)
  gc()

  # Estimation of the second order moment using the presmoothing bandwidth
  ## Smooth curves at s and at t
  dt_Xhat <- data.table::rbindlist(
    lapply(1:N, function(curve_index, s, t, presmooth_bw_s, presmooth_bw_t, kernel_smooth, data){
      Xhat_s <- estimate_nw(
        y = data[id_curve %in% curve_index][order(tobs), X],
        t = data[id_curve %in% curve_index][order(tobs), tobs],
        tnew = s,
        h = presmooth_bw_s,
        smooth_ker = kernel_smooth)
      Xhat_s <- Xhat_s[, yhat]

      Xhat_t <- estimate_nw(
        y = data[id_curve %in% curve_index][order(tobs), X],
        t = data[id_curve %in% curve_index][order(tobs), tobs],
        tnew = t,
        h = presmooth_bw_t,
        smooth_ker = kernel_smooth)
      Xhat_t <- Xhat_t[, yhat]

      dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t, "Xhat_s" = Xhat_s, "Xhat_t" = Xhat_t)
    }, s = s, t = t,
    presmooth_bw_s = hs, presmooth_bw_t = ht,
    kernel_smooth = smooth_ker, data = data))

  ##  Transform NaN to NA
  dt_Xhat[is.nan(Xhat_s), Xhat_s := NA]
  dt_Xhat[is.nan(Xhat_t), Xhat_t := NA]
  dt_EX2 <- dt_Xhat[, list("EX2_s" = mean(Xhat_s, na.rm = TRUE), "EX2_t" = mean(Xhat_t, na.rm = TRUE)),
                        by = c("s", "t")]
  rm(dt_Xhat)
  gc()

  # Estimation of the empirical autocovariance of X_0(s)X_{\ell}(s)
  ## The bandwidth is taking as the mean of the presmoothing bandwidths at s and t
  if (is.null(hs) | is.null(ht)) {
    h_XsXt <- NULL
  } else {
    h_XsXt <- (hs + ht) / 2
  }
  dt_XsXt_autocov <- estimate_empirical_XsXt_autocov(
    data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = s, t = t, cross_lag = lag,
    lag = 0:(N - lag - 1),
    h = h_XsXt,
    smooth_ker = smooth_ker)
  rm(hs, ht, h_XsXt)

  # Estimate the risk function
  dt_autocov_risk <- data.table::rbindlist(
    lapply(bw_grid, function(h, s, t, lag, Hs, Ls, Ht, Lt,
                             kernel_smooth, dt_sigma, N, data, dt_EX2, dt_XsXt_autocov){
      # compute quantities to estimate estimate the risk
      dt_risk <- data.table::rbindlist(lapply(1:N, function(curve_index, h, s, t, Hs, Ls, Ht, Lt, kernel_smooth, data){
        # Extract the observation points per curve
        # and compute the weight of the estimator
        ## By convention NaN = 0
        Tn <- data[id_curve == curve_index, tobs]

        ## For the argument s
        ker_s <- outer(X = s, Y = Tn, FUN = function(si, Tnm, Ker) Ker((Tnm - si) / h), Ker = kernel_smooth)
        wker_s <- apply(ker_s, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
        wker_s <- t(wker_s)

        ### Take the maximum and c_n(s,h), the sum of the weight
        wmax_s <- apply(wker_s, 1, max)
        cn_s <- apply(wker_s, 1, function(r) sum(abs(r)))

        Tn_s_2H_s <- outer(
          X = 1:length(s), Y = Tn,
          FUN = function(sid, Tnm, s, Hs) abs((Tnm - s[sid]) / h) ** (2 * Hs[sid]),
          s = s, Hs = Hs
        )

        bn2H_s <- diag(abs(wker_s) %*% t(Tn_s_2H_s))

        ### \pi_n(s;h)
        pin_s <- outer(X = s, Y = Tn, FUN = function(si, Tnm, h) as.numeric(abs(Tnm - si) <= h), h = h)
        pin_s <- as.numeric(rowSums(pin_s) >= 1)

        ## For argument t
        ### By convention NaN = 0
        ker_t <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, Ker) Ker((Tnm - ti) / h), Ker = kernel_smooth)
        wker_t <- apply(ker_t, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
        wker_t <- t(wker_t)

        ### Take the maximum and c_n(t,h), the sum of the weight
        wmax_t <- apply(wker_t, 1, max)
        cn_t <- apply(wker_t, 1, function(r) sum(abs(r)))

        Tn_t_2H_t <- outer(
          X = 1:length(t), Y = Tn,
          FUN = function(tid, Tnm, t, Ht) abs((Tnm - t[tid]) / h) ** (2 * Ht[tid]),
          t = t, Ht = Ht
        )

        bn2H_t <- diag(abs(wker_t) %*% t(Tn_t_2H_t))

        # \pi_n(t;h)
        pin_t <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, h) as.numeric(abs(Tnm - ti) <= h), h = h)
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
      }, h = h, s = s, t = t, Hs = Hs, Ht = Ht, kernel_smooth = kernel_smooth, data = data))

      # Split and merge by lag
      ## The argument s is associated to the curves n = 1,..., N-lag
      dt_risk_s <- dt_risk[, list(id_curve, s, t, wmax_s, cn_s, bn2H_s, pin_s)]
      dt_risk_s <- dt_risk_s[id_curve %in% 1:(N - lag)]
      dt_id_lag_s <- data.table::data.table(
        "id_curve" = sort(unique(dt_risk_s[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_risk_s <- data.table::merge.data.table(
        x = dt_id_lag_s,
        y = dt_risk_s,
        by = "id_curve")

      ## The argument t is associated to the curves n = 1 + lag,..., N
      dt_risk_t <- dt_risk[, list(id_curve, s, t, wmax_t, cn_t, bn2H_t, pin_t)]
      dt_risk_t <- dt_risk_t[id_curve %in% (1 + lag):N]
      dt_id_lag_t <- data.table::data.table(
        "id_curve" = sort(unique(dt_risk_t[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_risk_t <- data.table::merge.data.table(
        x = dt_id_lag_t,
        y = dt_risk_t,
        by = "id_curve")

      ## Merge
      dt_risk_merge <- data.table::merge.data.table(
        x = dt_risk_s,
        y = dt_risk_t,
        by = c("id_lag", "s", "t"))
      dt_risk_merge <- dt_risk_merge[order(id_curve.x)]

      ## Clean unnecessary object
      rm(dt_id_lag_s, dt_id_lag_t, dt_risk_s, dt_risk_t, dt_risk)
      gc()

      # Compute P_N(s, t; h)
      dt_risk_merge[, PNl := sum(pin_s * pin_t), by = c("s", "t")]

      # compute \mathbb B(s|t,h, 2H_s, 0) and \mathbb B(t|s,h, 2H_t, \ell)
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
      bias_term <- 3 * dt_risk[, EX2_t] * Ls * (h ** (2 * Hs)) * dt_risk[, Bs] +
        3 * dt_risk[, EX2_s] * Lt * (h ** (2 * Ht)) * dt_risk[, Bt]

      ## Variance term
      varriance_term <- 3 * (dt_risk[, sig_error_s] ** 2) * dt_risk[, EX2_t] * dt_risk[, Vgamma1] +
        3 * (dt_risk[, sig_error_t] ** 2) * dt_risk[, EX2_s] * dt_risk[, Vgamma2] +
        3 * (dt_risk[, sig_error_s] ** 2) * (dt_risk[, sig_error_t] ** 2) * dt_risk[, Vgamma]

      ## Dependence term
      ### Compute p_k(s,t;h)
      dt_p_k <- unique(dt_risk_merge[, list(id_curve.x, id_curve.y, id_lag, s, t, pin_s, pin_t, PNl)])
      dt_p <- data.table::rbindlist(lapply(1:(N-lag - 1), function(k, s, t, dt_p_k){
        data.table::rbindlist(lapply(1:length(s), function(sid, k, s, t, dt_p_k){
          PNl <- unique(dt_p_k[s == s[sid] & t == t[sid], PNl])

          # Compute \pi_i(s; h) and \pi_{i + k}(t; h)
          pin_s_vec <- dt_p_k[s == s[sid] & t == t[sid]][order(id_curve.x), pin_s]
          pin_s_i <- pin_s_vec[1:((N - lag) - k )]
          pin_s_i_plus_k <- pin_s_vec[(1 + k):(N - lag)]

          # Compute \pi_{i + \ell}}(t; h) and \pi_{i + \ell + k}(t;h)
          pin_t_vec <- dt_p_k[s == s[sid] & t == t[sid]][order(id_curve.y), pin_t]
          pin_t_i_plus_ell <- pin_t_vec[1:((N - lag) - k)]
          pin_t_i_plus_ell_plus_k <- pin_t_vec[(1 + k):(N - lag)]

          p_k <- sum((pin_s_i * pin_t_i_plus_ell) * (pin_s_i_plus_k * pin_t_i_plus_ell_plus_k)) / PNl
          dt_res <- data.table::data.table("s" = s[sid], "t" = t[sid], "lag" = k, "pk" = p_k)

          rm(pin_s_vec, pin_s_i, pin_s_i_plus_k, pin_t_vec, pin_t_i_plus_ell, pin_t_i_plus_ell_plus_k, PNl)
          gc()
          return(dt_res)
        }, k = k, s = s, t = t, dt_p_k = dt_p_k))
      }, s = s, t = t, dt_p_k = dt_p_k))

      dt_long_run_var <- data.table::merge.data.table(
        x = dt_XsXt_autocov[lag > 0],
        y = dt_p,
        by = c("s", "t", "lag")
      )
      dt_long_run_var <- dt_long_run_var[, list("long_run_var" = sum(2 * XsXt_autocov * pk)), by = c("s", "t")]
      dependence_coef <- dt_XsXt_autocov[lag == 0][order(s, t), XsXt_autocov] + dt_long_run_var[order(s, t), long_run_var]
      ### Note that dependence_coef <= abs(dependence_coef), thus
      dependence_coef <- abs(dependence_coef)
      dependence_term <- dependence_coef  / dt_risk[order(s, t), PNl]

      ## Final risk function
      autocov_risk <- 2 * (bias_term + varriance_term + dependence_term)

      # Result to returned
      dt_res <- data.table::data.table("s" = s, "t" = t, "h" = h, "Hs" = Hs, "Ls" = Ls, "Ht" = Ht, "Lt" = Lt,
                                       "bias_term" = bias_term, "varriance_term" = varriance_term,
                                       "dependence_term" = dependence_term, "autocov_risk" = autocov_risk)
      return(dt_res)
    }, s = s, t = t, lag = lag, Hs = Hs, Ls = Ls, Ht = Ht, Lt = Lt,
    dt_sigma = dt_sigma, kernel_smooth = smooth_ker, N = N, data = data,
    dt_EX2 = dt_EX2, dt_XsXt_autocov = dt_XsXt_autocov))

  return(dt_autocov_risk)
}

#' Estimate lag-\eqn{\ell} (\eqn{\ell > 0}) autocovariance function
#'
#' @inheritParams estimate_autocov_risk
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
estimate_autocov <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             s = c(1/5, 2/5, 4/5),
                             t = c(1/4, 1/2, 3/4),
                             lag = 1,
                             bw_grid = seq(0.005, 0.15, len = 45),
                             Delta = NULL, h = NULL, smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'s' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! (methods::is(s, "numeric") & all(data.table::between(s, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! length(s) == length(t))
    stop("Arguments 's' and 't' must be of equal length.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! (all(methods::is(bw_grid, "numeric") & data.table::between(bw_grid, 0, 1)) & length(bw_grid) > 1))
    stop("'bw_grid' must be a vector of positive values between 0 and 1.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(lag < 0)| (length(lag) > 1) | any(lag - floor(lag) > 0) | any(N <= lag))
    stop("'lag' must be a positive integer lower than the number of curves.")

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st)
  gc()

  # Estimate mean function at s and at t
  dt_mean_s <- estimate_mean(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol,
    t = s, bw_grid = bw_grid, Delta = Delta,
    h = h, smooth_ker = smooth_ker)
  data.table::setnames(x = dt_mean_s,
                       old = c("t", "locreg_bw", "H", "L", "optbw", "muhat"),
                       new = c("s", "locreg_bw_s", "Hs", "Ls", "optbw_muhat_s", "muhat_s"))
  dt_mean_t <- estimate_mean(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol,
    t = t, bw_grid = bw_grid, Delta = Delta,
    h = h, smooth_ker = smooth_ker)
  data.table::setnames(x = dt_mean_t,
                       old = c("t", "locreg_bw", "H", "L", "optbw", "muhat"),
                       new = c("t", "locreg_bw_t", "Ht", "Lt", "optbw_muhat_t", "muhat_t"))
  dt_mean <- cbind(dt_mean_s[, .SD, .SDcols = c("s", "locreg_bw_s", "Hs", "Ls", "optbw_muhat_s", "muhat_s")],
                   dt_mean_t[, .SD, .SDcols = c("t", "locreg_bw_t", "Ht", "Lt", "optbw_muhat_t", "muhat_t")])
  dt_mean[, locreg_bw := (locreg_bw_s + locreg_bw_t) / 2]
  dt_mean[, c("locreg_bw_s", "locreg_bw_t") := NULL]
  data.table::setcolorder(x = dt_mean,
                          neworder = c("s", "t", "locreg_bw", "Hs", "Ls", "Ht", "Lt",
                                       "optbw_muhat_s", "optbw_muhat_t", "muhat_s", "muhat_t"))
  dt_mean <- dt_mean[order(s,t)]
  rm(dt_mean_s, dt_mean_t) ; gc()

  # Estimate the risk function
  dt_autocov_risk <- estimate_autocov_risk(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol,
    s = s, t = t, lag = lag,
    bw_grid = bw_grid, smooth_ker = smooth_ker,
    Hs = dt_mean[, Hs], Ls = dt_mean[, H],
    Ht = dt_mean[, Ht], Lt = dt_mean[, Lt],
    Delta = NULL, h = dt_mean[, unique(locreg_bw)]
  )

  # Take the optimum of the risk function
  dt_autocov_risk[, optbw := h[which.min(autocov_risk)], by = c("s", "t")]
  dt_autocov_optbw <- unique(dt_autocov_risk[, list(s, t, optbw)])
  dt_autocov_optbw <- dt_autocov_optbw[order(s,t)]

  # Smooth curves with optimal bandwidth parameters
  dt_Xhat <- data.table::rbindlist(lapply(1:N, function(curve_index, s, t, data, dt_autocov_optbw, smooth_ker){
    # Extract data
    Tn <- data[id_curve == curve_index, tobs]
    Yn <- data[id_curve == curve_index, X]

    # \pi_n(s,h)
    pin_s <- sapply(X = s, function(si, Tn, dt_autocov_optbw){
      as.numeric(abs(Tn - si) <= dt_autocov_optbw[s == si, optbw])
    }, Tn = Tn, dt_autocov_optbw = dt_autocov_optbw)
    pin_s <- t(pin_s)
    pin_s <- as.numeric(rowSums(pin_s, na.rm = TRUE) >= 1)

    # \pi_{n + \ell}(t,h)
    pin_t <- sapply(X = t, function(ti, Tn, dt_autocov_optbw){
      as.numeric(abs(Tn - ti) <= dt_autocov_optbw[t == ti, optbw])
    }, Tn = Tn, dt_autocov_optbw = dt_autocov_optbw)
    pin_t <- t(pin_t)
    pin_t <- as.numeric(rowSums(pin_t, na.rm = TRUE) >= 1)

    # \widehat X(s;h)
    Xhat_s <- mapply(function(s, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = s, h = h, smooth_ker = ker)$yhat
    }, s = dt_autocov_optbw[, s], h = dt_autocov_optbw[, optbw],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    # \widehat X(t;h)
    Xhat_t <- mapply(function(t, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
    }, t = dt_autocov_optbw[, t], h = dt_autocov_optbw[, optbw],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t,
                                     "pin_s" = pin_s, "pin_t" = pin_t,
                                     "Xhat_s" = Xhat_s, "Xhat_t" = Xhat_t)
    return(dt_res)
  }, data = data, s = s, t = t, dt_autocov_optbw = dt_autocov_optbw, smooth_ker = smooth_ker))

  # Estimate \gamma_\ell(s,t, h_opt)
  dt_gamma <- data.table::rbindlist(lapply(1:length(s), function(sid, s, t, lag, dt_Xhat){
    N <- dt_Xhat[s == s[sid] & t == t[sid], .N]

    ## For the argument s
    Xhat_s_vec <- dt_Xhat[s == s[sid] & t == t[sid], Xhat_s]
    pin_s_vec <- dt_Xhat[s == s[sid] & t == t[sid], pin_s]
    Xhat_s_i <- Xhat_s_vec[1:(N - lag)]
    pin_s_i <- pin_s_vec[1:(N - lag)]

    ## For the argument t
    Xhat_t_vec <- dt_Xhat[s == s[sid] & t == t[sid], Xhat_t]
    pin_t_vec <- dt_Xhat[s == s[sid] & t == t[sid], pin_t]
    Xhat_t_i_plus_lag <- Xhat_t_vec[(1 + lag):N]
    pin_t_i_plus_lag <- pin_t_vec[(1 + lag):N]

    ## Estimate gamma
    XsXt_vec <- pin_s_i * pin_t_i_plus_lag * Xhat_s_i * Xhat_t_i_plus_lag

    XsXt_vec <- XsXt_vec[!is.nan(XsXt_vec)]
    gammahat <- mean(XsXt_vec)
    dt_res <- data.table::data.table("s" = s[sid], "t" = t[sid], "lag" = lag, "gammahat" = gammahat)
  }, s = s, t = t, lag = lag, dt_Xhat = dt_Xhat))

  # Estimate  lag-\ell autocovariance
  dt_autocov <- data.table::merge.data.table(
    x = dt_mean, y = dt_gamma, by = c("s", "t")
  )
  dt_autocov <- data.table::merge.data.table(
    x = dt_autocov, y = dt_autocov_optbw[, list("optbw_gamma" = optbw), by = c("s", "t")],
    by = c("s", "t")
  )
  dt_autocov[, autocovhat := gammahat - muhat_s * muhat_t, by = c("s", "t")]
  data.table::setcolorder(
    x = dt_autocov,
    neworder = c("s", "t", "locreg_bw", "Hs", "Ls", "Ht", "Lt", "lag",
                 "optbw_gamma", "optbw_muhat_s", "optbw_muhat_t",
                 "muhat_s", "muhat_t", "gammahat", "autocovhat"))
  return(dt_autocov)
}

