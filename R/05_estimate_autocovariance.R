#' Estimate the risk of the mean function
#'
#' @inheritParams .format_data
#' @param s \code{vector (numeric)}. First argument of the autocovariance function.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param lag \code{vector (positive integer)}. Lag of the autocovariance.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each pair (\code{s}, \code{t}).
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{s} and \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The points at which the risk function is estimated.}
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{H :}{ The estimates of the local exponent.}
#'            \item{L :}{ The estimates of the Hölder constant. It corresponds to $L_t^2$.}
#'            \item{bias_term :}{ The bias term of the risk function.}
#'            \item{varriance_term :}{ The variance term of the risk function.}
#'            \item{dependence_term :}{ The dependence term of the risk function.}
#'            \item{mean_risk :}{ The estimates of the risk function of the mean.}
#'         }
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
#'
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist between setnames
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
#' dt_mean_risk <- estimate_mean_risk(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
#'   Delta = NULL, h = NULL, smooth_ker = epanechnikov)
#'
#' # Plot mean risk
#' dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
#'       main = "t = 0.25", xlab = "h", ylab = "risk function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
#'       main = "t = 0.5", xlab = "h", ylab = "risk function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
#'       main = "t = 0.75", xlab = "h", ylab = "risk function")
#'   ),
#'   nrow = 3
#' )
#'
#' }
#'
#'
estimate_autocov_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
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
  if (! (all(methods::is(bw_grid, "numeric") & data.table::between(t, 0, 1)) & length(bw_grid) > 1))
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

  # Estimate local regularity parameters
  # This function controls the remaining arguments
  dt_locreg_s <- estimate_locreg(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = s, Delta = Delta, h = h,
    smooth_ker = smooth_ker)
  dt_locreg_t <- estimate_locreg(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t, Delta = Delta, h = h,
    smooth_ker = smooth_ker)

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

  # Estimation of the secon order moment using the presmoothing bandwidth
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
    presmooth_bw_s = dt_locreg_s[, unique(h)],
    presmooth_bw_t = dt_locreg_t[, unique(h)],
    kernel_smooth = smooth_ker, data = data))

  ##  Transform NaN to NA
  dt_Xhat[is.nan(Xhat_s), Xhat_s := NA]
  dt_Xhat[is.nan(Xhat_t), Xhat_t := NA]
  dt_EX2 <- dt_Xhat[, list("EX2_s" = mean(Xhat_s, na.rm = TRUE), "EX2_t" = mean(Xhat_t, na.rm = TRUE)),
                        by = c("s", "t")]
  rm(dt_Xhat)
  gc()

  # Estimate the risk function
  dt_autocov_risk <- data.table::rbindlist(
    lapply(bw_grid, function(h, s, t, lag, Hs, Ls, Ht, Lt,
                             kernel_smooth, dt_sigma, N, data, dt_EX2){
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
      varriance_term <- 3 * dt_risk[, sig_error_s] * dt_risk[, EX2_t] * dt_risk[, Vgamma1] +
        3 * dt_risk[, sig_error_t] * dt_risk[, EX2_s] * dt_risk[, Vgamma2] +
        3 * dt_risk[, sig_error_s] * dt_risk[, sig_error_t] * dt_risk[, Vgamma]

      ## Convergence term to true mean
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

      dt_lr_var <- data.table::merge.data.table(
        x = dt_autocov[lag > 0],
        y = dt_rho,
        by = c("t", "lag")
      )
      dt_lr_var <- dt_lr_var[, list("lr_var" = sum(2 * autocov * rho)), by = t]
      dependence_coef <- dt_autocov[lag == 0][order(t), autocov] + dt_lr_var[order(t), lr_var]
      ### Note that dependence_coef <= abs(dependence_coef), thus
      dependence_coef <- abs(dependence_coef)
      dependence_term <- 2 * dependence_coef  / dt_rk[, PN]

      ## Final risk function
      mean_risk <- bias_term + varriance_term + dependence_term

      # Result to returned
      dt_res <- data.table::data.table("t" = t, "h" = h, "H" = H, "L" = L,
                                       "bias_term" = bias_term, "varriance_term" = varriance_term,
                                       "dependence_term" = dependence_term, "mean_risk" = mean_risk)
      return(dt_res)
    }, s = s, t = t, lag = lag, Hs = dt_locreg_s[, H], Ls = dt_locreg_s[, L],
    Ht = dt_locreg_t[, H], Lt = dt_locreg_t[, L], dt_sigma = dt_sigma,
    kernel_smooth = smooth_ker, N = N, data = data, dt_EX2 = dt_EX2))

  return(dt_mean_risk)
}
