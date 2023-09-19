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
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()].
#'
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist between
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

  # Estimation of the variance using the presmoothing bandwidth
  dt_autocov <- estimate_empirical_autocov(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t, lag = 0:(N-1),
    h = dt_locreg[, unique(h)],
    smooth_ker = smooth_ker)

  # Estimate the risk function
  dt_autocov_risk <- data.table::rbindlist(lapply(bw_grid, function(h, s, t, lag, Hs, Ls, Ht, Lt, kernel_smooth, sig_error_s, sig_error_t, N, data, dt_autocov){
    # compute quantities to estimate estimate the risk
    dt_risk <- data.table::rbindlist(lapply(1:N, function(curve_index, h, s, t, lag, Hs, Ls, Ht, Lt, kernel_smooth, data){
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

    # Compute P_N(t, h)
    dt_risk[, PN := sum(pi_n), by = t]

    # compute \mathbb B(t,h, 2H) for each n = 1, ..., N
    dt_risk[, B := (sum(pi_n * cn * bn2H) / PN) ** 2, by = t]

    # Compute \mathbb V_mu(t,h) for each n = 1, ..., N
    dt_risk[, Vmu := sum(pi_n * cn * wmax) / PN ** 2, by = t]

    # Compute the risk for each h
    dt_rk <- unique(dt_risk[, list(t, B, Vmu, PN)])

    ## Biais term
    bias_term <- 2 * L * (h ** (2 * H)) * dt_rk[, B]

    ## Variance term
    varriance_term <- 2 * (sig_error ** 2) * dt_rk[, Vmu]

    ## Convergence term to true mean
    ### Compute \rho_\ell(t;h)
    dt_rho_ell <- unique(dt_risk[, list(id_curve, t, pi_n, PN)])
    dt_rho <- data.table::rbindlist(lapply(1:(N-1), function(ell, t, dt_rho_ell){
      data.table::rbindlist(lapply(t, function(ti, ell, dt_rho_ell){
        PN <- unique(dt_rho_ell[t == ti, PN])
        pi_vec <- dt_rho_ell[t == ti][order(id_curve), pi_n]
        pi_i <- pi_vec[1:(N - ell)]
        pi_i_plus_ell <- pi_vec[(1 + ell):N]
        rho_ell <- sum(pi_i * pi_i_plus_ell) / PN
        dt_res <- data.table::data.table("t" = ti, "lag" = ell, "rho" = rho_ell)
        return(dt_res)
      }, ell = ell, dt_rho_ell = dt_rho_ell))
    }, t = t, dt_rho_ell =dt_rho_ell))

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
  Ht = dt_locreg_t[, H], Lt = dt_locreg_t[, L],
  sig_error_s = dt_sigma_s[, sig], sig_error_t = dt_sigma_t[, sig],
  kernel_smooth = smooth_ker, N = N, data = data, dt_autocov = dt_autocov))

  return(dt_mean_risk)
}
