#' Estimate the the standard deviation of the observation error
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the standard deviation of the error.
#'
#' @return A data.table with two columns: \code{t} and \code{sig} corresponding to the estimated standard deviation.
#' @export
#'
#' @import data.table
#'
estimate_sigma <- function(data, idcol = NULL, tcol = "tobs", ycol = "X", t = c(1/4, 1/2, 3/4)) {
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  dt_sig <- data.table::rbindlist(lapply(data[, unique(id_curve)], function(idx, data, t) {
    # Compute get the index the two T_{n,k} closest to time for each X_n
    ## Get the rank i_t and j_t
    data_idx <- data[id_curve == idx]
    dif_mat <- outer(X = t, Y = data_idx[, tobs], function(ti, tobsi) abs(ti - tobsi))
    rank_mat <- apply(X = dif_mat, MARGIN = 1, function(r) rank(r))

    ## Compute [Y_{n,i_t} - Y_{n,j_t}]^2 / 2 for each curve X_n and each t
    Zn <- apply(X = rank_mat, MARGIN = 2, function(c, dt){
      first <- which(c == 1)
      second <- which(c == 2)
      Y1 <- dt[first][, X]
      Y2 <- dt[second][, X]
      Zn <- 1 / 2 * (Y1 - Y2) ** 2
    }, dt = data_idx)
    dt_res <- data.table::data.table("t" = t, "Zn" = Zn)
    rm(data_idx, Zn)
    return(dt_res)
  }, data = data, t = t))
  # Estimate sigma
  dt_sig[, sig := sqrt(mean(Zn)), by = t]
  dt_sig <- unique(dt_sig[, list(t, sig)])

  return(dt_sig)
}

#' Estimate empirical autocovariance function
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the empirical autocovariance function.
#' @param lag \code{vector (integer)}. Lag of the autocovariance.
#' @param h \code{numeric (positive vector or scalar)}. The smoothing bandwidth parameter.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A data.table with three columns: \code{t}, \code{lag} and \code{autocov} corresponding to the estimated autocovariance.
#' @export
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist
#'
#' @seealso [estimate_nw()]
#'
estimate_empirical_autocov <- function(data, idcol = NULL, tcol = "tobs", ycol = "X",
                                       t = c(1/4, 1/2, 3/4), lag = c(0, 1, 2), h = NULL,
                                       smooth_ker = epanechnikov){
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  if (any(N <= lag))
    stop("'lag' must be lower than the number of curves.")
  if (! all(methods::is(t, "numeric") & data.table::between(t, 0, 1)))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

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
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    if (N > 50) {
      sample_curves <- sample(x = 1:N, size = 30)
    } else {
      sample_curves <- 1:N
    }
    h <- median(unlist(lapply(sample_curves, function(i, lambdahat, data, smooth_ker){
      K <- 100
      b0 <- 2 / lambdahat
      bK <- lambdahat ** (- 1 / 3)
      a <- exp((log(bK) - log(b0)) / K)
      hgrid <- b0 * a ** (seq_len(K))
      hbest <- estimate_nw_bw(y = data[id_curve == i, X],
                              t = data[id_curve ==i, tobs],
                              h_grid = hgrid,
                              smooth_ker = smooth_ker)

    }, lambdahat = lambdahat, data = data, smooth_ker = smooth_ker)))
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Smooth curves
  dt_smooth <- data.table::rbindlist(lapply(1:N, function(curve_index, hvec, t, kernel_smooth, data){
    x_smooth <- estimate_nw(y = data[id_curve == curve_index, X],
                            t = data[id_curve == curve_index, tobs],
                            tnew = t,
                            h = hvec[curve_index],
                            smooth_ker = kernel_smooth)
    dt_res <- data.table::data.table("curve_index" = curve_index, "t" = t, "x" = x_smooth$yhat)
    return(dt_res)
  }, hvec = h, t = t, kernel_smooth = smooth_ker, data = data))

  # Estimate mean function
  dt_mean <- dt_smooth[!is.nan(x), list("muhat" = mean(x)), by = "t"]

  # Estimate autocovariance
  dt_autocov <- data.table::rbindlist(lapply(lag, function(lg, dt_smooth, dt_mean){
    data.table::rbindlist(lapply(dt_smooth[, unique(t)], function(ti, lg, dt_mean, dt_smooth){
      N <- dt_smooth[t == ti, .N]
      muhat <- dt_mean[t == ti, muhat]
      Xvec <- dt_smooth[t == ti, x]
      Xk <- Xvec[1:(N-lg)]
      Xk_plus_lag <- Xvec[(1 + lg):N]
      autocov_vec <- (Xk - muhat) * (Xk_plus_lag - muhat)
      autocov_vec <- autocov_vec[!is.nan(autocov_vec)]
      autocov <- mean(autocov_vec)
      dt_res <- data.table::data.table("t" = ti, "lag" = lg, "autocov" = autocov)
    }, lg = lg, dt_mean = dt_mean, dt_smooth = dt_smooth))
  }, dt_smooth = dt_smooth, dt_mean = dt_mean))
  return(dt_autocov)
}


#' Estimate empirical \eqn{X_0(s)X_{\ell}(t)} autocovariance function
#'
#' @inheritParams .format_data
#' @param s \code{vector (numeric)}. First argument in \eqn{X_0(s)X_{\ell}(t)}.
#' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{t}
#' @param t \code{vector (numeric)}. Second argument in \eqn{X_0(s)X_{\ell}(t)}.
#' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
#' It has to be of the same length as the \code{s}.
#' @param cross_lag \code{integer (positive integer)}. It corresponds to the lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.
#' @param lag \code{vector (integer)}. Lag of the autocovariance of the random variable \eqn{X_0(s)X_{\ell}(t)}.
#' @param h \code{numeric (positive vector or scalar)}. The smoothing bandwidth parameter.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A data.table with three columns: \code{s}, \code{t}, \code{cross_lag}, \code{lag} and \code{XsXt_autocov}
#' corresponding to the estimated autocovariance of the random variable \eqn{X_0(s)X_{\ell}(t)}.
#' @export
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist merge.data.table
#'
#' @seealso [estimate_nw()]
#'

estimate_empirical_XsXt_autocov <- function(data, idcol = NULL, tcol = "tobs", ycol = "X",
                                            s = c(1/5, 2/5, 4/5),
                                            t = c(1/4, 1/2, 3/4),
                                            cross_lag = 1,
                                            lag = c(0, 1, 2), h = NULL,
                                            smooth_ker = epanechnikov){
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
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

  # Sort by s and t
  dt_st <- data.table::data.table("s" = s, "t" = t)
  dt_st <- dt_st[order(s,t)]
  s <- dt_st[, s]
  t <- dt_st[, t]
  rm(dt_st)
  gc()

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
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    if (N > 50) {
      sample_curves <- sample(x = 1:N, size = 30)
    } else {
      sample_curves <- 1:N
    }
    h <- median(unlist(lapply(sample_curves, function(i, lambdahat, data, smooth_ker){
      K <- 100
      b0 <- 2 / lambdahat
      bK <- lambdahat ** (- 1 / 3)
      a <- exp((log(bK) - log(b0)) / K)
      hgrid <- b0 * a ** (seq_len(K))
      hbest <- estimate_nw_bw(y = data[id_curve == i, X],
                              t = data[id_curve ==i, tobs],
                              h_grid = hgrid,
                              smooth_ker = smooth_ker)

    }, lambdahat = lambdahat, data = data, smooth_ker = smooth_ker)))
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Smooth curves
  dt_smooth <- data.table::rbindlist(lapply(1:N, function(curve_index, hvec, s, t, kernel_smooth, data){
    Xhat_s <- estimate_nw(y = data[id_curve == curve_index, X],
                          t = data[id_curve == curve_index, tobs],
                          tnew = s,
                          h = hvec[curve_index],
                          smooth_ker = kernel_smooth)
    Xhat_t <- estimate_nw(y = data[id_curve == curve_index, X],
                          t = data[id_curve == curve_index, tobs],
                          tnew = t,
                          h = hvec[curve_index],
                          smooth_ker = kernel_smooth)
    dt_res <- data.table::data.table("id_curve" = curve_index, "s" = s, "t" = t,
                                     "Xhat_s" = Xhat_s[, yhat], "Xhat_t" = Xhat_t[, yhat])
    return(dt_res)
  }, hvec = h, s = s, t = t, kernel_smooth = smooth_ker, data = data))

  # Take into account the cross_lag
  ## The argument s is associated to the curves n = 1,..., N - cross_lag
  dt_smooth_s <- dt_smooth[, list(id_curve, s, t, Xhat_s)]
  dt_smooth_s <- dt_smooth_s[id_curve %in% 1:(N - cross_lag)]
  dt_id_lag_s <- data.table::data.table(
    "id_curve" = sort(unique(dt_smooth_s[, id_curve])),
    "id_lag" = paste0(1:(N - cross_lag), "_", (1 + cross_lag):N))
  dt_smooth_s <- data.table::merge.data.table(
    x = dt_id_lag_s,
    y = dt_smooth_s,
    by = "id_curve")

  ## The argument t is associated to the curves n = 1 + cross_lag,..., N
  dt_smooth_t <- dt_smooth[, list(id_curve, s, t, Xhat_t)]
  dt_smooth_t <- dt_smooth_t[id_curve %in% (1 + cross_lag):N]
  dt_id_lag_t <- data.table::data.table(
    "id_curve" = sort(unique(dt_smooth_t[, id_curve])),
    "id_lag" = paste0(1:(N - cross_lag), "_", (1 + cross_lag):N))
  dt_smooth_t <- data.table::merge.data.table(
    x = dt_id_lag_t,
    y = dt_smooth_t,
    by = "id_curve")

  ## Merge
  dt_smooth_merge <- data.table::merge.data.table(
    x = dt_smooth_s,
    y = dt_smooth_t,
    by = c("id_lag", "s", "t"))
  dt_smooth_merge <- dt_smooth_merge[order(id_curve.x)]

  ## Clean unnecessary object
  rm(dt_id_lag_s, dt_id_lag_t, dt_smooth_s, dt_smooth_t, dt_smooth)
  gc()

  # Estimate mean function
  dt_gamma_cross_lag <- dt_smooth_merge[!(is.nan(Xhat_s) | is.nan(Xhat_t)),
                                        list("gamma_cross_lag" = mean(Xhat_s * Xhat_t)),
                                        by = c("s", "t")]

  # Estimate autocovariance
  dt_autocov <- data.table::rbindlist(lapply(lag, function(lg, dt_smooth_merge, dt_gamma_cross_lag){
    data.table::rbindlist(lapply(1:length(s), function(sid, lg, dt_gamma_cross_lag, dt_smooth_merge){
      N <- dt_smooth_merge[s == s[sid] & t == t[sid], .N]
      gamma_cross_lag <- dt_gamma_cross_lag[s == s[sid] & t == t[sid], gamma_cross_lag]

      ## For the argument s
      Xhat_s_vec <- dt_smooth_merge[s == s[sid] & t == t[sid], Xhat_s]
      Xhat_s_i <- Xhat_s_vec[1:(N-lg)]
      Xhat_s_i_plus_lag <- Xhat_s_vec[(1 + lg):N]

      ## For the argument t
      Xhat_t_vec <- dt_smooth_merge[s == s[sid] & t == t[sid], Xhat_t]
      Xhat_t_i_plus_cross_lag <- Xhat_t_vec[1:(N-lg)]
      Xhat_t_i_plus_cross_lag_plus_lag <- Xhat_s_vec[(1 + lg):N]

      XsXt_autocov_vec <- (Xhat_s_i * Xhat_t_i_plus_cross_lag - gamma_cross_lag) *
        (Xhat_s_i_plus_lag * Xhat_t_i_plus_cross_lag_plus_lag - gamma_cross_lag)

      XsXt_autocov_vec <- XsXt_autocov_vec[!is.nan(XsXt_autocov_vec)]
      XsXt_autocov <- mean(XsXt_autocov_vec)
      dt_res <- data.table::data.table("s" = s[sid], "t" = t[sid], "cross_lag" = cross_lag, "lag" = lg, "XsXt_autocov" = XsXt_autocov)
    }, lg = lg, dt_gamma_cross_lag = dt_gamma_cross_lag, dt_smooth_merge = dt_smooth_merge))
  }, dt_smooth_merge = dt_smooth_merge, dt_gamma_cross_lag = dt_gamma_cross_lag))
  return(dt_autocov)
}

