#' Estimate the risk of the mean function
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each \code{t}.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it will be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the estimated risk function for mean function estimation.
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()].
#'
#' @importFrom methods is
#' @importFrom data.table data.table rbindlist
#'
estimate_mean_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                               t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
                               Delta = NULL, h = NULL, smooth_ker = epanechnikov){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate local regularity parameters
  # This function controls the remaining arguments
  dt_locreg <- estimate_locreg(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t, Delta = Delta, h = h,
    smooth_ker = smooth_ker)

  # Estimation of the observation error standard deviation
  dt_sigma <- estimate_sigma(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t)

  # Estimation of the variance using the presmoothing bandwidth
  dt_autocov <- estimate_empirical_autocov(
    data = data, idcol = "id_curve",
    tcol = "tobs", ycol = "X", t = t, lag = 0:(N-1),
    h = dt_locreg[, unique(h)],
    smooth_ker = smooth_ker)

  # Estimate the risk function
  dt_mean_risk <- rbindlist(lapply(bw_grid, function(h, t, H, L, kernel_smooth, data, sig_error, N, dt_autocov){
    # compute quantities to estimate estimate the risk
    dt_risk <- data.table::rbindlist(lapply(1:N, function(curve_index, t, h, H, kernel_smooth, data){
      # Compute the weight of the estimator
      Tn <- data[id_curve == curve_index, tobs]
      ker <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, Ker) Ker((Tnm - ti) / h), Ker = kernel_smooth)
      ## By convention NaN = 0
      wker <- apply(ker, 1, function(r) r / ifelse(sum(r) != 0, sum(r), 1))
      wker <- t(wker)

      # Take the maximum and c_n(t,h), the sum of the weight
      wmax <- apply(wker, 1, max)
      cn <- apply(wker, 1, function(r) sum(abs(r)))

      Tn_t_2H <- outer(
        X = 1:length(t), Y = Tn,
        FUN = function(tid, Tnm, t, H) abs((Tnm - t[tid]) / h) ** (2 * H[tid]),
        t = t, H = H
      )

      bn2H <- diag(abs(wker) %*% t(Tn_t_2H))

      # \pi_n(t;h)
      pi_n <- outer(X = t, Y = Tn, FUN = function(ti, Tnm, h) as.numeric(abs(Tnm - ti) <= h), h = h)
      pi_n <- as.numeric(rowSums(pi_n) >= 1)

      # data.table return
      dt_res <- data.table::data.table(
        "id_curve" = curve_index, "t" = t, "wmax" = wmax,
        "cn" = cn, "bn2H" = bn2H, "pi_n" = pi_n
      )
      return(dt_res)
    }, t = t, h = h, H = H, kernel_smooth = kernel_smooth, data = data))

    # Compute P_N(t, h)
    dt_risk[, PN := sum(pi_n), by = t]

    # compute \mathbb B(t,h, 2H) for each n = 1, ..., N
    dt_risk[, B := (sum(pi_n * cn * bn2H) / PN) ** 2, by = t]

    # Compute \mathbb V_mu(t,h) for each n = 1, ..., N
    dt_risk[, Vmu := sum(pi_n * cn * wmax) / PN ** 2, by = t]

    # Compute the risk for each h
    dt_rk <- unique(dt_risk[, .(t, B, Vmu, PN)])

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
        pi_vec <- dt_rho_ell[t == ti, pi_n]
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
    dependence_term <- 2 * dependence_coef  / dt_rk[, PN]

    ## Final risk function
    mean_risk <- bias_term + varriance_term + dependence_term

    # Result to returned
    dt_res <- data.table::data.table("t" = t, "h" = h, "mean_risk" = mean_risk,
                                     "bias_term" = bias_term, "varriance_term" = varriance_term,
                                     "dependence_term" = dependence_term, "H" = H, "L" = L)
    return(dt_res)
  }, t = t, H = dt_locreg[, H], L = dt_locreg[, L], kernel_smooth = smooth_ker, data = data,
  sig_error = dt_sigma[, sig], N = N, dt_autocov = dt_autocov))

  return(dt_mean_risk)
}


#' Estimate mean function
#'
#' @inheritParams estimate_mean_risk
#'
#' @return A \code{data.table} containing the estimated mean function.
#' @export
#'
#' @examples
estimate_mean <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
                          Delta = NULL, h = NULL, smooth_ker = epanechnikov){

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate the risk function
  dt_mean_risk <- estimate_mean_risk(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol,
    t = t, bw_grid = bw_grid,
    Delta = Delta, h = Delta, smooth_ker = smooth_ker)

  # Take the optimum of the risk function
  dt_mean_optbw <- dt_mean_risk[, list("optbw" = h[which.min(mean_risk)]), by = t]

  # Smooth curves with optimal bandwidth parameters
  dt_Xhat <- data.table::rbindlist(lapply(1:N, function(curve_index, t, data, dt_mean_optbw, smooth_ker){
    # \pi_n(t,h)
    Tn <- data[id_curve == curve_index, tobs]
    Yn <- data[id_curve == curve_index, X]
    pi_n <- sapply( X = t, function(ti, Tn, dt_mean_optbw){
      as.numeric(abs(Tn - ti) <= dt_mean_optbw[t == ti, optbw])
    }, Tn = Tn, dt_mean_optbw = dt_mean_optbw)
    pi_n <- t(pi_n)
    pi_n <- as.numeric(rowSums(pi_n, na.rm = TRUE) >= 1)

    # \widehat X(t;h)
    Xhat <- mapply(function(t, h, Yn, Tn, ker){
      estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
    }, t = dt_mean_optbw[, t], h = dt_mean_optbw[, optbw],
    MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

    dt_res <- data.table("id_curve" = curve_index, "t" = t, "Xhat" = Xhat, "pi_n" = pi_n)
    return(dt_res)
  }, data = data, t = t, dt_mean_optbw = dt_mean_optbw, smooth_ker = smooth_ker))

  # Estimate mean function \widehat \mu(t)
  dt_Xhat[is.nan(Xhat) & pi_n == 0, Xhat := 0]
  dt_Xhat[, PN := sum(pi_n), by = t]
  dt_Xhat <- dt_Xhat[!is.nan(Xhat)]
  dt_Xhat[, "muhat" := sum(pi_n * Xhat) / PN, by = t]

  # Add the optimal bandwidth
  dt_muhat <- unique(dt_Xhat[, .(t, muhat, PN)])
  dt_muhat <- data.table::merge.data.table(x = dt_muhat, y = dt_mean_optbw, by = "t")

  return(dt_muhat)
}

# mean function estimator : Rubìn et Paranaretos ----
# Following the Rubìn and Panaretos Equation (B.1), we define S_r_fun and Q_r_fun

#' Sr function. See Rubìn and Panaretos (2020) Equation (B.1)
#'
#' @inheritParams .format_data
#' @param t \code{numeric (positive scalar)}. Observation point at which we want to estimate the mean function.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each \code{t}.
#' @param r \code{numeric (integer)}. It is used as exponent.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{numeric} scalar.
#' @export
#'
.Sr_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    t = 1/2, r = 1, h, smooth_ker = epanechnikov){

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Extract all {T_{n,k}}
  Tn <- data[order(tobs), tobs]
  rm(data)

  # Compute S_r
  Sr <- sum( ((Tn - t) ** r) * (1 / h) * smooth_ker((Tn - t) / h)) / N

  return(Sr)
}

#' Qr function. See Rubìn and Panaretos (2020) Equation (B.1)
#'
#' @inheritParams .format_data
#' @param t \code{numeric (positive scalar)}. Observation point at which we want to estimate the mean function.
#' @param r \code{numeric (integer)}. It is used as exponent.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{numeric} scalar.
#' @export
#'
.Qr_fun <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                    t = 1/2, r = 0, h, smooth_ker = epanechnikov){
  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Extract all {T_{n,k}} and {Y_{n,k}}
  Tn <- data[order(tobs), tobs]
  Yn <- data[order(tobs), X]
  rm(data)

  # Compute Q_r
  Qr <- sum( ((Tn - t) ** r) * Yn * (1 / h) * smooth_ker((Tn - t) / h)) / N

  return(Qr)
}

#' Estimate mean function using mean Rubìn et Paranaretos (2020) method
#'
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#' @return A \code{data.table} containing among other the estimates of the mean function a each \code{t}.
#' @export
#'
estimate_mean_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             t = c(1/4, 1/2, 3/4), h, smooth_ker = epanechnikov){
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  muhat <- sapply(t, function(ti, data, idcol, tcol, ycol, h, ker){
    Q0 <- .Qr_fun(data = data, idcol = idcol, tcol = tcol, ycol = ycol, t = ti, r = 0, h = h, smooth_ker = ker)
    Q1 <- .Qr_fun(data = data, idcol = idcol, tcol = tcol, ycol = ycol, t = ti, r = 1, h = h, smooth_ker = ker)
    S0 <- .Sr_fun(data = data, idcol = idcol, tcol = tcol, ycol = ycol, t = ti, r = 0, h = h, smooth_ker = ker)
    S1 <- .Sr_fun(data = data, idcol = idcol, tcol = tcol, ycol = ycol, t = ti, r = 1, h = h, smooth_ker = ker)
    S2 <- .Sr_fun(data = data, idcol = idcol, tcol = tcol, ycol = ycol, t = ti, r = 2, h = h, smooth_ker = ker)
    muhat_ti <- (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)
    return(muhat_ti)
  }, data = data, idcol = idcol, tcol = tcol, ycol = ycol, h = h, ker = smooth_ker, simplify = TRUE)

  ## Return result
  dt_res <- data.table::data.table("t" = t, "h" = h, "muhat_RP" = muhat)
  return(dt_res)
}


# Bandwidth selection for R-P method
## NB : the bandwidth is global
get_mean_cvbw_rp <- function(data, Kfold = 10, Bmu_grid = seq(0.001, 0.15, len = 45), kernel_smooth = epanechnikov){
  # Input :
  #                data : sample of FTS
  #               Kfold : the number of fold for the cross validation
  #            Bmu_grid : the bandwidth grid
  #       kernel_smooth : the smoothing kernel
  #
  #              Output : data.table containing, among other, the estimates of the risk at each Bmu0 in Bmu_grid

  # Create Kfold folds
  fold <- caret::createFolds(y = 1:length(data), k = Kfold, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(lapply(Bmu_grid, function(Bmu0, data, fold, kernel_smooth){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, Bmu0, kernel_smooth){
        # split train - test
        dt_test <- data.table::rbindlist(data[f])
        dt_train <- data[setdiff(unlist(fold), f)]

        # Estimation of mean on fold\f and test on f
        dt_mu <- estimate_mean_rp(data = dt_train,
                                  x = dt_test[, t],
                                  Bmu = Bmu0,
                                  kernel_smooth = kernel_smooth)

        Sqerror <- (dt_test[, x] - dt_mu[, muhat_RP]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, Bmu0 = Bmu0, kernel_smooth = kernel_smooth, simplify = TRUE),
      error = function(e){
        message("Error in estimating the mean function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold, na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("Bmu0" = Bmu0, "cv_error" = cv_err)
    return(dt_res)

  }, data = data, fold = fold, kernel_smooth = kernel_smooth))

  return(dt_bw)
}
