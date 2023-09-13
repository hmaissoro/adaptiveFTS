#' Estimate the risk of the mean function
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each \code{t}.
#' @param Delta \code{numeric (positive)}. The length of the neighbourhood of \code{t} around which the local regularity is to be estimated.
#' Default \code{Delta = NULL} and thus it we be estimated from the data.
#' @param h \code{numeric (positive vector or scalar)}. The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
#' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
#' If \code{h} is a \code{scalar}, then all curves will be smoothed with the same bandwidth.
#' Otherwise, if \code{h} is a \code{vector}, its length must be equal to the number of curves in \code{data}
#' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return
#' @export
#'
#' @importFrom methods is
#'
#' @examples
estimate_mean_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                               t = 1/2, bw_grid = seq(0.005, 0.15, len = 45),
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
  dt_var_x <- data.table::rbindlist(lapply(1:N, function(curve_index, presmooth_bw, t, kernel_smooth, data){
      x_smooth <- estimate_nw(y = data[id_curve == curve_index, X],
                              t = data[id_curve == curve_index, tobs],
                              tnew = t,
                              h = presmooth_bw,
                              smooth_ker = kernel_smooth)
      dt_smooth <- data.table("curve_index" = curve_index, "t" = t, "x" = x_smooth$yhat)
      return(dt_smooth)
    }, presmooth_bw = dt_locreg[, unique(h)], t = t, kernel_smooth = smooth_ker, data = data))
  dt_var_x <- dt_var_x[!is.nan(x)]
  dt_var_x <- dt_var_x[, .("var_x" = var(x)), by = t]

  # Estimate the risk function
  dt_Rmu <- rbindlist(lapply(bw_grid, function(h, t, H, L, kernel_smooth, data, sig_error, var_x, N){
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
      dt_res <- data.table(
        "curve_index" = curve_index, "t" = t, "wmax" = wmax,
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
    dt_risk <- unique(dt_risk[, .(t, B, Vmu, PN)])

    ## Biais term
    bias_term <- 2 * L * (h ** (2 * H)) * dt_risk[, B]

    ## Variance term
    varriance_term <- 2 * (sig_error ** 2) * dt_risk[, Vmu]

    ## Convergence term to true mean
    ker_far_norm <- 0.5
    dependence_coef <- 2 * max(t ** H) * ker_far_norm / (1 - ker_far_norm) ** 2
    dependence_term <- 2 * sqrt(var_x) * (sqrt(var_x) + 2 * dependence_coef)  / dt_risk[, PN]

    ## Final risk function
    Rmu <- bias_term + varriance_term + dependence_term

    # Result to returned
    dt_res <- data.table("t" = t, "h" = h, "Rmu" = Rmu)
    return(dt_res)
  }, t = t, H = dt_locreg[, H], L = dt_locreg[, L], kernel_smooth = smooth_ker, data = data,
  sig_error = dt_sigma[, sig], var_x = dt_var_x[, var_x], N = N))

  return(dt_Rmu)

}

estimate_mean <- function(data, t, Delta, bw_grid = seq(0.005, 0.15,len = 45), kernel_smooth = epanechnikov){
  # Input :
  #           data : data, sample of FTS
  #              t : where we aim to estimate the local regularity and mean function
  #          Delta : define the length of the interval around which the local regularity is estimated
  #         bw_grid : the bandwidth grid
  #         kernel : kernel of the N-W smoother
  # Output : data.table containing estimates of the mean function

  # Estimate the risk function
  dt_risk <- mean_risk(data = data, t = t, bw_grid = bw_grid, Delta = Delta,
                       kernel_smooth = kernel_smooth)

  ## Take the optimum bw
  dt_hmu <- dt_risk[, .("hmu" = h[which.min(Rmu)]), by = t]

  # Estimate
  dt_Xhat <- data.table::rbindlist(lapply(1:length(data), function(curve_index, t, data, dt_hmu){

    # w_n(t,h)
    Tn <- data[[curve_index]]$t
    Xn <- data[[curve_index]]$x
    pi_n <- sapply(
      X = t,
      function(ti, Tn, dt_hmu) as.numeric(abs(Tn - ti) <= dt_hmu[t == ti, hmu]),
      Tn = Tn, dt_hmu = dt_hmu
    )
    pi_n <- t(pi_n)
    pi_n <- as.numeric(rowSums(pi_n, na.rm = TRUE) >= 1)

    Xhat <- mapply(function(t, h, Xn, Tn, ker){
      nw(x = Tn, y = Xn, xout = t, h = h, kernel = ker)$yhat
    }, t = dt_hmu[, t], h = dt_hmu[, hmu],
    MoreArgs = list(Xn = Xn, Tn = Tn, ker = kernel_smooth))

    dt_res <- data.table("curve_index" = curve_index, "t" = t, "Xhat" = Xhat, "pi_n" = pi_n)
    return(dt_res)
  }, data = data, t = t, dt_hmu = dt_hmu))

  # Estimate mean function \widehat \mu(t)
  dt_Xhat[, P_N := sum(pi_n), by = t]
  dt_Xhat[, "muhat" := sum(pi_n * Xhat) / P_N, by = t]

  # Add the optimal bandwidth
  dt_muhat <- unique(dt_Xhat[, .(t, muhat, P_N)])
  dt_muhat <- data.table::merge.data.table(x = dt_muhat, y = dt_hmu, by = "t")

  return(dt_muhat)
}

# mean function estimator : Rubìn et Paranaretos ----
# Following the Rubìn and Panaretos Equation (B.1), we define S_r_fun and Q_r_fun

.Sr_fun <- function(data, x, r = 1, Bmu, kernel_smooth = epanechnikov){
  # Input :
  #          data : sample of FTS
  #             x : a point at which we want to estimate the mean function
  #             r : an integer, used as exponent
  #           Bmu : a scalar bandwidth
  #        kernel_smooth : the smoothing kernel
  #
  # Output : a real number

  # Extract X_all and N (the number of curves)
  N <- length(data)
  data_all <- data.table::rbindlist(data)
  data_all <- data_all[order(t)]
  X_all <- data_all[, t]

  rm(data_all)

  # Compute S_r
  Sr <- sum( ((X_all - x) ** r) * (1 / Bmu) * kernel_smooth((X_all - x) / Bmu)) / N

  return(Sr)
}

.Qr_fun <- function(data, x, r = 0, Bmu, kernel_smooth = epanechnikov){
  # Input :
  #          data : sample of FTS
  #             x : a point at which we cant to estimate the mean function
  #             r : an integer, used as exponent
  #           Bmu : a scalar bandwidth
  #        kernel : the smoothing kernel
  #
  # Output :  a vector of real values

  # Extract X_all, Y_all and N
  N <- length(data)
  data_all <- data.table::rbindlist(data)
  data_all <- data_all[order(t)]
  X_all <- data_all[, t]
  Y_all <- data_all[, x]

  rm(data_all)

  # Compute Q_r
  Qr <- sum( ((X_all - x) ** r) * Y_all * (1 / Bmu) * kernel_smooth((X_all - x) / Bmu)) / N

  return(Qr)
}

estimate_mean_rp <- function(data, x, Bmu, kernel_smooth = epanechnikov){
  # Input :
  #          data : sample of FTS
  #             x : vector of points at which we cant to estimate the mean function
  #           Bmu : a scalar bandwidth
  # kernel_smooth : the smoothing kernel
  #
  # Output : data.table containing among other the estimates of the mean function a each x
  #

  muhat <- sapply(x, function(xi, data = data, Bmu = Bmu, kernel_smooth){
    Q0 <- .Qr_fun(data = data, x = xi, r = 0, Bmu = Bmu, kernel_smooth = kernel_smooth)
    Q1 <- .Qr_fun(data = data, x = xi, r = 1, Bmu = Bmu, kernel_smooth = kernel_smooth)
    S0 <- .Sr_fun(data = data, x = xi, r = 0, Bmu = Bmu, kernel_smooth = kernel_smooth)
    S1 <- .Sr_fun(data = data, x = xi, r = 1, Bmu = Bmu, kernel_smooth = kernel_smooth)
    S2 <- .Sr_fun(data = data, x = xi, r = 2, Bmu = Bmu, kernel_smooth = kernel_smooth)
    muhat_xi <- (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)
    return(muhat_xi)
  }, data = data, Bmu = Bmu, kernel_smooth = kernel_smooth, simplify = TRUE)

  ## Return result
  dt_res <- data.table::data.table("t" = x, "Bmu" = Bmu, "muhat_RP" = muhat)
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
