#' Estimate the Risk of the Mean Function Estimator
#'
#' Estimates the risk \eqn{R_\mu(t;h)} of the adaptive mean function estimator
#' over a grid of candidate bandwidths, as described in Section 4.1 of
#' Maissoro, Patilea and Vimond (2025). Minimising it over
#' \code{h} at each \code{t} is what \link{estimate_mean} does to select its
#' bandwidth.
#'
#' @details
#' The risk bound splits into three terms, returned separately so that the
#' selected bandwidth can be traced back to what drove it: a bias term growing
#' with \eqn{h^{2H_t}} through the local regularity, a variance term decreasing in
#' \eqn{h} through the number of usable points, and a dependence term reflecting
#' the serial dependence between curves. The local regularity parameters are
#' estimated internally at each \code{t}, so \code{Ht} and \code{Lt2} are reported
#' alongside the risk.
#'
#' Left to \code{NULL}, \code{bw_grid} is a 20-point geometric grid running from
#' \eqn{4(N\widehat\lambda)^{-0.9}} to \eqn{4(N\widehat\lambda)^{-1/3}}, where
#' \eqn{N} is the number of curves and \eqn{\widehat\lambda} the average number of
#' observation points per curve.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Points of \eqn{[0, 1]} at which the risk is
#' estimated.
#' @param bw_grid \code{vector (numeric)}. Candidate bandwidths, from which
#' \link{estimate_mean} picks the risk-minimising one at each \code{t}. Default
#' \code{NULL} builds the grid from the data; see Details.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#' @inheritParams estimate_locreg
#'
#' @return A \code{data.table} with one row per (\code{t}, \code{h}) pair and
#' columns:
#' \itemize{
#'   \item \code{t}: the point at which the risk is estimated.
#'   \item \code{h}: the candidate bandwidth.
#'   \item \code{PN}: the number of curves contributing to the estimate at
#'     \code{t}, \eqn{P_N(t;h)}.
#'   \item \code{locreg_bw}: the bandwidth used to estimate the local regularity.
#'   \item \code{Ht}: the estimated local exponent \eqn{H_t}.
#'   \item \code{Lt2}: the estimated squared Hölder constant \eqn{L_t^2}.
#'   \item \code{bias_term}, \code{variance_term}, \code{dependence_term}: the
#'     three components of the risk.
#'   \item \code{mean_risk}: the estimated risk.
#' }
#'
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()].
#'
#' @import data.table
#' @importFrom methods is
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#' @examples
#' data("data_far")
#'
#' dt_mean_risk <- estimate_mean_risk(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.02, 0.15, length.out = 8),
#'   kernel_name = "epanechnikov")
#'
#' # The risk-minimising bandwidth at each t.
#' dt_mean_risk[, list(h = h[which.min(mean_risk)]), by = "t"]
#'
estimate_mean_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                               t = c(1/4, 1/2, 3/4), bw_grid = NULL,
                               presmooth_bw = NULL, Delta = NULL,
                               presmooth_bw_grid = NULL, presmooth_nsubset = NULL,
                               kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  .check_locreg_args(presmooth_bw = presmooth_bw, Delta = Delta,
                     presmooth_bw_grid = presmooth_bw_grid,
                     presmooth_nsubset = presmooth_nsubset)

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
    bw_grid <- .default_bw_grid(data)
  }

  # Estimate risk function using C++  function
  dt_mean_risk <- estimate_mean_risk_cpp(data = data, t = t, bw_grid = bw_grid, kernel_name = kernel_name,
                                         presmooth_bw = presmooth_bw, Delta = Delta,
                                         presmooth_bw_grid = presmooth_bw_grid,
                                         presmooth_nsubset = presmooth_nsubset)
  dt_mean_risk <- data.table::as.data.table(dt_mean_risk)
  data.table::setnames(x = dt_mean_risk,
                       new = c("t", "h", "PN", "locreg_bw", "Ht", "Lt2", "bias_term",
                               "variance_term", "dependence_term", "mean_risk"))

    return(.as_adaptive_est(dt_mean_risk, "mean_risk",
                            meta = list(kernel = kernel_name, N = N,
                                        n_bw = length(bw_grid))))
}


#' Estimate the Mean Function
#'
#' Estimates the mean function of the underlying process with the adaptive
#' estimator of Maissoro, Patilea and Vimond (2025), using at
#' each point the bandwidth that minimises the estimated risk.
#'
#' @details
#' Unless \code{bw} is supplied, the bandwidth is selected point by point by
#' minimising the risk of \link{estimate_mean_risk} over \code{bw_grid}, so
#' neighbouring points may be smoothed differently according to the local
#' regularity of the process. Supplying \code{bw} skips the risk estimation
#' entirely, which is worth doing when the same bandwidths are reused across
#' many calls.
#'
#' @inheritParams estimate_mean_risk
#' @param bw \code{vector (numeric)}. Bandwidth to use at each point of \code{t},
#' recycled if a scalar. Default \code{NULL} selects it by minimising the risk
#' estimated by \link{estimate_mean_risk}.
#'
#' @return A \code{data.table} with one row per point of \code{t} and columns:
#' \itemize{
#'   \item \code{t}: the point at which the mean function is estimated.
#'   \item \code{optbw}: the bandwidth used at \code{t}.
#'   \item \code{Ht}: the estimated local exponent \eqn{H_t}.
#'   \item \code{Lt2}: the estimated squared Hölder constant \eqn{L_t^2}.
#'   \item \code{PN}: the number of curves contributing to the estimate at
#'     \code{t}.
#'   \item \code{muhat}: the estimated mean function.
#' }
#'
#' @export
#'
#' @seealso [estimate_mean_risk()], [estimate_locreg()], [estimate_autocov()].
#'
#' @import data.table
#'
#' @references
#' Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
#' Weakly Dependent Functional Time Series. \emph{Journal of Time Series
#' Analysis}. \doi{10.1111/jtsa.70006}
#'
#' @examples
#' data("data_far")
#'
#' dt_mean <- estimate_mean(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.02, 0.15, length.out = 8),
#'   kernel_name = "epanechnikov")
#' dt_mean
#'
#' summary(dt_mean)
#'
estimate_mean <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = c(1/4, 1/2, 3/4), bw = NULL, bw_grid = NULL,
                          presmooth_bw = NULL, Delta = NULL,
                          presmooth_bw_grid = NULL, presmooth_nsubset = NULL,
                          kernel_name = "epanechnikov"){
  # Control on t, bw and kernel_name arguments
  # NB : The remaining arguments are controlled using the format_data and estimate_mean_risk functions, if required.
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  .check_locreg_args(presmooth_bw = presmooth_bw, Delta = Delta,
                     presmooth_bw_grid = presmooth_bw_grid,
                     presmooth_nsubset = presmooth_nsubset)

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate mean function using C++  function
  dt_muhat <- estimate_mean_cpp(data = data, t = t, bw = bw, bw_grid = bw_grid, kernel_name = kernel_name,
                                presmooth_bw = presmooth_bw, Delta = Delta,
                                presmooth_bw_grid = presmooth_bw_grid,
                                presmooth_nsubset = presmooth_nsubset)
  dt_muhat <- data.table::as.data.table(dt_muhat)
  data.table::setnames(x = dt_muhat, new = c("t", "optbw", "Ht", "Lt2", "PN", "muhat"))
  return(.as_adaptive_est(dt_muhat, "mean_est",
                          meta = list(kernel = kernel_name, N = N)))
}


#' Estimate the Mean Function by the Rubìn-Panaretos Method
#'
#' Estimates the mean function with the local-linear estimator of
#' Rubìn and Panaretos (2020), which pools the observation points
#' of all curves and smooths them with a single bandwidth. It is provided for
#' comparison with the adaptive estimator of \link{estimate_mean}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Points of \eqn{[0, 1]} at which the mean
#' function is estimated.
#' @param bw \code{numeric (positive scalar)}. Bandwidth of the estimator, common
#' to every point of \code{t}. See \link{estimate_mean_bw_rp} to select it by
#' cross-validation.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#'
#' @return A \code{data.table} with one row per point of \code{t} and columns:
#' \itemize{
#'   \item \code{t}: the point at which the mean function is estimated.
#'   \item \code{bw}: the bandwidth used.
#'   \item \code{muhat_RP}: the estimated mean function.
#' }
#' @export
#'
#' @seealso [estimate_mean_bw_rp()], [estimate_mean()].
#'
#' @import data.table
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @examples
#' data("data_far")
#'
#' dt_mean_rp <- estimate_mean_rp(
#'   data = data_far[data_far$id_curve <= 20, ],
#'   idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw = 5/70, kernel_name = "epanechnikov")
#' dt_mean_rp
#'
estimate_mean_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             t = c(1/4, 1/2, 3/4), bw, kernel_name = "epanechnikov"){
  smooth_ker <- .select_kernel(kernel_name)

  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]
  lambdahat <- mean(data[, .N, by = "id_curve"][, N])

  if (N <= 250 && lambdahat < 200) {
    # Estimate mean function
    Tn <- data[order(tobs), tobs]
    Yn <- data[order(tobs), X]
    data_curve <- kronecker(
      X = matrix(data = rep(1, length(t)), ncol = 1),
      Y = cbind(Tn, Yn)
    )
    colnames(data_curve) <- c("Tn", "Yn")
    data_curve <- data.table::as.data.table(data_curve)
    tvec <- rep(t, each = length(Tn))
    data_curve[, t := tvec]
    rm(Yn, Tn, tvec, data) ; gc() ; gc()

    data_curve[, Tn_minus_t := (Tn - t)]
    data_curve[, Tn_minus_t_over_bw := Tn_minus_t / bw]

    # Compute mean Q and S function
    dt_res_by_curve <- data_curve[
      ,
      .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
        "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
        "S0" = sum((Tn_minus_t ** 0) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
        "S1" = sum((Tn_minus_t ** 1) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
        "S2" = sum((Tn_minus_t ** 2) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N),
      by = "t"
    ]
    rm(data_curve) ; gc() ; gc()

    # Estimate mean
    dt_res <- dt_res_by_curve[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
    dt_res[, "bw" := bw]
    data.table::setcolorder(x = dt_res, neworder = c("t", "bw", "muhat_RP"))
    rm(dt_res_by_curve) ; gc() ; gc()

  } else {
    # Split t if N x \lambda >> 0
    if (N <= 450 && lambdahat <= 300) {
      N_t_by_list <- 40 * 300 / 10
      t_list <- split(t, ceiling(seq_along(t) / N_t_by_list))
    } else if (N <= 1000 & lambdahat <= 50){
      N_t_by_list <- 1000 * 40 / 50
      t_list <- split(t, ceiling(seq_along(t) / N_t_by_list))
    } else {
      t_list <- t
    }

    # Estimate mean function
    dt_res <- data.table::rbindlist(lapply(t_list, function(t_list_i, data, N, bw){
      Tn <- data[order(tobs), tobs]
      Yn <- data[order(tobs), X]
      if (length(t_list_i) > 1) {
        data_curve <- kronecker(
          X = matrix(data = rep(1, length(t_list_i)), ncol = 1),
          Y = cbind(Tn, Yn)
        )
        colnames(data_curve) <- c("Tn", "Yn")
        data_curve <- data.table::as.data.table(data_curve)
        tvec <- rep(t_list_i, each = length(Tn))
        data_curve[, t := tvec]
        rm(Yn, Tn, tvec) ; gc() ; gc()
      } else {
        data_curve <- data.table::data.table("Tn" = Tn, "Yn" = Yn, "t" = t_list_i)
      }

      data_curve[, Tn_minus_t := (Tn - t)]
      data_curve[, Tn_minus_t_over_bw := Tn_minus_t / bw]

      # Compute mean Q and S function
      dt_res_by_t_list_i <- data_curve[
        ,
        .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
          "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
          "S0" = sum((Tn_minus_t ** 0) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
          "S1" = sum((Tn_minus_t ** 1) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N,
          "S2" = sum((Tn_minus_t ** 2) * (1 / bw) * smooth_ker(Tn_minus_t_over_bw)) / N),
        by = "t"
      ]
      rm(data_curve) ; gc() ; gc()
      # Estimate mean
      dt_res_by_t_list_i <- dt_res_by_t_list_i[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
      dt_res_by_t_list_i[, "bw" := bw]
      data.table::setcolorder(x = dt_res_by_t_list_i, neworder = c("t", "bw", "muhat_RP"))

      # Return
      return(dt_res_by_t_list_i)
    }, data = data, N = N, bw = bw))
  }
  return(dt_res)
}

#' Select the Bandwidth of the Rubìn-Panaretos Mean Estimator
#'
#' Selects the bandwidth of \link{estimate_mean_rp} by \eqn{K}-fold
#' cross-validation over the curves, as described in
#' Rubìn and Panaretos (2020). Curves are split into folds; each
#' fold is predicted from the mean function estimated on the others, and the
#' squared prediction errors are averaged.
#'
#' @inheritParams format_data
#' @param n_folds \code{integer (positive)}. Number of cross-validation folds.
#' @param bw_grid \code{vector (numeric)}. Candidate bandwidths.
#' @param kernel_name \code{string}. Kernel of the smoothing estimator, one of
#' "epanechnikov" (default), "biweight", "triweight", "tricube", "triangular" and
#' "uniform".
#'
#' @return A \code{data.table} with one row per candidate bandwidth and columns:
#' \itemize{
#'   \item \code{bw}: the candidate bandwidth.
#'   \item \code{cv_error}: the cross-validation error at \code{bw}. The
#'     bandwidth minimising it is the one to pass to \link{estimate_mean_rp}.
#' }
#' @export
#' @seealso [estimate_mean_rp()].
#'
#' @import data.table
#'
#' @references
#' Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
#' series: estimation and prediction. \emph{Electronic Journal of Statistics},
#' 14(1), 1137--1210. \doi{10.1214/20-EJS1690}
#'
#' @examples
#' \donttest{
#' data("data_far")
#' dt_small <- data_far[data_far$id_curve <= 10, ]
#'
#' dt_bw <- estimate_mean_bw_rp(
#'   data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   n_folds = 5, bw_grid = seq(0.02, 0.15, length.out = 5),
#'   kernel_name = "epanechnikov")
#'
#' dt_mean_rp <- estimate_mean_rp(
#'   data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw = dt_bw[, bw[which.min(cv_error)]],
#'   kernel_name = "epanechnikov")
#' dt_mean_rp
#' }
#'
estimate_mean_bw_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                n_folds = 10, bw_grid = seq(0.001, 0.15, len = 45),
                                kernel_name = "epanechnikov"){
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  # Create n_folds folds
  fold <- .create_folds(y = unique(data[, id_curve]), k = n_folds, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(lapply(bw_grid, function(Bmu0, data, fold, kernel_name){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, Bmu0, kernel_name){
        # split train - test
        dt_test <- data[id_curve %in% unlist(f)]
        dt_test <- dt_test[order(tobs)]
        dt_train <- data[id_curve %in% setdiff(unlist(fold), unlist(f))]
        dt_train <- dt_train[order(tobs)]

        # Estimation of mean on fold\f and test on f
        dt_mu <- estimate_mean_rp(
          data = dt_train, idcol = "id_curve", tcol = "tobs", ycol = "X",
          t = dt_test[, tobs], bw = Bmu0, kernel_name = kernel_name)

        Sqerror <- (dt_test[, X] - dt_mu[, muhat_RP]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, Bmu0 = Bmu0, kernel_name = kernel_name, simplify = TRUE),
      error = function(e){
        message("Error in estimating the mean function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold[!is.nan(err_fold)], na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("bw" = Bmu0, "cv_error" = cv_err)
    return(dt_res)

  }, data = data, fold = fold, kernel_name = kernel_name))
  rm(data, fold) ; gc() ; gc()

  return(dt_bw)
}
