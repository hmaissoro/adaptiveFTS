#' Estimate the Risk Function of the Mean Function
#'
#' This function estimates the risk function \eqn{R_\mu(t;h)} for the mean function estimation as described in
#' Section 4.1 of \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points where the mean function of the underlying process is estimated.
#' @param bw_grid \code{vector (numeric)}. A bandwidth grid from which the best smoothing parameter is selected for each \code{t}.
#' Default is \code{NULL}, in which case it is defined as an exponential grid of \eqn{N \lambda}.
#' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
#' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item{\code{t} :}{ The observation points where the risk function is estimated.}
#'   \item{\code{h} :}{ The candidate bandwidth values tested.}
#'   \item{\code{PN} :}{ The number of curves used to estimate the mean at each \code{t}, corresponding to \eqn{P_N(t;h)}.}
#'   \item{\code{locreg_bw} :}{ The bandwidth used to estimate the local regularity parameters.}
#'   \item{\code{Ht} :}{ Estimates of the local exponent at each \code{t}, corresponding to \eqn{H_t}.}
#'   \item{\code{Lt} :}{ Estimates of the Hölder constant at each \code{t}, corresponding to \eqn{L_t^2}.}
#'   \item{\code{bias_term} :}{ The bias term component of the risk function.}
#'   \item{\code{variance_term} :}{ The variance term component of the risk function.}
#'   \item{\code{dependence_term} :}{ The dependence term component of the risk function.}
#'   \item{\code{mean_risk} :}{ The estimated risk function for the mean.}
#' }
#'
#' @export
#' @seealso [estimate_mean()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
#'
#' @import data.table
#' @importFrom Rdpack reprompt
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
#' # Estimate the risk function for mean estimation
#' dt_mean_risk <- estimate_mean_risk(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = NULL,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # Plot the mean risk function at different points
#' dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
#'       main = "t = 0.25", xlab = "h", ylab = "Risk Function"
#'     ),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
#'       main = "t = 0.5", xlab = "h", ylab = "Risk Function"
#'     ),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
#'       main = "t = 0.75", xlab = "h", ylab = "Risk Function"
#'     )
#'   ),
#'   nrow = 3
#' )
#' }
#'
estimate_mean_risk <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                               t = c(1/4, 1/2, 3/4), bw_grid = NULL,
                               kernel_name = "epanechnikov"){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

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

  # Estimate risk function using C++  function
  dt_mean_risk <- estimate_mean_risk_cpp(data = data, t = t, bw_grid = bw_grid, kernel_name = kernel_name)
  dt_mean_risk <- data.table::as.data.table(dt_mean_risk)
  data.table::setnames(x = dt_mean_risk,
                       new = c("t", "h", "PN", "locreg_bw", "Ht", "Lt", "bias_term",
                               "variance_term", "dependence_term", "mean_risk"))

    return(dt_mean_risk)
}


#' Estimate Mean Function
#'
#' This function estimates the mean function of an underlying process using the adaptive estimator described in
#' \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
#'
#' @inheritParams estimate_mean_risk
#' @param optbw \code{vector (numeric)}. Optimal bandwidth parameters for mean function estimation at each \code{t}.
#' If \code{optbw = NULL} (default), it will be estimated using the \link{estimate_mean_risk} function.
#'
#' @return A \code{data.table} containing the following columns:
#' \itemize{
#'   \item{\code{t} :}{ The observation points at which the mean function is estimated.}
#'   \item{\code{optbw} :}{ The optimal bandwidth used to estimate the mean function at each \code{t}.}
#'   \item{\code{Ht} :}{ Local exponent estimates for each \code{t}, corresponding to \eqn{H_t}.}
#'   \item{\code{Lt} :}{ Estimates of the Hölder constant for each \code{t}, corresponding to \eqn{L_t^2}.}
#'   \item{\code{PN} :}{ The number of selected curves used in the estimation for each \code{t}.}
#'   \item{\code{muhat} :}{ Estimated values of the mean function at each \code{t}.}
#' }
#'
#' @export
#'
#' @seealso [estimate_mean_risk()], [estimate_locreg()], [estimate_sigma()], [estimate_nw()], [estimate_empirical_autocov()].
#'
#' @import data.table
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Load data
#' data("data_far")
#'
#' # Estimate risk function for the mean
#' dt_mean_risk <- estimate_mean_risk(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = NULL,
#'   kernel_name = "epanechnikov"
#' )
#'
#' # Visualize mean risk at various observation points
#' dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
#' manipulateWidget::combineWidgets(
#'   list = list(
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
#'       main = "t = 0.25", xlab = "h", ylab = "Risk Function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
#'       main = "t = 0.5", xlab = "h", ylab = "Risk Function"),
#'     dygraphs::dygraph(
#'       data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
#'       main = "t = 0.75", xlab = "h", ylab = "Risk Function")
#'   ),
#'   nrow = 3
#' )
#'
#' # Estimate mean function with optimal bandwidths
#' dt_mean <- estimate_mean(
#'   data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
#'   kernel_name = "epanechnikov"
#' )
#'
#' # Display rounded estimates of the mean function
#' DT::datatable(data = dt_mean[, lapply(.SD, function(X) round(X, 3))])
#' }
#'
estimate_mean <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          t = c(1/4, 1/2, 3/4), optbw = NULL, bw_grid = NULL,
                          kernel_name = "epanechnikov"){
  # Control on t, optbw and smooth_ker arguments
  # NB : The remaining arguments are controlled using the format_data and estimate_mean_risk functions, if required.
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")

  # Check the name of the kernel
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform")
  )

  # Control and format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Estimate mean function using C++  function
  dt_muhat <- estimate_mean_cpp(data = data, t = t, optbw = optbw, bw_grid = bw_grid, kernel_name = kernel_name)
  dt_muhat <- data.table::as.data.table(dt_muhat)
  data.table::setnames(x = dt_muhat, new = c("t", "optbw", "Ht", "Lt", "PN", "muhat"))
  return(dt_muhat)
}


#' Estimate mean function using \insertCite{rubin2020;textual}{adaptiveFTS} method.
#'
#' This function estimates the mean function of a set of curves using the method proposed
#' by \insertCite{rubin2020;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the mean function of the underlying process.
#' @param h \code{numeric (positive scalar)}. The bandwidth of the estimator.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{t :}{ The Observation points at which the mean function is estimated.}
#'            \item{h :}{ The bandwidth parameter.}
#'            \item{muhat_RP :}{ The estimates of the mean function using Rubìn and Panaretos (2020) method.}
#'         }
#' @export
#'
#' @seealso [estimate_mean_bw_rp()]
#'
#' @import data.table
#' @import Rdpack
#' @importFrom fastmatrix kronecker.prod
#'
#' @references
#' \insertRef{rubin2020}{adaptiveFTS}
#'
#' @examples
#' \dontrun{
#' # Generate a FAR A process
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
#' # Estimate mean function using Rubìn and Panaretos (2020) method
#' dt_mean_rp <- estimate_mean_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), h = 5/70, smooth_ker = epanechnikov)
#'
#' DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])
#'
#' }
#'
estimate_mean_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                             t = c(1/4, 1/2, 3/4), h, smooth_ker = epanechnikov){
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]
  lambdahat <- mean(data[, .N, by = "id_curve"][, N])

  if (N <= 250 && lambdahat < 200) {
    # Estimate mean function
    Tn <- data[order(tobs), tobs]
    Yn <- data[order(tobs), X]
    data_curve <- fastmatrix::kronecker.prod(
      x = matrix(data = rep(1, length(t)), ncol = 1),
      y = cbind(Tn, Yn)
    )
    colnames(data_curve) <- c("Tn", "Yn")
    data_curve <- data.table::as.data.table(data_curve)
    tvec <- rep(t, each = length(Tn))
    data_curve[, t := tvec]
    rm(Yn, Tn, tvec, data) ; gc() ; gc()

    data_curve[, Tn_minus_t := (Tn - t)]
    data_curve[, Tn_minus_t_over_h := Tn_minus_t / h]

    # Compute mean Q and S function
    dt_res_by_curve <- data_curve[
      ,
      .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S0" = sum((Tn_minus_t ** 0) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S1" = sum((Tn_minus_t ** 1) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
        "S2" = sum((Tn_minus_t ** 2) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N),
      by = "t"
    ]
    rm(data_curve) ; gc() ; gc()

    # Estimate mean
    dt_res <- dt_res_by_curve[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
    dt_res[, "h" := h]
    data.table::setcolorder(x = dt_res, neworder = c("t", "h", "muhat_RP"))
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
    dt_res <- data.table::rbindlist(lapply(t_list, function(t_list_i, data, N, h){
      Tn <- data[order(tobs), tobs]
      Yn <- data[order(tobs), X]
      if (length(t_list_i) > 1) {
        data_curve <- fastmatrix::kronecker.prod(
          x = matrix(data = rep(1, length(t_list_i)), ncol = 1),
          y = cbind(Tn, Yn)
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
      data_curve[, Tn_minus_t_over_h := Tn_minus_t / h]

      # Compute mean Q and S function
      dt_res_by_t_list_i <- data_curve[
        ,
        .("Q0" = sum((Tn_minus_t ** 0) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "Q1" = sum((Tn_minus_t ** 1) * Yn * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S0" = sum((Tn_minus_t ** 0) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S1" = sum((Tn_minus_t ** 1) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N,
          "S2" = sum((Tn_minus_t ** 2) * (1 / h) * smooth_ker(Tn_minus_t_over_h)) / N),
        by = "t"
      ]
      rm(data_curve) ; gc() ; gc()
      # Estimate mean
      dt_res_by_t_list_i <- dt_res_by_t_list_i[, .("muhat_RP" = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 ** 2)), by = "t"]
      dt_res_by_t_list_i[, "h" := h]
      data.table::setcolorder(x = dt_res_by_t_list_i, neworder = c("t", "h", "muhat_RP"))

      # Return
      return(dt_res_by_t_list_i)
    }, data = data, N = N, h = h))
  }
  return(dt_res)
}

#' Bandwidth estimation using cross-validation for the \insertCite{rubin2020;textual}{adaptiveFTS} mean function estimator.
#'
#' This function estimates the optimal bandwidth for the mean function estimator
#' using cross-validation, as described in \insertCite{rubin2020;textual}{adaptiveFTS}.
#'
#' @inheritParams format_data
#' @param Kfold \code{integer (positive)}. Number of fold for the cross-validation.
#' @param bw_grid \code{vector (numeric)}. The bandwidth grid.
#' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator. Default \code{smooth_ker = epanechnikov}.
#'
#' @return A \code{data.table} containing the following columns.
#'          \itemize{
#'            \item{h :}{ The candidate bandwidth.}
#'            \item{cv_error :}{ The estimates of the Cross-Validation error for each \code{h}.}
#'         }
#' @export
#' @seealso [estimate_mean_rp()]
#'
#' @import data.table
#' @import Rdpack
#' @importFrom caret createFolds
#'
#' @references
#' \insertRef{rubin2020}{adaptiveFTS}
#'
#' @examples
#' \dontrun{
#' # Generate a FAR A process
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
#' ## Estimate the bandwidth by Cross-Validation
#' dt_bw_mean_rp <- estimate_mean_bw_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
#'   smooth_ker = epanechnikov)
#'
#' ## Plot the Cross-Validation error
#' dygraphs::dygraph(dt_bw_mean_rp)
#'
#' ## Select the best bandwidth
#' optbw <- dt_bw_mean_rp[, h[which.min(cv_error)]]
#'
#' ## Estimate the mean function
#' dt_mean_rp <- estimate_mean_rp(
#'   data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
#'   t = c(1/4, 1/2, 3/4), h = optbw, smooth_ker = epanechnikov)
#'
#' DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])
#'
#' }
#'
estimate_mean_bw_rp <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
                                smooth_ker = epanechnikov){
  # Format data
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  # Create Kfold folds
  fold <- caret::createFolds(y = unique(data[, id_curve]), k = Kfold, list = TRUE)

  # Get risk for each bandwidth in the grid
  dt_bw <- data.table::rbindlist(lapply(bw_grid, function(Bmu0, data, fold, kernel_smooth){

    # Compute the cross-validation error for each f in fold
    err_fold <- tryCatch(
      expr = sapply(fold, function(f, data, Bmu0, kernel_smooth){
        # split train - test
        dt_test <- data[id_curve %in% unlist(f)]
        dt_test <- dt_test[order(tobs)]
        dt_train <- data[id_curve %in% setdiff(unlist(fold), unlist(f))]
        dt_train <- dt_train[order(tobs)]

        # Estimation of mean on fold\f and test on f
        dt_mu <- estimate_mean_rp(
          data = dt_train, idcol = "id_curve", tcol = "tobs", ycol = "X",
          t = dt_test[, tobs], h = Bmu0, smooth_ker = kernel_smooth)

        Sqerror <- (dt_test[, X] - dt_mu[, muhat_RP]) ** 2
        err <- sum(Sqerror)
        return(err)
      }, data = data, Bmu0 = Bmu0, kernel_smooth = kernel_smooth, simplify = TRUE),
      error = function(e){
        message("Error in estimating the mean function:")
        print(e)
        return(NA)

      })

    # Cross-validaiton error
    cv_err <- mean(err_fold[!is.nan(err_fold)], na.rm = TRUE)

    # Return the result
    dt_res <- data.table::data.table("h" = Bmu0, "cv_error" = cv_err)
    return(dt_res)

  }, data = data, fold = fold, kernel_smooth = smooth_ker))
  rm(data, fold) ; gc() ; gc()

  return(dt_bw)
}

