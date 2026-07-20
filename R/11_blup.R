## Adaptive functional BLUP: fit / predict engine. The design-weighted,
## Tikhonov-regularised BLUP conditions on the previous curve only (single lag-1
## block); the covariance assembly is kept separable so the multi-lag case can be
## added later without a rewrite.

#' Check whether all curves share the same observation design
#'
#' @param data A prepared functional data.table.
#' @param idcol,tcol Identifier and observation-time column names.
#' @return `TRUE` if every curve shares the same sorted design, else `FALSE`.
#' @keywords internal
.is_common_design <- function(data, idcol = "id_curve", tcol = "tobs") {
  design_by_curve <- data[, .(tobs_list = list(sort(unique(get(tcol))))), by = idcol]
  reference_design <- design_by_curve$tobs_list[[1]]
  return(all(vapply(design_by_curve$tobs_list,
                    function(tt) identical(tt, reference_design), logical(1))))
}

#' Fit the adaptive functional BLUP
#'
#' Estimates every component of the adaptive Best Linear Unbiased Predictor that
#' does not depend on the prediction points.
#'
#' @details
#' The fit conditions on a single curve (its immediate successor is what
#' `predict()` reconstructs) and caches the adaptive bandwidths, so that
#' `predict()` only re-runs the cheap plug-in estimates at the requested
#' prediction points. The covariance assembly is a single lag-1 block.
#'
#' @inheritParams format_data
#' @param id_lag Integer id of the conditioning curve. Its successor is the
#'   curve to be predicted. Default `NULL` uses the last curve in `data`.
#' @param tikhonov Tikhonov regularisation parameter \eqn{\alpha} added to the
#'   variance matrix. Default `1e-6`.
#' @param bw_grid Bandwidth grid for the adaptive mean/(auto)covariance risk.
#'   Default `NULL` sets a geometric grid from the data.
#' @param kernel_name Kernel name. Default `"epanechnikov"`.
#' @param homoscedastic If `TRUE` (default) a constant noise variance (median of
#'   the pointwise estimates) is used; otherwise the t-varying estimates.
#' @param density_bw Optional fixed design-density bandwidth reused for every
#'   `estimate_density` call (independent design only). Default `NULL` selects it
#'   once via [get_density_optimal_bw()].
#' @param sub_grid_length Number of points per axis of the coarse sub-grid on
#'   which the adaptive bandwidths are selected. Default `10`.
#'
#' @return An object of class `blup_fit`: a list whose main elements are:
#'   \itemize{
#'     \item `Tn0`, `Yn0`: the conditioning-curve design points and values.
#'     \item `rho`, `root_Dn0`: the design weights and their square-root matrix.
#'     \item `muhat_Tn0`, `c0hat`, `sigma2`: the mean, covariance operator and
#'       noise level of the conditioning curve.
#'     \item `V`, `resid`: the regularised variance matrix and the conditioning
#'       residual.
#'     \item `opt_mean`, `opt_cov`, `opt_autocov`: the cached adaptive bandwidths.
#'     \item `data`, `kernel_name`, `is_common_design`, `density_bw`, `bw_grid`,
#'       `tikhonov`, `homoscedastic`: the information needed by
#'       [predict.blup_fit()].
#'   }
#'
#' @seealso [predict.blup_fit()], [get_density_optimal_bw()].
#' @export
#' @import data.table
#' @importFrom methods is
blup_fit <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     id_lag = NULL, tikhonov = 1e-6,
                     bw_grid = NULL, kernel_name = "epanechnikov",
                     homoscedastic = TRUE, density_bw = NULL,
                     sub_grid_length = 10L) {

  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform"))

  n0 <- if (is.null(id_lag)) data[, max(id_curve)] else as.integer(id_lag)
  Tn0 <- data[id_curve == n0, sort(unique(tobs))]
  Yn0 <- data[id_curve == n0][order(tobs), X]
  Mn0 <- length(Tn0)

  is_common_design <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")

  if (is_common_design) {
    rho <- rep(1 / Mn0, Mn0)
  } else {
    if (is.null(density_bw))
      density_bw <- get_density_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        kernel_name = kernel_name, lower = 0, upper = 1)
    ghat <- estimate_density(x = Tn0, h = density_bw, kernel_name = kernel_name,
                             lower = 0, upper = 1)$estimate
    rho <- 1 / (Mn0 * pmax(ghat, 1e-6))
    rho <- rho / sum(rho)
  }
  root_Dn0 <- diag(sqrt(rho))

  if (is.null(bw_grid)) {
    N <- data[, length(unique(id_curve))]
    lambdahat <- data[, .(Mn = .N), by = id_curve][, mean(Mn)]
    K <- 15
    b0 <- ifelse(is_common_design, 0.5 / lambdahat, 0.5 / ((N * lambdahat) ** (1 / 2)))
    bK <- 0.05
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
  }

  # opt_mean/opt_cov/opt_autocov are reused by predict and by
  # select_tikhonov_parameter through blup_*_at_cpp.
  cpp <- blup_fit_cpp(
    data = data, id_lag = as.integer(n0), bw_grid = as.numeric(bw_grid),
    rho = rho, homoscedastic = homoscedastic,
    tikhonov = tikhonov, sub_grid_length = as.integer(sub_grid_length),
    kernel_name = kernel_name)

  return(structure(
    list(
      data = data,
      kernel_name = kernel_name,
      is_common_design = is_common_design,
      id_lag = n0,
      Tn0 = as.vector(cpp$Tn0),
      Yn0 = as.vector(cpp$Yn0),
      Mn0 = Mn0,
      rho = rho,
      root_Dn0 = root_Dn0,
      density_bw = density_bw,
      bw_grid = bw_grid,
      opt_mean = cpp$opt_mean,
      opt_cov = cpp$opt_cov,
      opt_autocov = cpp$opt_autocov,
      muhat_Tn0 = as.vector(cpp$muhat_Tn0),
      c0hat = cpp$c0hat,
      sigma2 = if (homoscedastic) as.numeric(cpp$sigma2) else as.vector(cpp$sigma2),
      homoscedastic = homoscedastic,
      tikhonov = tikhonov,
      V = cpp$V,
      resid = cpp$resid
    ),
    class = "blup_fit"
  ))
}

#' Predict with an adaptive functional BLUP fit
#'
#' Evaluates the adaptive Best Linear Unbiased Predictor of the curve following
#' the conditioning curve at the requested prediction points. For `horizon > 1`
#' (multi-step-ahead), each intermediate curve is predicted on the target grid
#' `t` and fed back as the new conditioning curve, so the `(n_0 + i)`-th
#' prediction is built from the `(n_0 + i - 1)`-th predicted curve, and all
#' intermediate horizons are returned.
#'
#' @details
#' During multi-step prediction the estimation data is left unchanged: the mean,
#' the (auto)covariance operators, the adaptive bandwidths and the noise level
#' are estimated once and held fixed, and each predicted curve enters only as the
#' conditioning values of the next step. Injecting a predicted (denoised) curve
#' back into the estimation sample would bias those plug-in estimates, so it is
#' deliberately avoided.
#'
#' @param object A `blup_fit` object.
#' @param t Numeric vector of prediction points in \eqn{[0, 1]}. Default is the
#'   conditioning-curve design points.
#' @param newdata Optional numeric vector of conditioning-curve values at the
#'   fit's design points (`object$Tn0`), overriding `object$Yn0`. Used internally
#'   for the multi-step recursion; must have length `object$Mn0`.
#' @param horizon Integer prediction horizon (steps ahead). Default `1`. For
#'   `horizon > 1`, each intermediate curve is predicted on the target grid `t`
#'   and fed back as the new conditioning curve; the conditioning quantities are
#'   recomputed for that grid using the cached adaptive bandwidths (works for
#'   both designs).
#' @param ... Unused; for S3 compatibility.
#'
#' @return A `data.table` with one block of rows per horizon (`horizon * length(t)`
#'   rows in total):
#'   \itemize{
#'     \item `horizon`: the prediction step, from `1` to `horizon`.
#'     \item `t`: the prediction points.
#'     \item `muhat`: the mean estimate at `t`.
#'     \item `prediction`: the adaptive BLUP at `t` for that horizon.
#'   }
#'
#' @seealso [blup_fit()].
#' @export
#' @import data.table
predict.blup_fit <- function(object, t = object$Tn0, newdata = NULL, horizon = 1L, ...) {
  if (!(methods::is(t, "numeric") && all(t >= 0 & t <= 1)))
    stop("'t' must be a numeric vector with values between 0 and 1.")
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("'horizon' must be a positive integer.")

  Yn0 <- if (is.null(newdata)) object$Yn0 else as.numeric(newdata)
  if (length(Yn0) != object$Mn0)
    stop("'newdata' must have length equal to the fit's design (", object$Mn0, ").")

  # density_bw is unused under the common design; pass a dummy positive value.
  density_bw <- if (is.null(object$density_bw)) 0.1 else object$density_bw

  out <- blup_predict_cpp(
    data = object$data, opt_mean = object$opt_mean, opt_cov = object$opt_cov,
    opt_autocov = object$opt_autocov, Tn0 = as.numeric(object$Tn0),
    muhat_Tn0 = as.numeric(object$muhat_Tn0), V = object$V, root_D = object$root_Dn0,
    Yn0 = as.numeric(Yn0), density_bw = density_bw,
    is_common = object$is_common_design, homoscedastic = object$homoscedastic,
    tikhonov = object$tikhonov, t = as.numeric(t), horizon = horizon,
    kernel_name = object$kernel_name)

  return(data.table::data.table(
    horizon = as.integer(out[, 1]), t = out[, 2],
    muhat = out[, 3], prediction = out[, 4]))
}

#' Fit and predict the adaptive functional BLUP in one call
#'
#' Convenience wrapper that fits the adaptive BLUP on `data` and immediately
#' predicts the curve following the conditioning curve at `t`. Equivalent to
#' `predict(blup_fit(data, ...), t = t, h = h)`.
#'
#' @inheritParams blup_fit
#' @param t Numeric vector of prediction points in \eqn{[0, 1]}.
#' @param horizon Integer prediction horizon (steps ahead). Default `1`.
#'
#' @return A `data.table` with columns `horizon`, `t`, `muhat` and `prediction`
#'   (one block of rows per horizon); see [predict.blup_fit()].
#'
#' @seealso [blup_fit()], [predict.blup_fit()], [select_tikhonov_parameter()].
#' @export
#' @import data.table
#' @importFrom stats predict
blup <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                 t = seq(0.01, 0.99, length.out = 99), id_lag = NULL, horizon = 1L,
                 tikhonov = 1e-6, bw_grid = NULL,
                 kernel_name = "epanechnikov", homoscedastic = TRUE,
                 density_bw = NULL, sub_grid_length = 10L) {
  fit <- blup_fit(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol, id_lag = id_lag,
    tikhonov = tikhonov, bw_grid = bw_grid,
    kernel_name = kernel_name, homoscedastic = homoscedastic,
    density_bw = density_bw, sub_grid_length = sub_grid_length)
  return(predict(fit, t = t, horizon = horizon))
}

#' Select the Tikhonov regularisation parameter
#'
#' Selects the Tikhonov regularisation parameter \eqn{\alpha} of the adaptive
#' BLUP. Currently only `method = "cv"` is implemented, a one-step-ahead
#' cross-validation.
#'
#' @details
#' Each of the last `n_val` curves is predicted from its immediate predecessor
#' and scored by the design-weighted squared prediction error at its observation
#' points, \eqn{\sum_i \varrho_{n,i}\,(Y_{n,i} - \widehat X_n(T_{n,i};\alpha))^2},
#' where \eqn{\varrho_{n,i}} is the design weight of the held-out (target) curve.
#' Two regimes:
#' \itemize{
#'   \item \strong{Common design} — the operators are estimated once on the
#'     training block and only the conditioning values vary across the
#'     validation set.
#'   \item \strong{Independent design} — a rolling origin: for each target curve
#'     the plug-in estimates are refreshed on all curves observed up to its
#'     predecessor, with the adaptive bandwidths held fixed at those selected
#'     once on the initial block (cached in the fit).
#' }
#' Only the \eqn{M \times M} system depends on \eqn{\alpha}; it is solved for the
#' whole grid from a single eigendecomposition per fold.
#'
#' @inheritParams blup_fit
#' @param method Selection method. Currently only `"cv"` (one-step-ahead
#'   cross-validation) is available.
#' @param tikhonov_grid Candidate values. If `NULL`, a 25-point grid
#'   \eqn{\{e^{-5}, \ldots, e^3\}} is used.
#' @param n_val Number of trailing curves used for one-step-ahead validation.
#'
#' @return A list with:
#'   \itemize{
#'     \item `tikhonov_star`: the selected Tikhonov parameter.
#'     \item `tikhonov_grid`: the candidate grid.
#'     \item `cv_curve`: the mean cross-validation score per candidate.
#'     \item `cv_matrix`: the per-fold cross-validation scores (fold by candidate).
#'     \item `val_ids`: the ids of the validation curves.
#'   }
#'
#' @seealso [blup_fit()], [predict.blup_fit()].
#' @export
#' @import data.table
select_tikhonov_parameter <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                      method = c("cv"), tikhonov_grid = NULL, n_val = 30L,
                                      bw_grid = NULL, kernel_name = "epanechnikov",
                                      homoscedastic = TRUE, density_bw = NULL) {
  method <- match.arg(method)

  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform"))

  ids <- data[, sort(unique(id_curve))]
  n <- length(ids)
  if (n_val >= n) stop("'n_val' must be smaller than the number of curves.")
  fit_ids <- ids[seq_len(n - n_val)]
  val_pos <- (n - n_val + 1L):n
  data_fit <- data[id_curve %in% fit_ids]
  is_common <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")

  if (!is_common && is.null(density_bw))
    density_bw <- get_density_optimal_bw(
      data = data_fit, idcol = "id_curve", tcol = "tobs", ycol = "X",
      kernel_name = kernel_name, lower = 0, upper = 1)

  fit <- blup_fit(
    data = data_fit, id_lag = max(fit_ids), tikhonov = 1e-6,
    bw_grid = bw_grid, kernel_name = kernel_name, homoscedastic = homoscedastic,
    density_bw = density_bw)

  # Common design: the operators are constant across folds.
  if (is_common) {
    Tn0 <- fit$Tn0
    c1_common <- blup_autocov_at_cpp(data_fit, fit$opt_autocov, Tn0, Tn0, 1L, FALSE, kernel_name)
    A0_common <- fit$root_Dn0 %*% fit$c0hat %*% fit$root_Dn0 + diag(fit$sigma2 * fit$rho)
    C1rD_common <- t(c1_common) %*% fit$root_Dn0
  }

  ## Phase 1: assemble the Tikhonov-free pieces for each validation curve.
  folds <- vector("list", n_val)
  for (k in seq_len(n_val)) {
    id_targ <- ids[val_pos[k]]
    id_prev <- ids[val_pos[k] - 1L]
    Y_prev <- data[id_curve == id_prev][order(tobs), X]
    Y_targ <- data[id_curve == id_targ][order(tobs), X]

    if (is_common) {
      folds[[k]] <- list(
        A0 = A0_common, C1rD = C1rD_common, mu_pred = fit$muhat_Tn0,
        resid = fit$root_Dn0 %*% matrix(Y_prev - fit$muhat_Tn0, ncol = 1),
        Y_targ = Y_targ, rho_targ = rep(1 / length(Y_targ), length(Y_targ)))
    } else {
      # Rolling origin: refresh the plug-in estimates on the grown window with
      # the bandwidths held fixed at those cached in `fit`.
      data_roll <- data[id_curve <= id_prev]
      Tprev <- data[id_curve == id_prev, sort(unique(tobs))]
      Ttarg <- data[id_curve == id_targ, sort(unique(tobs))]
      c0 <- blup_autocov_at_cpp(data_roll, fit$opt_cov, Tprev, Tprev, 0L, TRUE, kernel_name)
      c0 <- psd_project_cpp(c0)
      c1 <- blup_autocov_at_cpp(data_roll, fit$opt_autocov, Tprev, Ttarg, 1L, FALSE, kernel_name)
      mu_prev <- blup_mean_at_cpp(data_roll, fit$opt_mean, Tprev, kernel_name)
      mu_targ <- blup_mean_at_cpp(data_roll, fit$opt_mean, Ttarg, kernel_name)
      sig2 <- adaptiveFTS::estimate_sigma(
        data = data_roll, idcol = "id_curve", tcol = "tobs", ycol = "X", t = Tprev)[, sig ** 2]
      if (homoscedastic) sig2 <- stats::median(sig2, na.rm = TRUE)
      ghat <- estimate_density(x = Tprev, h = density_bw, kernel_name = kernel_name,
                               lower = 0, upper = 1)$estimate
      rho <- 1 / (length(Tprev) * pmax(ghat, 1e-6))
      rho <- rho / sum(rho)
      root_D <- diag(sqrt(rho))
      ghat_targ <- estimate_density(x = Ttarg, h = density_bw, kernel_name = kernel_name,
                                    lower = 0, upper = 1)$estimate
      rho_targ <- 1 / (length(Ttarg) * pmax(ghat_targ, 1e-6))
      rho_targ <- rho_targ / sum(rho_targ)
      folds[[k]] <- list(
        A0 = root_D %*% c0 %*% root_D + diag(sig2 * rho), C1rD = t(c1) %*% root_D,
        mu_pred = mu_targ, resid = root_D %*% matrix(Y_prev - mu_prev, ncol = 1),
        Y_targ = Y_targ, rho_targ = rho_targ)
    }
  }

  if (is.null(tikhonov_grid)) tikhonov_grid <- exp(seq(-5, 3, length.out = 25))

  ## Phase 2: only the (A0 + tikhonov * I)^{-1} step depends on the Tikhonov
  ## parameter; A0 is symmetric, so eigendecompose once per fold and reuse it
  ## across the whole grid.
  cv_matrix <- matrix(NA_real_, nrow = n_val, ncol = length(tikhonov_grid))
  for (k in seq_len(n_val)) {
    f <- folds[[k]]
    eg <- eigen(f$A0, symmetric = TRUE)
    lambda <- eg$values
    z <- as.vector(crossprod(eg$vectors, f$resid))
    W <- f$C1rD %*% eg$vectors
    for (l in seq_along(tikhonov_grid)) {
      pred <- tryCatch(
        f$mu_pred + as.vector(W %*% (z / (lambda + tikhonov_grid[l]))),
        error = function(e) rep(NA_real_, length(f$Y_targ)))
      cv_matrix[k, l] <- sum(f$rho_targ * (f$Y_targ - pred) ^ 2)
    }
  }

  cv_curve <- colMeans(cv_matrix, na.rm = TRUE)
  l_star <- which.min(cv_curve)
  tikhonov_star <- tikhonov_grid[l_star]
  if (l_star %in% c(1L, length(tikhonov_grid)))
    warning("tikhonov_star at a grid boundary; widen tikhonov_grid.")

  return(list(
    tikhonov_star = tikhonov_star,
    tikhonov_grid = tikhonov_grid,
    cv_curve = cv_curve,
    cv_matrix = cv_matrix,
    val_ids = ids[val_pos]
  ))
}
