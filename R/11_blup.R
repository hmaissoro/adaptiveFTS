## =============================================================================
## Adaptive functional BLUP: fit / predict engine
##
## `blup_fit()` estimates every component that does not depend on the
## prediction points (the design-weighted covariance operator of the
## conditioning curve, its mean, the noise level, the design weights, the
## adaptive bandwidths and the regularised variance matrix) and returns a
## `blup_fit` object. `predict.blup_fit()` then evaluates the one-step-ahead
## adaptive BLUP at arbitrary prediction points, and loops for h-step-ahead
## prediction by feeding each predicted curve back as the new conditioning
## curve. The split mirrors `stats::lm()` / `stats::predict()`.
##
## The method is the design-weighted, Tikhonov-regularised adaptive BLUP of
## Maissoro, Patilea and Vimond: it conditions on the previous curve only
## (single lag-1 block). The internal covariance assembly is kept separable so
## the multi-lag generalisation can be added later without a rewrite.
## =============================================================================

#' Check whether all curves share the same observation design
#'
#' @param data A prepared functional data.table.
#' @param idcol,tcol Identifier and observation-time column names.
#' @return `TRUE` if every curve shares the same sorted design, else `FALSE`.
#' @keywords internal
.is_common_design <- function(data, idcol = "id_curve", tcol = "tobs") {
  design_by_curve <- data[, .(tobs_list = list(sort(unique(get(tcol))))), by = idcol]
  reference_design <- design_by_curve$tobs_list[[1]]
  all(vapply(design_by_curve$tobs_list, function(tt) identical(tt, reference_design), logical(1)))
}

#' Nearest-neighbour index in one dimension
#'
#' Returns, for each query point, the index of the nearest reference point.
#' Base-R replacement for `RANN::nn2(..., k = 1)` in 1-D; ties are broken by the
#' first (smallest) index, matching `which.min`.
#'
#' @param ref Numeric vector of reference locations.
#' @param query Numeric vector of query locations.
#' @return Integer vector of nearest-reference indices, one per query.
#' @keywords internal
.nn1_1d <- function(ref, query) {
  vapply(query, function(q) which.min(abs(ref - q)), integer(1))
}

#' Nearest-neighbour index in two dimensions
#'
#' Returns, for each query row, the index of the nearest reference row under the
#' Euclidean metric. Base-R replacement for `RANN::nn2(..., k = 1)` in 2-D.
#'
#' @param ref_s,ref_t Numeric vectors of reference coordinates.
#' @param q_s,q_t Numeric vectors of query coordinates.
#' @return Integer vector of nearest-reference indices, one per query.
#' @keywords internal
.nn1_2d <- function(ref_s, ref_t, q_s, q_t) {
  vapply(seq_along(q_s), function(i) {
    which.min((ref_s - q_s[i]) ^ 2 + (ref_t - q_t[i]) ^ 2)
  }, integer(1))
}

#' Adaptive optimal mean bandwidths on a sub-grid
#'
#' Runs `estimate_mean_risk` on a coarse sub-grid and returns, for each sub-grid
#' point, the bandwidth minimising the risk. Cached in the fit and reused at
#' predict time via nearest-neighbour matching.
#' @keywords internal
.mean_optbw <- function(data, sub_grid_vec, bw_grid, kernel_name) {
  dt_risk_mean <- adaptiveFTS::estimate_mean_risk(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = sub_grid_vec, bw_grid = bw_grid, kernel_name = kernel_name)
  dt_risk_mean[, .("optbw" = h[which.min(mean_risk)]), by = t]
}

#' Adaptive optimal (auto)covariance bandwidths on a sub-grid
#' @keywords internal
.autocov_optbw <- function(data, sub_grid, lag, bw_grid, kernel_name) {
  dt_risk <- adaptiveFTS::estimate_autocov_risk(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = sub_grid$s, t = sub_grid$t, lag = lag, bw_grid = bw_grid,
    use_same_bw = FALSE, center = TRUE, kernel_name = kernel_name)
  dt_risk[
    ,
    .("optbw_s" = hs[which.min(autocov_risk)], "optbw_t" = ht[which.min(autocov_risk)]),
    by = c("s", "t")
  ]
}

#' Evaluate the mean at new locations using cached bandwidths
#'
#' Reuses the adaptive optimal bandwidths selected in `blup_fit` (`dt_optbw_mean`),
#' matched by nearest neighbour to `t`, so only the plug-in `estimate_mean` is
#' re-run; the risk minimisation is not repeated.
#' @keywords internal
.mean_at <- function(dt_optbw_mean, data, t, kernel_name) {
  idx <- .nn1_1d(dt_optbw_mean$t, t)
  optbw_t <- dt_optbw_mean$optbw[idx]
  dt <- adaptiveFTS::estimate_mean(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t, optbw = optbw_t, bw_grid = NULL, kernel_name = kernel_name)
  dt[order(t), muhat]
}

#' Evaluate a (auto)covariance block at new locations using cached bandwidths
#'
#' Reuses the adaptive optimal bandwidths selected in `blup_fit`
#' (`dt_optbw_cov` for `lag = 0`, `dt_optbw_autocov` for `lag = 1`), matched by
#' nearest neighbour, so only the plug-in `estimate_autocov` is re-run.
#' @return A `length(s)` x `length(t)` matrix of \eqn{\hat c_{lag}(s_i, t_j)}.
#' @keywords internal
.autocov_at <- function(dt_optbw, data, s, t, lag, kernel_name) {
  grid <- data.table::as.data.table(expand.grid("s" = s, "t" = t))
  idx <- .nn1_2d(dt_optbw$s, dt_optbw$t, grid$s, grid$t)
  grid$optbw_s <- dt_optbw$optbw_s[idx]
  grid$optbw_t <- dt_optbw$optbw_t[idx]
  dt <- adaptiveFTS::estimate_autocov(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = grid[, s], t = grid[, t], lag = lag,
    optbw_s = grid[, optbw_s], optbw_t = grid[, optbw_t],
    bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
    correct_diagonal = (lag == 0L), kernel_name = kernel_name)
  dt_dcast <- data.table::dcast(dt[order(s, t)], formula = s ~ t, value.var = "autocov")
  m <- as.matrix(dt_dcast[, .SD, .SDcols = !"s"])
  colnames(m) <- NULL
  m
}

#' Fit the adaptive functional BLUP
#'
#' Estimates every component of the adaptive Best Linear Unbiased Predictor that
#' does not depend on the prediction points, conditioning on a single curve (its
#' immediate successor is what `predict()` reconstructs). The returned object
#' caches the adaptive bandwidths so that `predict()` only re-runs the cheap
#' plug-in estimates at the requested prediction points.
#'
#' @inheritParams format_data
#' @param id_lag Integer id of the conditioning curve. Its successor is the
#'   curve to be predicted. Default `NULL` uses the last curve in `data`.
#' @param tikhonov_reg_param Tikhonov regularisation parameter \eqn{\alpha}
#'   added to the variance matrix. Default `1e-6`.
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
#' @return An object of class `blup_fit`: a list with the conditioning-curve
#'   design and values, the design weights, the mean and covariance estimates,
#'   the regularised variance matrix, the cached adaptive bandwidths and the
#'   information needed by [predict.blup_fit()].
#'
#' @seealso [predict.blup_fit()], [get_density_optimal_bw()].
#' @export
#' @import data.table
#' @importFrom methods is
blup_fit <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                     id_lag = NULL, tikhonov_reg_param = 1e-6,
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

  # Design weights rho_{n0,i}
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

  # Default bandwidth grid (common-design geometric grid)
  if (is.null(bw_grid)) {
    N <- data[, length(unique(id_curve))]
    lambdahat <- data[, .(Mn = .N), by = id_curve][, mean(Mn)]
    K <- 15
    b0 <- ifelse(is_common_design, 0.5 / lambdahat, 0.5 / ((N * lambdahat) ** (1 / 2)))
    bK <- 0.05
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
  }

  sub_grid_vec <- seq(0.05, 0.95, length.out = sub_grid_length)
  sub_grid <- expand.grid("s" = sub_grid_vec, "t" = sub_grid_vec)

  # Adaptive bandwidths (cached for predict-time reuse)
  dt_optbw_mean <- .mean_optbw(data, sub_grid_vec, bw_grid, kernel_name)
  dt_optbw_cov <- .autocov_optbw(data, sub_grid, lag = 0, bw_grid, kernel_name)
  dt_optbw_autocov <- .autocov_optbw(data, sub_grid, lag = 1, bw_grid, kernel_name)

  # Mean at the conditioning-curve design points
  muhat_Tn0 <- .mean_at(dt_optbw_mean, data, Tn0, kernel_name)

  # Covariance C0 at Tn0 x Tn0 (symmetrised)
  c0hat <- .autocov_at(dt_optbw_cov, data, Tn0, Tn0, lag = 0, kernel_name)
  c0hat <- (c0hat + t(c0hat)) / 2

  # Noise variance at Tn0
  sigma2 <- adaptiveFTS::estimate_sigma(
    data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = Tn0)[, sig ** 2]
  if (homoscedastic) sigma2 <- stats::median(sigma2, na.rm = TRUE)

  # Regularised variance matrix and the alpha-independent conditioning residual
  V <- root_Dn0 %*% c0hat %*% root_Dn0 + diag(sigma2 * rho) + tikhonov_reg_param * diag(Mn0)
  resid <- root_Dn0 %*% matrix(data = Yn0 - muhat_Tn0, ncol = 1)

  structure(
    list(
      data = data,
      kernel_name = kernel_name,
      is_common_design = is_common_design,
      id_lag = n0,
      Tn0 = Tn0,
      Yn0 = Yn0,
      Mn0 = Mn0,
      rho = rho,
      root_Dn0 = root_Dn0,
      density_bw = density_bw,
      bw_grid = bw_grid,
      dt_optbw_mean = dt_optbw_mean,
      dt_optbw_cov = dt_optbw_cov,
      dt_optbw_autocov = dt_optbw_autocov,
      muhat_Tn0 = muhat_Tn0,
      c0hat = c0hat,
      sigma2 = sigma2,
      homoscedastic = homoscedastic,
      tikhonov_reg_param = tikhonov_reg_param,
      V = V,
      resid = resid
    ),
    class = "blup_fit"
  )
}

#' Predict with an adaptive functional BLUP fit
#'
#' Evaluates the adaptive Best Linear Unbiased Predictor of the curve following
#' the conditioning curve at the requested prediction points. For `h > 1`
#' (h-step-ahead), each intermediate curve is predicted at the (common) design
#' points and fed back as the new conditioning curve, so the `(n_0 + i)`-th
#' prediction is built from the `(n_0 + i - 1)`-th predicted curve.
#'
#' @param object A `blup_fit` object.
#' @param t Numeric vector of prediction points in \eqn{[0, 1]}. Default is the
#'   conditioning-curve design points.
#' @param newdata Optional numeric vector of conditioning-curve values at the
#'   fit's design points (`object$Tn0`), overriding `object$Yn0`. Used internally
#'   for the h-step recursion; must have length `object$Mn0`.
#' @param h Integer prediction horizon (steps ahead). Default `1`. For `h > 1`,
#'   each intermediate curve is predicted on the target grid `t` and fed back as
#'   the new conditioning curve; the conditioning quantities are recomputed for
#'   that grid using the cached adaptive bandwidths (works for both designs).
#' @param ... Unused; for S3 compatibility.
#'
#' @return A `data.table` with columns `t`, `muhat` (mean estimate) and
#'   `prediction` (the adaptive BLUP).
#'
#' @seealso [blup_fit()].
#' @export
#' @import data.table
predict.blup_fit <- function(object, t = object$Tn0, newdata = NULL, h = 1L, ...) {
  if (!(methods::is(t, "numeric") && all(t >= 0 & t <= 1)))
    stop("'t' must be a numeric vector with values between 0 and 1.")
  h <- as.integer(h)
  if (h < 1L) stop("'h' must be a positive integer.")

  # Conditioning-curve values (possibly overridden for the recursion)
  Yn0 <- if (is.null(newdata)) object$Yn0 else as.numeric(newdata)
  if (length(Yn0) != object$Mn0)
    stop("'newdata' must have length equal to the fit's design (", object$Mn0, ").")

  data <- object$data
  kern <- object$kernel_name
  alpha <- object$tikhonov_reg_param

  # Conditioning-dependent quantities for an arbitrary design `Td` and values
  # `Yd`, reusing the cached adaptive bandwidths and the historical data for the
  # plug-in estimates. Used when a predicted curve (on grid `t`) becomes the new
  # conditioning curve during the h-step recursion.
  condition <- function(Td, Yd) {
    Md <- length(Td)
    if (object$is_common_design) {
      rho <- rep(1 / Md, Md)
    } else {
      ghat <- estimate_density(x = Td, h = object$density_bw, kernel_name = kern,
                               lower = 0, upper = 1)$estimate
      rho <- 1 / (Md * pmax(ghat, 1e-6))
      rho <- rho / sum(rho)
    }
    root_D <- diag(sqrt(rho))
    muhat_Td <- .mean_at(object$dt_optbw_mean, data, Td, kern)
    c0 <- .autocov_at(object$dt_optbw_cov, data, Td, Td, lag = 0, kern)
    c0 <- (c0 + t(c0)) / 2
    sig2 <- adaptiveFTS::estimate_sigma(
      data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = Td)[, sig ** 2]
    if (object$homoscedastic) sig2 <- stats::median(sig2, na.rm = TRUE)
    V <- root_D %*% c0 %*% root_D + diag(sig2 * rho) + alpha * diag(Md)
    resid <- root_D %*% matrix(Yd - muhat_Td, ncol = 1)
    list(Td = Td, root_D = root_D, V = V, resid = resid)
  }

  # One-step BLUP at `tpred` given conditioning quantities `cond`.
  blup_at <- function(cond, tpred) {
    muhat_t <- .mean_at(object$dt_optbw_mean, data, tpred, kern)
    c1 <- .autocov_at(object$dt_optbw_autocov, data, cond$Td, tpred, lag = 1, kern)
    pred <- muhat_t + t(c1) %*% cond$root_D %*% solve(cond$V, cond$resid)
    list(muhat = muhat_t, prediction = as.vector(pred))
  }

  # Step 1 conditions on the fitted (real) curve: reuse the cached V and root_D
  # so the h = 1 result is bit-identical to blup_fit's assembly.
  cond <- list(Td = object$Tn0, root_D = object$root_Dn0, V = object$V,
               resid = object$root_Dn0 %*% matrix(Yn0 - object$muhat_Tn0, ncol = 1))

  if (h == 1L) {
    res <- blup_at(cond, t)
  } else {
    xprev <- blup_at(cond, t)$prediction
    if (h > 2L)
      for (step in seq_len(h - 2L)) {
        cond <- condition(t, xprev)
        xprev <- blup_at(cond, t)$prediction
      }
    cond <- condition(t, xprev)
    res <- blup_at(cond, t)
  }

  data.table::data.table(t = t, muhat = res$muhat, prediction = res$prediction)
}

#' Fit and predict the adaptive functional BLUP in one call
#'
#' Convenience wrapper that fits the adaptive BLUP on `data` and immediately
#' predicts the curve following the conditioning curve at `t`. Equivalent to
#' `predict(blup_fit(data, ...), t = t, h = h)`.
#'
#' @inheritParams blup_fit
#' @param t Numeric vector of prediction points in \eqn{[0, 1]}.
#' @param h Integer prediction horizon (steps ahead). Default `1`.
#'
#' @return A `data.table` with columns `t`, `muhat` and `prediction`.
#'
#' @seealso [blup_fit()], [predict.blup_fit()], [cv_blup_alpha()].
#' @export
#' @import data.table
blup <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                 t = seq(0.01, 0.99, length.out = 99), id_lag = NULL, h = 1L,
                 tikhonov_reg_param = 1e-6, bw_grid = NULL,
                 kernel_name = "epanechnikov", homoscedastic = TRUE,
                 density_bw = NULL, sub_grid_length = 10L) {
  fit <- blup_fit(
    data = data, idcol = idcol, tcol = tcol, ycol = ycol, id_lag = id_lag,
    tikhonov_reg_param = tikhonov_reg_param, bw_grid = bw_grid,
    kernel_name = kernel_name, homoscedastic = homoscedastic,
    density_bw = density_bw, sub_grid_length = sub_grid_length)
  predict(fit, t = t, h = h)
}

#' One-step-ahead cross-validation for the Tikhonov parameter
#'
#' Selects the Tikhonov regularisation parameter \eqn{\alpha} by one-step-ahead
#' cross-validation: each of the last `n_val` curves is predicted from its
#' immediate predecessor and scored by the design-weighted squared prediction
#' error at its observation points,
#' \eqn{\sum_i \varrho_{n,i}\,(Y_{n,i} - \widehat X_n(T_{n,i};\alpha))^2}, where
#' \eqn{\varrho_{n,i}} is the design weight of the held-out (target) curve. Two
#' regimes:
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
#' @param alpha_grid Candidate values. If `NULL`, a 25-point grid
#'   \eqn{\{e^{-5}, \ldots, e^0\}} is used.
#' @param n_val Number of trailing curves used for one-step-ahead validation.
#'
#' @return A list with `alpha_star`, `alpha_grid`, `cv_curve`, `cv_matrix`
#'   (fold by alpha) and `val_ids`.
#'
#' @seealso [blup_fit()].
#' @export
#' @import data.table
cv_blup_alpha <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          alpha_grid = NULL, n_val = 30L, bw_grid = NULL,
                          kernel_name = "epanechnikov", homoscedastic = TRUE,
                          density_bw = NULL) {

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

  # Select the design-density bandwidth once on the initial training block.
  if (!is_common && is.null(density_bw))
    density_bw <- get_density_optimal_bw(
      data = data_fit, idcol = "id_curve", tcol = "tobs", ycol = "X",
      kernel_name = kernel_name, lower = 0, upper = 1)

  # Estimate every alpha-free component once on the training block.
  fit <- blup_fit(
    data = data_fit, id_lag = max(fit_ids), tikhonov_reg_param = 1e-6,
    bw_grid = bw_grid, kernel_name = kernel_name, homoscedastic = homoscedastic,
    density_bw = density_bw)

  # Common-design invariants (the operators do not change across folds).
  if (is_common) {
    Tn0 <- fit$Tn0
    c1_common <- .autocov_at(fit$dt_optbw_autocov, data_fit, Tn0, Tn0, lag = 1, kernel_name)
    A0_common <- fit$root_Dn0 %*% fit$c0hat %*% fit$root_Dn0 + diag(fit$sigma2 * fit$rho)
    C1rD_common <- t(c1_common) %*% fit$root_Dn0
  }

  ## Phase 1: assemble the alpha-free pieces for each validation curve.
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
      c0 <- .autocov_at(fit$dt_optbw_cov, data_roll, Tprev, Tprev, lag = 0, kernel_name)
      c0 <- (c0 + t(c0)) / 2
      c1 <- .autocov_at(fit$dt_optbw_autocov, data_roll, Tprev, Ttarg, lag = 1, kernel_name)
      mu_prev <- .mean_at(fit$dt_optbw_mean, data_roll, Tprev, kernel_name)
      mu_targ <- .mean_at(fit$dt_optbw_mean, data_roll, Ttarg, kernel_name)
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

  if (is.null(alpha_grid)) alpha_grid <- exp(seq(-5, 0, length.out = 25))

  ## Phase 2: only (A0 + alpha I)^{-1} depends on alpha; A0 is symmetric, so
  ## eigendecompose once per fold and reuse across the whole grid.
  cv_matrix <- matrix(NA_real_, nrow = n_val, ncol = length(alpha_grid))
  for (k in seq_len(n_val)) {
    f <- folds[[k]]
    eg <- eigen(f$A0, symmetric = TRUE)
    lambda <- eg$values
    z <- as.vector(crossprod(eg$vectors, f$resid))
    W <- f$C1rD %*% eg$vectors
    for (l in seq_along(alpha_grid)) {
      pred <- tryCatch(
        f$mu_pred + as.vector(W %*% (z / (lambda + alpha_grid[l]))),
        error = function(e) rep(NA_real_, length(f$Y_targ)))
      cv_matrix[k, l] <- sum(f$rho_targ * (f$Y_targ - pred) ^ 2)
    }
  }

  cv_curve <- colMeans(cv_matrix, na.rm = TRUE)
  l_star <- which.min(cv_curve)
  alpha_star <- alpha_grid[l_star]
  if (l_star %in% c(1L, length(alpha_grid)))
    warning("alpha_star at a grid boundary; widen alpha_grid.")

  list(
    alpha_star = alpha_star,
    alpha_grid = alpha_grid,
    cv_curve = cv_curve,
    cv_matrix = cv_matrix,
    val_ids = ids[val_pos]
  )
}
