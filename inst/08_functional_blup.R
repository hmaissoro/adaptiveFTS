library(data.table)
library(ggplot2)
Rcpp::sourceCpp("./src/08_estimate_curve_cpp.cpp")
source("./R/10_density_estimator.R")


#' Check whether all curves share the same observation design
#'
#' Internal helper used to verify that every curve in a data table has the
#' same sorted set of observation times.
#'
#' @param data A data.table containing the functional data.
#' @param idcol Character scalar giving the curve identifier column.
#' @param tcol Character scalar giving the observation time column.
#'
#' @return A logical scalar. `TRUE` if all curves share the same design,
#'   otherwise `FALSE`.
#' @keywords internal
.is_common_design <- function(data, idcol = "id_curve", tcol = "tobs") {
    design_by_curve <- data[, .(tobs_list = list(sort(unique(get(tcol))))), by = idcol]
    reference_design <- design_by_curve$tobs_list[[1]]
    all(vapply(design_by_curve$tobs_list, function(tt) identical(tt, reference_design), logical(1)))
}

# Import the data
data("data_far")


# Prepare the data
data_prepared <- adaptiveFTS::format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")
data_train <- data_prepared[-which(id_curve == 150)]
data_test <- data_prepared[which(id_curve == 150)]
rm(data_prepared)

# Prediction points of the new curve X_{n_0 + 1}
t0 <- data_test[, sort(unique(tobs))]

# Set bandwidth grid for adaptive estimation of the mean and (auto)covariance functions
is_common_design <- .is_common_design(data = data_train, idcol = "id_curve", tcol = "tobs")
N <- data_train[, length(unique(id_curve))]
lambdahat <- data_train[, .(Mn = .N), by = id_curve][, mean(Mn)]
K <- 15
b0 <- ifelse(is_common_design, 0.5 / lambdahat, 0.5 / ((N * lambdahat) ** (1/2)))
bK <- 0.05
a <- exp((log(bK) - log(b0)) / K)
bw_grid_blup <- b0 * a ** (seq_len(K))

# Adaptive BLUP for the new curve X_{n_0 + 1}
predict_next_curve <- function(
    data = data_train, prediction_points = t0, tikhonov_reg_param = 1e-6,
    bw_grid = bw_grid_blup, kernel_name = "epanechnikov",
    homoscedastic = TRUE, density_bw = NULL) {

    n0 <- data[, max(id_curve)]
    Tn0 <- data[id_curve == n0, sort(unique(tobs))]
    Yn0 <- data[id_curve == n0][order(tobs), X]
    Mn0 <- length(Tn0)

    # Compute the weights \rho_{n0,i} for i = 1, ..., M_{n0}
    is_common_design <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")

    if (is_common_design) {
        rho <- rep(1 / Mn0, Mn0)
    } else {
        # Select the design-density bandwidth once on a subset of curves, then
        # reuse the exact same value for every estimate_density() call.
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

    # Set bandwidth grid : recall that we are in the common design case
    if (is.null(bw_grid)) {
        N <- data[, length(unique(id_curve))]
        lambdahat <- data[, .(Mn = .N), by = id_curve][, mean(Mn)]
        K <- 15
        b0 <- ifelse(is_common_design, 0.5 / lambdahat, 0.5 / ((N * lambdahat) ** (1/2)))
        bK <- 0.05
        a <- exp((log(bK) - log(b0)) / K)
        bw_grid <- b0 * a ** (seq_len(K))
    }

    # Build estimation locations for the covariance and autocovariance functions
    grid_cov <- expand.grid("s" = Tn0, "t" = Tn0)
    grid_cov <- data.table::as.data.table(grid_cov)
    grid_autocov <- expand.grid("s" = Tn0, "t" = prediction_points)
    grid_autocov <- data.table::as.data.table(grid_autocov)

    ## To reduce computation time, we estimate the (auto)covariance function 
    ## on a sub-grid and then use nearest neighbor matching to find the optimal 
    ## bandwidths for the full grid.
    sub_grid_vec <- seq(0.05, 0.95, length.out = 10)
    sub_grid <- expand.grid("s" = sub_grid_vec, "t" = sub_grid_vec)

    # Estimate the mean function at the prediction points
    dt_risk_mean <- adaptiveFTS::estimate_mean_risk(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = sub_grid_vec, bw_grid = bw_grid, kernel_name = kernel_name)
    dt_optbw_mean <- dt_risk_mean[
        ,
        .("optbw" = h[which.min(mean_risk)]),
        by = t
    ]

    knn_mean_pred <- RANN::nn2(
        data = matrix(dt_optbw_mean$t, ncol = 1),
        query = matrix(prediction_points, ncol = 1), k = 1
    )
    optbw_mean_prediction_points <- dt_optbw_mean$optbw[knn_mean_pred$nn.idx]
    dt_muhat_prediction_points <- adaptiveFTS::estimate_mean(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = prediction_points, optbw = optbw_mean_prediction_points, bw_grid = bw_grid,
        kernel_name = kernel_name)
    muhat_prediction_points <- dt_muhat_prediction_points[order(t), muhat]

    # Estimate the mean function at Tn0 points
    knn_Tn0 <- RANN::nn2(
        data = matrix(dt_optbw_mean$t, ncol = 1),
        query = matrix(Tn0, ncol = 1), k = 1
    )
    optbw_Tn0 <- dt_optbw_mean$optbw[knn_Tn0$nn.idx]
    dt_muhat_Tn0 <- adaptiveFTS::estimate_mean(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = Tn0, optbw = optbw_Tn0, bw_grid = bw_grid,
        kernel_name = kernel_name)
    muhat_Tn0 <- dt_muhat_Tn0[order(t), muhat]

    ## Estimate coavariance risk function
    dt_risk_cov <- adaptiveFTS::estimate_autocov_risk(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = sub_grid$s, t = sub_grid$t, lag = 0, bw_grid = bw_grid,
        use_same_bw = FALSE, center = TRUE, kernel_name = kernel_name)
    dt_optbw_cov <- dt_risk_cov[
        ,
        .("optbw_s" = hs[which.min(autocov_risk)], "optbw_t" = ht[which.min(autocov_risk)]),
        by = c("s", "t")
    ]

    ## Estimate lag-1 autocovariance risk function
    dt_risk_autocov <- adaptiveFTS::estimate_autocov_risk(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = sub_grid$s, t = sub_grid$t, lag = 1, bw_grid = bw_grid,
        use_same_bw = FALSE, center = TRUE, kernel_name = kernel_name)
    dt_optbw_autocov <- dt_risk_autocov[
        ,
        .("optbw_s" = hs[which.min(autocov_risk)], "optbw_t" = ht[which.min(autocov_risk)]),
        by = c("s", "t")
    ]

    

    ## Matching (auto)cov grid using RANN to find nearest neighbors
    ## Covariance
    knn_cov <- RANN::nn2(data = dt_optbw_cov[, 1:2], query = grid_cov[, 1:2], k = 1)
    grid_cov$optbw_s <- dt_optbw_cov$optbw_s[knn_cov$nn.idx]
    grid_cov$optbw_t <- dt_optbw_cov$optbw_t[knn_cov$nn.idx]

    ## lag-1 Autocovariance
    knn_autocov <- RANN::nn2(data = dt_optbw_autocov[, 1:2], query = grid_autocov[, 1:2], k = 1)
    grid_autocov$optbw_s <- dt_optbw_autocov$optbw_s[knn_autocov$nn.idx]
    grid_autocov$optbw_t <- dt_optbw_autocov$optbw_t[knn_autocov$nn.idx]

    ## Covariance
    dt_cov <- adaptiveFTS::estimate_autocov(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = grid_cov[, s], t = grid_cov[, t], lag = 0,
        optbw_s = grid_cov[, optbw_s], optbw_t = grid_cov[, optbw_t],
        bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
        correct_diagonal = TRUE, kernel_name = kernel_name)

    dt_cov_dcast <- data.table::dcast(data = dt_cov[order(s,t)], formula = s ~ t, value.var = "autocov")
    c0hat <- as.matrix(dt_cov_dcast[, .SD, .SDcols = ! "s"])
    colnames(c0hat) <- NULL
    c0hat <- (c0hat + t(c0hat)) / 2

    ## lag-1 Autocovariance
    dt_autocov <- adaptiveFTS::estimate_autocov(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = grid_autocov[, s], t = grid_autocov[, t], lag = 1,
        optbw_s = grid_autocov[, optbw_s], optbw_t = grid_autocov[, optbw_t],
        bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
        correct_diagonal = FALSE, kernel_name = kernel_name)

    dt_autocov_dcast <- data.table::dcast(data = dt_autocov[order(s,t)], formula = s ~ t, value.var = "autocov")
    c1hat <- as.matrix(dt_autocov_dcast[, .SD, .SDcols = ! "s"])
    colnames(c1hat) <- NULL

    ## Estimate sigma^2, the error variance, at Tn0 points. If homoscedastic
    ## (default) use a single constant, the median over Tn0; otherwise keep the
    ## t-varying vector. The noise term is diag(sigma2 * rho) either way.
    sigma2 <- adaptiveFTS::estimate_sigma(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X", t = Tn0)[, sig ** 2]
    if (homoscedastic) sigma2 <- stats::median(sigma2, na.rm = TRUE)

    ## Compute the BLUP
    ### Build the variance matrix
    V <- root_Dn0 %*% c0hat %*% root_Dn0 + diag(sigma2 * rho) + tikhonov_reg_param * diag(Mn0)

    ### Compute the blup
    Yn0_centred_weighted <- root_Dn0 %*% matrix(data = Yn0 - muhat_Tn0, ncol = 1)
    blup <- muhat_prediction_points + t(c1hat) %*% root_Dn0 %*% solve(V, Yn0_centred_weighted)

    dt_blup <- data.table::data.table("t" = prediction_points, "blup" = as.vector(blup))

    res <- list(
        "blup" = dt_blup,
        "muhat_Tn0" = muhat_Tn0,
        "muhat_prediction_points" = muhat_prediction_points,
        "c0hat" = c0hat,
        "c1hat" = c1hat,
        "sigma2" = sigma2,
        "homoscedastic" = homoscedastic,
        "rho" = rho,
        "V" = V,
        "tikhonov_reg_param" = tikhonov_reg_param,
        "Tn0" = Tn0,
        "Yn0" = Yn0,
        "bw_grid" = bw_grid,
        "density_bw" = density_bw,
        "dt_risk_mean" = dt_risk_mean,
        "dt_optbw_mean" = dt_optbw_mean,
        "dt_risk_cov" = dt_risk_cov,
        "dt_risk_autocov" = dt_risk_autocov,
        "dt_optbw_cov" = dt_optbw_cov,
        "dt_optbw_autocov" = dt_optbw_autocov,
        "grid_cov" = grid_cov,
        "grid_autocov" = grid_autocov
    )

    return(res)
}


## Example / sanity check

res_blup_one <- predict_next_curve(
    data = data_train, prediction_points = t0,
    tikhonov_reg_param = 1e-6,
    bw_grid = bw_grid_blup, kernel_name = "epanechnikov"
)

dt_blup <- merge(
    data_test[, .(tobs, X)], res_blup_one$blup, by.x = "tobs", by.y = "t", all.x = TRUE
)

dt_graph <- rbind(
  dt_blup[, .("t" = tobs, "Quantity" = "blup", value = blup)],
  dt_blup[, .("t" = tobs, "Quantity" = "Xtrue", value = X)]
)

ggplot(data = dt_graph, mapping = aes(x = t, y = value, group = Quantity, colour = Quantity)) +
  geom_line() +
  xlab("t") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.width= unit(0.8, 'cm'),
        legend.position = "top")

#==========================================================================
## Cross-validation to select the optimal Tikhonov regularization parameter
#==========================================================================


#' Re-evaluate the mean at new locations using cached bandwidths
#'
#' Reuses the adaptive optimal bandwidths already selected inside
#' `predict_next_curve` (`dt_optbw_mean`), matched by nearest neighbour to
#' the requested locations, so that only the plug-in `estimate_mean` is re-run;
#' the risk minimisation is not repeated.
#'
#' @param fit Output list of `predict_next_curve`.
#' @param data Training data used for the plug-in estimate.
#' @param t Evaluation locations.
#' @param kernel_name Kernel name.
#'
#' @return A vector of mean estimates at locations `t`.
#' @keywords internal
mean_at <- function(fit, data, t, kernel_name = "epanechnikov") {
    optbw <- fit$dt_optbw_mean
    knn <- RANN::nn2(data = matrix(optbw$t, ncol = 1), query = matrix(t, ncol = 1), k = 1)
    optbw_t <- optbw$optbw[knn$nn.idx]

    dt <- adaptiveFTS::estimate_mean(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = t, optbw = optbw_t, bw_grid = NULL,
        kernel_name = kernel_name)

    return(dt[order(t), muhat])
}

#' Re-evaluate a (auto)covariance block at new locations using cached bandwidths
#'
#' Reuses the adaptive optimal bandwidths already selected inside
#' `predict_next_curve` (`dt_optbw_cov` for `lag = 0`, `dt_optbw_autocov` for
#' `lag = 1`), matched by nearest neighbour to the requested locations, so that
#' only the plug-in `estimate_autocov` is re-run; the risk minimisation is not
#' repeated.
#'
#' @param fit Output list of `predict_next_curve`.
#' @param data Training data used for the plug-in estimate.
#' @param s,t Evaluation locations (rows indexed by `s`, columns by `t`).
#' @param lag 0 for the covariance, 1 for the lag-1 autocovariance.
#' @param kernel_name Kernel name.
#'
#' @return A `length(s)` x `length(t)` matrix of \eqn{\hat c_{lag}(s_i, t_j)}.
#' @keywords internal
autocov_at <- function(fit, data, s, t, lag = 1L, kernel_name = "epanechnikov") {
    optbw <- if (lag == 0L) fit$dt_optbw_cov else fit$dt_optbw_autocov
    grid <- data.table::as.data.table(expand.grid("s" = s, "t" = t))
    knn <- RANN::nn2(data = optbw[, 1:2], query = grid[, 1:2], k = 1)
    grid$optbw_s <- optbw$optbw_s[knn$nn.idx]
    grid$optbw_t <- optbw$optbw_t[knn$nn.idx]

    dt <- adaptiveFTS::estimate_autocov(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = grid[, s], t = grid[, t], lag = lag,
        optbw_s = grid[, optbw_s], optbw_t = grid[, optbw_t],
        bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
        correct_diagonal = (lag == 0L), kernel_name = kernel_name)

    dt_dcast <- data.table::dcast(dt[order(s, t)], formula = s ~ t, value.var = "autocov")
    m <- as.matrix(dt_dcast[, .SD, .SDcols = ! "s"])
    colnames(m) <- NULL
    return(m)
}


#' One-step-ahead cross-validation for the Tikhonov regularization parameter
#'
#' Selects \eqn{\alpha} by one-step-ahead cross-validation: each of the last
#' `n_val` curves is predicted from its immediate predecessor and scored by the
#' \emph{design-weighted} squared prediction error at its observation points,
#' \eqn{\sum_i \varrho_{n,i} (Y_{n,i} - \widehat X_n(T_{n,i};\alpha))^2}, where
#' \eqn{\varrho_{n,i}} is the design weight of the held-out (target) curve
#' (\eqn{1/|T_n|} under the common design; the importance-sampling correction
#' \eqn{\{M_n \widehat g(T_{n,i})\}^{-1}} under the independent design). This
#' makes the score target the \eqn{\LL^2(I)} risk rather than the \eqn{\LL^2(g)}
#' risk; the \eqn{\alpha}-independent innovation variance shifts every score
#' equally and does not move the argmin. Two regimes:
#' \itemize{
#'   \item \strong{Common design} — the design points never change, so the mean
#'     and (auto)covariance operators are estimated once on the training block
#'     (holdout scheme of Zhao (2026)) and only the conditioning values
#'     \eqn{Y_{v-1}} vary across the validation set.
#'   \item \strong{Independent design} — a \emph{rolling origin}: for each target
#'     curve the operators are re-estimated on all curves observed up to its
#'     predecessor (the training block grown with every earlier validation
#'     curve), so the estimator sees exactly the information available at
#'     prediction time. The adaptive bandwidths are held fixed at those selected
#'     once on the initial block (cached in `fit`); only the plug-in estimates
#'     are refreshed on the growing window (`mean_at`, `autocov_at`).
#' }
#' In both regimes only the \eqn{M \times M} system depends on \eqn{\alpha},
#' which Phase 2 solves for the whole grid from a single eigendecomposition.
#'
#' @param data A prepared functional data.table (`id_curve`, `tobs`, `X`).
#' @param alpha_grid Candidate values. If `NULL`, a 25-point grid
#'   \eqn{\{e^{-5}, e^{-5 + 5/24}, \ldots, e^0\}} is used.
#' @param n_val Number of trailing curves used for one-step-ahead validation.
#' @param bw_grid,kernel_name Passed through to `predict_next_curve`.
#' @param homoscedastic Passed through to `predict_next_curve`: if `TRUE`
#'   (default) a constant \eqn{\sigma^2} (median of the pointwise estimates) is
#'   used, otherwise the t-varying estimates are kept.
#'
#' @return A list with `alpha_star`, `alpha_grid`, `cv_curve`, `cv_matrix`
#'   (fold by alpha), and `val_ids`.
#' @export
cv_alpha_blup <- function(data = data_train, alpha_grid = NULL, n_val = 30L,
                          bw_grid = bw_grid_blup, kernel_name = "epanechnikov",
                          homoscedastic = TRUE, density_bw = NULL) {

    ids <- data[, sort(unique(id_curve))]
    n <- length(ids)
    fit_ids <- ids[seq_len(n - n_val)]
    val_pos <- (n - n_val + 1L):n
    data_fit <- data[id_curve %in% fit_ids]
    is_common <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")

    ## Select the design-density bandwidth once on the initial training block and
    ## reuse it for every estimate_density() call below (independent design only).
    if (!is_common && is.null(density_bw))
        density_bw <- get_density_optimal_bw(
            data = data_fit, idcol = "id_curve", tcol = "tobs", ycol = "X",
            kernel_name = kernel_name, lower = 0, upper = 1)

    ## Estimate every alpha-free component once on the training block
    fit <- predict_next_curve(
        data = data_fit, prediction_points = data_fit[id_curve == max(fit_ids), sort(unique(tobs))],
        tikhonov_reg_param = 1e-6, bw_grid = bw_grid, kernel_name = kernel_name,
        homoscedastic = homoscedastic, density_bw = density_bw)

    ## Phase 1 : assemble the alpha-free pieces for each validation curve
    folds <- vector("list", n_val)
    for (k in seq_len(n_val)) {
        id_targ <- ids[val_pos[k]]
        id_prev <- ids[val_pos[k] - 1L]
        Y_prev <- data[id_curve == id_prev][order(tobs), X]
        Y_targ <- data[id_curve == id_targ][order(tobs), X]

        if (is_common) {
            root_Dn0 <- diag(sqrt(fit$rho))
            ## Design weights of the held-out (target) curve for the validation loss
            rho_targ <- rep(1 / length(Y_targ), length(Y_targ))
            folds[[k]] <- list(
                A0 = root_Dn0 %*% fit$c0hat %*% root_Dn0 + diag(fit$sigma2 * fit$rho),
                C1rD = t(fit$c1hat) %*% root_Dn0,
                mu_pred = fit$muhat_prediction_points,
                resid = root_Dn0 %*% matrix(Y_prev - fit$muhat_Tn0, ncol = 1),
                Y_targ = Y_targ,
                rho_targ = rho_targ,
                Mn0 = length(fit$rho))
        } else {
            ## Rolling origin: re-estimate the operators on every curve observed
            ## up to the prediction origin (id_prev), i.e. the initial training
            ## block grown with all earlier validation curves. Equivalent to the
            ## rbind-growth scheme, expressed as a prefix of the sorted ids. The
            ## bandwidths stay fixed at those selected once on the initial block
            ## (cached in `fit`), so only the plug-in estimates are refreshed.
            data_roll <- data[id_curve <= id_prev]
            Tprev <- data[id_curve == id_prev, sort(unique(tobs))]
            Ttarg <- data[id_curve == id_targ, sort(unique(tobs))]
            c0hat <- autocov_at(fit, data_roll, Tprev, Tprev, lag = 0, kernel_name)
            c0hat <- (c0hat + t(c0hat)) / 2
            c1hat <- autocov_at(fit, data_roll, Tprev, Ttarg, lag = 1, kernel_name)
            mu_prev <- mean_at(fit, data_roll, Tprev, kernel_name)
            mu_targ <- mean_at(fit, data_roll, Ttarg, kernel_name)
            sigma2 <- adaptiveFTS::estimate_sigma(
                data = data_roll, idcol = "id_curve", tcol = "tobs", ycol = "X", t = Tprev)[, sig ** 2]
            if (homoscedastic) sigma2 <- stats::median(sigma2, na.rm = TRUE)
            ghat <- estimate_density(x = Tprev, h = density_bw, kernel_name = kernel_name,
                                     lower = 0, upper = 1)$estimate
            rho <- 1 / (length(Tprev) * pmax(ghat, 1e-6))
            rho <- rho / sum(rho)
            root_Dn0 <- diag(sqrt(rho))

            ## Design weights of the held-out (target) curve for the validation loss
            ghat_targ <- estimate_density(x = Ttarg, h = density_bw, kernel_name = kernel_name,
                                          lower = 0, upper = 1)$estimate
            rho_targ <- 1 / (length(Ttarg) * pmax(ghat_targ, 1e-6))
            rho_targ <- rho_targ / sum(rho_targ)

            folds[[k]] <- list(
                A0 = root_Dn0 %*% c0hat %*% root_Dn0 + diag(sigma2 * rho),
                C1rD = t(c1hat) %*% root_Dn0,
                mu_pred = mu_targ,
                resid = root_Dn0 %*% matrix(Y_prev - mu_prev, ncol = 1),
                Y_targ = Y_targ,
                rho_targ = rho_targ,
                Mn0 = length(rho))
        }
    }

    if (is.null(alpha_grid)) {
        alpha_grid <- exp(seq(-5, 0, length.out = 25))
    }

    ## Phase 2 : only the (A0 + alpha I)^{-1} step depends on alpha. A0 is
    ## symmetric (root_Dn0 %*% c0hat %*% root_Dn0 + a diagonal), so eigendecompose
    ## it once per fold, A0 = Q Lambda Q^T, and reuse it for every alpha:
    ##   (A0 + alpha I)^{-1} resid = Q diag(1 / (lambda + alpha)) Q^T resid.
    ## Each alpha then costs an O(M) rescale instead of an O(M^3) solve, and the
    ## conditioning is explicit in (lambda + alpha).
    cv_matrix <- matrix(NA_real_, nrow = n_val, ncol = length(alpha_grid))
    for (k in seq_len(n_val)) {
        f <- folds[[k]]
        eg <- eigen(f$A0, symmetric = TRUE)
        lambda <- eg$values
        z <- as.vector(crossprod(eg$vectors, f$resid))   # Q^T resid
        W <- f$C1rD %*% eg$vectors                        # C1rD Q, reused across alpha
        for (l in seq_along(alpha_grid)) {
            pred <- tryCatch(
                f$mu_pred + as.vector(W %*% (z / (lambda + alpha_grid[l]))),
                error = function(e) rep(NA_real_, length(f$Y_targ)))
            ## Design-weighted validation loss: sum_i varrho_{n,i} (Y_{n,i} - Xhat_i)^2.
            ## Common design collapses to the plain MSE; independent design applies
            ## the importance-sampling correction so the score targets the L^2(I) risk.
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


## =============================================================================
## Example
## =============================================================================
if (FALSE) {

    cv <- cv_alpha_blup(
        data = data_train,
        alpha_grid = exp(seq(0.3, 1, length.out = 25)),
        n_val = 30L,
        bw_grid = bw_grid_blup
    )
    cv$alpha_star

    dt_cv <- data.table::data.table(alpha = cv$alpha_grid, cv = cv$cv_curve)
    ggplot(dt_cv, aes(x = alpha, y = cv)) +
        geom_line() + geom_point() +
        geom_vline(xintercept = cv$alpha_star, linetype = 2, colour = "red") +
        scale_x_log10() +
        labs(x = expression(alpha), y = expression(CV(alpha)),
             title = "Holdout CV for the Tikhonov parameter") +
        theme_minimal()

    fit <- predict_next_curve(
        data = data_train, prediction_points = t0,
        tikhonov_reg_param = cv$alpha_star, bw_grid = bw_grid_blup)

    dt_blup_cv <- merge(
        data_test[, .(tobs, X)], fit$blup, by.x = "tobs", by.y = "t", all.x = TRUE
    )

    dt_graph_cv <- rbind(
        dt_blup_cv[, .("t" = tobs, "Quantity" = "blup", value = blup)],
        dt_blup_cv[, .("t" = tobs, "Quantity" = "Xtrue", value = X)]
    )

    ggplot(data = dt_graph_cv,
           mapping = aes(x = t, y = value, group = Quantity, colour = Quantity)) +
        geom_line() +
        xlab("t") +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
              axis.title = element_text(size = 12),
              axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.key.width = unit(0.8, 'cm'),
              legend.position = "top")
}
