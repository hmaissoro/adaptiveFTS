library(data.table)
library(ggplot2)
Rcpp::sourceCpp("./src/08_estimate_curve_cpp.cpp")
source("./inst/09_parzen_rosenblatt_density_estimator.R")


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
    bw_grid = bw_grid_blup, kernel_name = "epanechnikov") {

    n0 <- data[, max(id_curve)]
    Tn0 <- data[id_curve == n0, sort(unique(tobs))]
    Yn0 <- data[id_curve == n0][order(tobs), X]
    Mn0 <- length(Tn0)

    # Compute the weights \rho_{n0,i} for i = 1, ..., M_{n0}
    is_common_design <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")

    if (is_common_design) {
        rho <- rep(1 / Mn0, Mn0)
    } else {
        ghat <- estimate_density(x = Tn0, kernel_name = kernel_name, lower = 0, upper = 1)$estimate
        rho <- 1 / (Mn0 * ghat)
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

    # Estimate the mean function at the prediction points
    dt_muhat_prediction_points <- adaptiveFTS::estimate_mean(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = prediction_points, optbw = NULL, bw_grid = bw_grid,
        kernel_name = kernel_name)
    muhat_prediction_points <- dt_muhat_prediction_points[order(t), muhat]

    # Estimate the mean function at Tn0 points
    dt_muhat_Tn0 <- adaptiveFTS::estimate_mean(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = Tn0, optbw = NULL, bw_grid = bw_grid,
        kernel_name = kernel_name)
    muhat_Tn0 <- dt_muhat_Tn0[order(t), muhat]

    # Build estimation locations for the covariance and autocovariance functions
    grid_cov <- expand.grid("s" = Tn0, "t" = Tn0)
    grid_cov <- data.table::as.data.table(grid_cov)
    grid_autocov <- expand.grid("s" = Tn0, "t" = prediction_points)
    grid_autocov <- data.table::as.data.table(grid_autocov)

    # Estimate the (auto)covariance function
    ## To reduce computation time, we estimate the (auto)covariance function 
    ## on a sub-grid and then use nearest neighbor matching to find the optimal 
    ## bandwidths for the full grid.
    sub_grid_vec <- seq(0.05, 0.95, length.out = 10)
    sub_grid <- expand.grid("s" = sub_grid_vec, "t" = sub_grid_vec)

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

    ## Estimate sigma, the error standard deviation, at Tn0 points
    dt_sigma <- adaptiveFTS::estimate_sigma(
        data = data, idcol = "id_curve", 
        tcol = "tobs", ycol = "X", t = Tn0)
    sigma2 <- dt_sigma[, sig ** 2]

    ## Compute the BLUP
    ### Build the variance matrix
    V <- root_Dn0 %*% c0hat %*% root_Dn0 + diag(sigma2 * rho) + tikhonov_reg_param * diag(Mn0)
    invV <- solve(V) 

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
        "rho" = rho,
        "V" = V,
        "invV" = invV,
        "tikhonov_reg_param" = tikhonov_reg_param,
        "Tn0" = Tn0,
        "Yn0" = Yn0,
        "bw_grid" = bw_grid,
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

