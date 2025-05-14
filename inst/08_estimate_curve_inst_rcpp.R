library(ggplot2)
library(data.table)
Rcpp::sourceCpp("./src/08_estimate_curve_cpp.cpp")

# Import the data
data("data_far")

# Prepare the data
data_prepared <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")
data_to_use <- data_prepared[-which(id_curve == 150 & tobs > 0.5)]

# Estimate the BLUP
t0 <- data_prepared[id_curve == 150, tobs]
res_blup_one <- estimate_curve_cpp(
  data = data_to_use, t = t0, id_curve = 150,
  kernel_name = "epanechnikov",
  bw_grid = NULL, use_same_bw = TRUE,
  center = TRUE, correct_diagonal = TRUE)

dt_res <- data.table::as.data.table(res_blup_one$res_blup)
names(dt_res) <- c("t", "muhat", "blup")
dt_res <- data.table::merge.data.table(
  x = dt_res, y = data_prepared[id_curve == 150, .("t" = tobs, "Xtrue" = X)],
  by = "t")

## compute MSE
dt_res[, .(mean((Xtrue - blup) ** 2))]
dt_graph <- rbind(
  dt_res[, .(t, "Quantity" = "blup", value = blup)],
  dt_res[, .(t, "Quantity" = "Xtrue", value = Xtrue)]
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

# Confidence band ----
data <- copy(data_to_use)
kernel_name <- "epanechnikov"
id_curve_to_predict <- 150
Tn0_lag_vec <- data[id_curve == (id_curve_to_predict - 1), sort(tobs)]
Tn0_vec <- data[id_curve == id_curve_to_predict, sort(tobs)]

## Compute the \mathbb S matrix

## Estimate covariance matrice
tcb <- seq(0.001, 0.999, len = 100)
grid_cov_fun <- expand.grid("s" = tcb, "t" = tcb)

### Optimal bw param for covariance function estimation
mat_opt_cov_bw <- get_nearest_best_autocov_bw(
  mat_opt_param = res_blup_one$opt_cov_param,
  snew = grid_cov_fun$s, tnew = grid_cov_fun$t)

### Estimate covariance function
mat_cov_melt <- estimate_autocov_cpp(
  data = data, s = mat_opt_cov_bw[, 1], t = mat_opt_cov_bw[, 2], lag = 0,
  optbw_s = mat_opt_cov_bw[, 3], optbw_t = mat_opt_cov_bw[, 4],
  bw_grid = NULL, use_same_bw = FALSE, center = TRUE, correct_diagonal = FALSE,
  kernel_name = kernel_name)
dt_cov_cb <- data.table::as.data.table(mat_cov_melt[, c(1, 2, 14)])
names(dt_cov_cb) <- c("s", "t", "cov")
mat_cov <- data.table::dcast(data = dt_cov_cb[order(s,t)], formula = s ~ t, value.var = "cov")
mat_cov <- as.matrix(mat_cov[, -1])
colnames(mat_cov) <- NULL

## Estimate Cov_MT(Y_{n_0, 1}, X_{n_0}(t))
### Optimal bw param for Cov_MT(Y_{n_0, 1}, X_{n_0}(t))
grid_autocov <- expand.grid("s" = Tn0_lag_vec, "t" = tcb)
grid_cov <- expand.grid("s" = Tn0_vec, "t" = tcb)

mat_opt_covMT_autocov <- get_nearest_best_autocov_bw(
  mat_opt_param = res_blup_one$opt_autocov_param,
  snew = grid_autocov$s, tnew = grid_autocov$t)

mat_opt_covMT_cov <- get_nearest_best_autocov_bw(
  mat_opt_param = res_blup_one$opt_cov_param,
  snew = grid_cov$s, tnew = grid_cov$t)

### Estimate Cov_MT(Y_{n_0, 1}, X_{n_0}(t))
mat_covMT_autocov_melt <- estimate_autocov_cpp(
  data = data, s = mat_opt_covMT_autocov[, 1], t = mat_opt_covMT_autocov[, 2], lag = 1,
  optbw_s = mat_opt_covMT_autocov[, 3], optbw_t = mat_opt_covMT_autocov[, 4],
  bw_grid = NULL, use_same_bw = FALSE, center = TRUE, correct_diagonal = FALSE,
  kernel_name = kernel_name)

mat_covMT_cov_melt <- estimate_autocov_cpp(
  data = data, s = mat_opt_covMT_cov[, 1], t = mat_opt_covMT_cov[, 2], lag = 0,
  optbw_s = mat_opt_covMT_cov[, 3], optbw_t = mat_opt_covMT_cov[, 4],
  bw_grid = NULL, use_same_bw = FALSE, center = TRUE, correct_diagonal = TRUE,
  kernel_name = kernel_name)

### Reshape covariance and autocovariance
dt_covMT_autocov <- data.table::as.data.table(mat_covMT_autocov_melt[, c(1, 2, 14)])
names(dt_covMT_autocov) <- c("s", "t", "autocov")
mat_covMT_autocov <- data.table::dcast(data = dt_covMT_autocov[order(s,t)], formula = s ~ t, value.var = "autocov")
mat_covMT_autocov <- as.matrix(mat_covMT_autocov[, -1])
colnames(mat_covMT_autocov) <- NULL

dt_covMT_cov <- data.table::as.data.table(mat_covMT_cov_melt[, c(1, 2, 14)])
names(dt_covMT_cov) <- c("s", "t", "cov")
mat_covMT_cov <- data.table::dcast(data = dt_covMT_cov[order(s,t)], formula = s ~ t, value.var = "cov")
mat_covMT_cov <- as.matrix(mat_covMT_cov[, -1])
colnames(mat_covMT_cov) <- NULL

### Build the matrice Cov_MT(Y_{n_0, 1}, X_{n_0}(t))
mat_covMT <- rbind(mat_covMT_autocov, mat_covMT_cov)

## Compute matrix \mathbb S
mat_varY <- res_blup_one$mat_VarY
mat_VarY_inverse <- solve(mat_varY)

## /!\ For the auto covariance matrix we need to add sigma ?
mat_SS <- mat_cov - t(mat_covMT) %*% mat_VarY_inverse %*% mat_covMT


## Estimate pointwise confidence band
SS_vec_sqrt <- diag(mat_SS) ** (1/2)
dt_raw_cb <- data.table::data.table("t" = tcb, "band" = SS_vec_sqrt * qnorm(p = 0.975))
dt_raw_cb <- dt_raw_cb[!is.nan(band)]
dt_band_smooth <- estimate_nw(y = dt_raw_cb[, t], t = dt_raw_cb[, band], tnew = t0, h = 0.1074619 * 1.7, kernel_name = kernel_name)

dt_res_with_cb <- merge(x = dt_res, y = dt_band_smooth[, .("t" = tnew, "band_smooth" = yhat)], by = "t")

dt_res_with_cb


dt_graph_cb <- rbind(
  dt_res_with_cb[, .(t, "Quantity" = "blup", value = blup)],
  dt_res_with_cb[, .(t, "Quantity" = "Xtrue", value = Xtrue)],
  dt_res_with_cb[, .(t, "Quantity" = "lower_band", value = blup - band_smooth)],
  dt_res_with_cb[, .(t, "Quantity" = "upper_band", value = blup + band_smooth)]
)

ggplot(data = dt_graph_cb, mapping = aes(x = t, y = value, group = Quantity, colour = Quantity)) +
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

## Estimate simultaneous confidence band
## Compute matrix \rho(s,t)
SS_vec <- diag(mat_SS) ** (1/2)
SS_vec <- matrix(data = SS_vec, ncol = 1)
mat_SS_ss_tt <- SS_vec %*% t(SS_vec)
dim(mat_SS_ss_tt)

mat_rho <- mat_SS / mat_SS_ss_tt
mat_rho[is.nan(mat_rho)] <- 0
isSymmetric(mat_rho)
mat_rho <- 0.5 * (mat_rho + t(mat_rho))
isSymmetric(mat_rho)

Zt <- MASS::mvrnorm(n = 1, mu = rep(0, ncol(mat_rho)), Sigma = mat_rho)

solve(mat_rho)


dim(res_blup_one$cov_pred_pred)
dim(res_blup_one$covY_Xn0)
dim(res_blup_one$covY_Xn0)

## plot confidence band



# N = 400 lambda = 300 ----
dt_res_blup_N400lambda300 <- data.table::rbindlist(parallel::mclapply(1:20, function(id){
  dt <- readRDS(paste0("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N400lambda300/dt_mc_FAR_mfBm_N=400_lambda=300_id_mc=", id, "_fts_model_2.RDS"))
  dt_common <- dt[ttag == "tcommon"]
  dt <- dt[ttag == "trandom"]
  # dt[! id_curve == 400 & tobs > 0.5]
  data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
  data_to_use <- data_prepared[-which(id_curve == 400 & tobs > 0.5)]

  res_blup_one <- estimate_curve_cpp(data = data_to_use, id_curve = 400, t = t0, kernel_name = "epanechnikov")

  dt_res <- data.table::as.data.table(res_blup_one$res_blup)
  names(dt_res) <- c("t", "muhat", "blup")
  dt_res <- data.table::merge.data.table(x = dt_res, y = dt_common[id_curve == 400, .("t" = tobs, "Xtrue" = X)], by = "t")
  dt_res[, id_mc := id]
  return(dt_res)
}, mc.cores = 20))

saveRDS(dt_res_blup_N400lambda300, "../../curve_reconstruction/adaptive_blup_emp_study/blup/dt_res_blup_N=400_lambda=300_fts_model_2.RDS")

dt_res_blup[, error := Xtrue - blup]
dt_res_blup[, t := as.factor(t)]
ggplt <- ggplot(data = dt_res_blup[error < 10 & error > -10], mapping = aes(x = t, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red") +
  ggtitle(latex2exp::TeX("N = 150, $\\lambda$ = 40")) +
  ylab(latex2exp::TeX("$X_{n_0}(t) - \\hat{X}_{n_0}(t)$")) +
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
ggsave(plot = ggplt, filename = file.path("../../curve_reconstruction/adaptive_blup_emp_study/blup/blup_N=150_lambda=40_fts_model_2.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")



