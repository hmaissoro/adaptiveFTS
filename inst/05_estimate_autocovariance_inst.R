# Adaptive lag-l autocovariance estimation, with one and with two bandwidths.

library(data.table)

data("data_far")

# Estimate risk function
## Using one bandwidth
dt_autocov_risk <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
  lag = 3, bw_grid = NULL, common_bw = TRUE,
  center_curves = TRUE, kernel_name = "epanechnikov")

dt_autocov_risk[, .(s, t, hs, ht, autocov_risk)]

summary(dt_autocov_risk)
plot(dt_autocov_risk)

## Using two bandwidths
dt_autocov_risk_2bw <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
  lag = 3, bw_grid = NULL, common_bw = FALSE,
  center_curves = TRUE, kernel_name = "epanechnikov")

dt_autocov_risk_2bw[, .(s, t, hs, ht, autocov_risk)]

# Lag-0 covariance surface on a grid
tgrid <- seq(0.2, 0.8, len = 10)
cov_grid <- data.table::CJ(s = tgrid, t = tgrid)
dt_cov_surface <- estimate_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = cov_grid[, s], t = cov_grid[, t], lag = 0,
  bw_grid = NULL, common_bw = TRUE, center_curves = TRUE,
  kernel_name = "epanechnikov")

summary(dt_cov_surface)
plot(dt_cov_surface)
