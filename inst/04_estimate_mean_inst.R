# Adaptive mean function estimation, compared with the Rubin-Panaretos estimator.

library(data.table)
library(ggplot2)

data("data_far")

# Estimate risk function
dt_mean_risk <- estimate_mean_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL, kernel_name = "epanechnikov")

# Risk against the bandwidth, one panel per observation point
summary(dt_mean_risk)
plot(dt_mean_risk)

# Estimate mean function
dt_mean <- estimate_mean(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL, kernel_name = "epanechnikov")

dt_mean[, lapply(.SD, function(X) round(X, 5))]
summary(dt_mean)
plot(dt_mean)

# Estimate mean function using Rubìn and Panaretos (2020) method
## Estimate the bandwidth by Cross-Validation
dt_bw_mean_rp <- estimate_mean_bw_rp(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  n_folds = 10, bw_grid = seq(0.001, 0.15, len = 45),
  kernel_name = "epanechnikov")

## Plot the Cross-Validation error
ggplot(dt_bw_mean_rp, aes(x = h, y = cv_error)) +
  geom_line(colour = "#1B4F72") +
  labs(x = "h", y = "cross-validation error") +
  theme_minimal()

## Select the best bandwidth
best_bw <- dt_bw_mean_rp[, h[which.min(cv_error)]]

## Estimate the mean function
dt_mean_rp <- estimate_mean_rp(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), h = best_bw, kernel_name = "epanechnikov")

dt_mean_rp[, lapply(.SD, function(X) round(X, 5))]
