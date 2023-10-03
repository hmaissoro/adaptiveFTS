# Generate a FAR A process
dt_far <- simulate_far(N = 50, lambda = 70,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = NULL,
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = get_real_data_far_kenel,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)
data <- dt_far[ttag == "trandom"]

# Add noise
dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]

# Estimate risk function
dt_mean_risk <- estimate_mean_risk(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  Delta = NULL, h = NULL, smooth_ker = epanechnikov)

# Plot mean risk function
dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")

manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(dt_dcast[, .(h, "t = 0.25" = `0.25`)], main = "t = 0.25", xlab = "h", ylab = "risk function"),
    dygraphs::dygraph(dt_dcast[, .(h, "t = 0.5" = `0.5`)], main = "t = 0.5", xlab = "h", ylab = "risk function"),
    dygraphs::dygraph(dt_dcast[, .(h, "t = 0.75" = `0.75`)], main = "t = 0.75", xlab = "h", ylab = "risk function")
  ),
  nrow = 3
)

# Estimate mean function
dt_mean <- estimate_mean(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  Delta = NULL, h = NULL, smooth_ker = epanechnikov)

# Table of the estimates of the mean function
DT::datatable(data = dt_mean[, lapply(.SD, function(X) round(X, 5))])

# Estimate mean function using Rubìn and Panaretos (2020) method
## Estimate the bandwidth by Cross-Validation
dt_bw_mean_rp <- estimate_mean_bw_rp(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  Kfold = 10, bw_grid = seq(0.001, 0.15, len = 45),
  smooth_ker = epanechnikov)

## Plot the Cross-Validation error
dygraphs::dygraph(dt_bw_mean_rp)

## Select the best bandwidth
optbw <- dt_bw_mean_rp[, h[which.min(cv_error)]]

## Estimate the mean function
dt_mean_rp <- estimate_mean_rp(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), h = optbw, smooth_ker = epanechnikov)

DT::datatable(data = dt_mean_rp[, lapply(.SD, function(X) round(X, 5))])





