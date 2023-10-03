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

# Add noise
dt_far[, X := X + rnorm(n = .N, mean = 0, sd = 0.9 ** (0.1)), by = id_curve]

# Estimate risk function
dt_autocov_risk <- estimate_autocov_risk(
  data = dt_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 3,
  bw_grid = seq(0.005, 0.15, len = 45),
  smooth_ker = epanechnikov,
  Hs = NULL, Ls = NULL,
  Ht = NULL, Lt = NULL,
  Delta = NULL, h = NULL
)

# Plot mean risk function
dt_dcast <- data.table::dcast(data = dt_autocov_risk,
                              formula = h ~ s + t ,
                              value.var = "autocov_risk")

manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.2, 0.25)" = `0.2_0.25`)],
                      main = "lag = 3 - (s, t) = (0.2, 0.25)",
                      xlab = "h",
                      ylab = "risk function"),
    dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.4, 0.5)" = `0.4_0.5`)],
                      main = "lag = 3 - (s, t) = (0.4, 0.5)",
                      xlab = "h",
                      ylab = "risk function"),
    dygraphs::dygraph(data = dt_dcast[, .(h, "(s, t) = (0.8, 0.75)" = `0.8_0.75`)],
                      main = "lag = 3 - (s, t) = (0.8, 0.75)",
                      xlab = "h",
                      ylab = "risk function")
  ),
  nrow = 3
)
