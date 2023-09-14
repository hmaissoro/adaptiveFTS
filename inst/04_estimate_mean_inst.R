# Generate a FAR A process
dt_far <- simulate_far(N = 50, lambda = 70,
                       tdesign = "random",
                       tdistribution = runif,
                       tcommon = seq(0.2, 0.8, len = 50),
                       hurst_fun = hurst_logistic,
                       L = 4,
                       far_kernel = get_real_data_far_kenel,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)
data <- dt_far[ttag == "trandom"]

# Estimate risk function
dt_mean_risk <- estimate_mean_risk(
  data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  Delta = NULL, h = NULL, smooth_ker = epanechnikov)

dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
dygraphs::dygraph(dt_dcast)

