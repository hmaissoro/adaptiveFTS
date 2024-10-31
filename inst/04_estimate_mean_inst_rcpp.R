Rcpp::sourceCpp("./src/04_estimate_mean_cpp.cpp")

# Import the data
data("data_far")

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 25)


# Estimate risk function of the mean function
res_mean_risk_cpp <- estimate_mean_risk_cpp(
  data = data_far,
  t = t0, bw_grid = NULL,
  kernel_name = "epanechnikov")

res_mean_risk_cpp

# Estimate mean function
res_mean_estim <- estimate_mean_cpp(
  data = data_far, t = t0, optbw = NULL,
  bw_grid = NULL, kernel_name = "epanechnikov")
dt_mean_estim <- data.table::as.data.table(res_mean_estim)
names(dt_mean_estim) <- c("t", "optbw", "Ht", "Lt", "PN", "muhat")

dygraphs::dygraph(dt_mean_estim[, .(t, muhat)])
