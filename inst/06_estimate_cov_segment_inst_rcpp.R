library(data.table)
Rcpp::sourceCpp("./src/06_estimate_cov_segment_cpp.cpp")

# Import the data
data("data_far")

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 25)


# Estimate risk function of the mean function
res_cov_segment_risk_cpp <- estimate_cov_segment_risk_cpp(
  data = data_far,
  t = t0, bw_grid = NULL,
  center = TRUE,
  kernel_name = "epanechnikov")

dt_cov <- as.data.table(res_cov_segment_risk_cpp)
dygraphs::dygraph(dt_cov[V1 == t0[1], .(V2, V10)][V10 <50])

# Estimate autocv segment function
res_cov_segment_estim <- estimate_cov_segment_cpp(
  data = data_far, t = t0, optbw = NULL,
  bw_grid = NULL, center = TRUE,
  kernel_name = "epanechnikov")

dt_covseg_estim <- data.table::as.data.table(res_cov_segment_estim)
names(dt_covseg_estim) <- c("t", "optbw", "Ht", "Lt", "PN", "covseghat", "corr_term", "covseghat_corrected")

dygraphs::dygraph(dt_covseg_estim[, .(t, covseghat, covseghat_corrected)])
