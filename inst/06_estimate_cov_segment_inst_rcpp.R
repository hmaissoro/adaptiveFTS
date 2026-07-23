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
plot(dt_cov[V1 == t0[1], .(V2, V10)][V10 < 50], type = "l",
     xlab = "h", ylab = "risk function")

# Estimate autocv segment function
res_cov_segment_estim <- estimate_cov_segment_cpp(
  data = data_far, t = t0, optbw = NULL,
  bw_grid = NULL, center = TRUE,
  kernel_name = "epanechnikov")

dt_covseg_estim <- data.table::as.data.table(res_cov_segment_estim)
names(dt_covseg_estim) <- c("t", "optbw", "Ht", "Lt2", "PN", "covseghat", "corr_term", "covseghat_corrected")

matplot(dt_covseg_estim[, t], dt_covseg_estim[, .(covseghat, covseghat_corrected)],
        type = "l", lty = 1, xlab = "t", ylab = "covariance segment")
