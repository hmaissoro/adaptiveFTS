Rcpp::sourceCpp("./src/05_estimate_autocov_cpp.cpp")

# Import the data
data("data_far")

# Estimation parameters
t0 <- c(0.2, 0.4, 0.5, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]
rm(dt_st) ; gc()

# One bandwidth
res_gammahat_risk_cpp <- estimate_autocov_risk_cpp(
  data = data_far, s = s0, t = t0, lag = 1, bw_grid = NULL,
  use_same_bw = TRUE, center = FALSE,
  kernel_name = "epanechnikov")
res_gammahat_risk_cpp

# two bandwidth
res_gammahat_risk_2bw_cpp <- estimate_autocov_risk_cpp(
  data = data_far, s = s0, t = t0, lag = 1, bw_grid = NULL,
  use_same_bw = FALSE, center = TRUE, kernel_name = "epanechnikov")
res_gammahat_risk_2bw_cpp


# autoCovariance function

res_cov <- estimate_autocov_cpp(
  data = data_far, s = s0, t = t0, lag = 0,
  use_same_bw = TRUE, center = FALSE, optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
  correct_diagonal = TRUE, kernel_name = "epanechnikov")
res_cov

res_gammahat_risk_2bw_cpp <- estimate_autocov_risk_cpp(
  data = data_far, s = s0, t = t0, lag = 1, bw_grid = NULL,
  use_same_bw = TRUE, center = TRUE, kernel_name = "epanechnikov")
