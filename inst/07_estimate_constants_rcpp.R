source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
source("./R/07_estimate_constants.R")
Rcpp::sourceCpp("./src/07_estimate_constants_cpp.cpp")

# Load the data
data("data_far")

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 10)
presmooth_bw <- rep(0.026, data_far[, length(unique(id_curve))])

# Estimate the empirical autocov
dt_autocov_rcpp <- estimate_empirical_autocov_cpp(
  data = data_far, t = t0, lag = c(0, 1, 2),
  h = presmooth_bw, kernel_name = "epanechnikov")
dt_autocov_rcpp

# Estimate XsXt autocovariance
t0 <- seq(0.1, 0.7, len = 10)
s0 <- seq(0.3, 0.9, len = 10)
dt_XsXt_autocov_rcpp <- estimate_empirical_XsXt_autocov_cpp(
  data = data_far, t = t0, s = s0, lag = 0:148,
  cross_lag = 1, h = presmooth_bw, kernel_name = "epanechnikov", center = FALSE)
dt_XsXt_autocov_rcpp

# Estimate sigma ----
dt_sig_cpp <- estimate_sigma_cpp(data = data_far, t = t0)
dt_sig_cpp

# Estimate moment ----
mat_mom <- estimate_empirical_mom_cpp(
  data = data_far, t = t0,
  h = presmooth_bw,  mom_order = 1,
  center = FALSE, kernel_name = "epanechnikov")
mat_mom


## Estimate \mathbb{D}(t,h_t) numerator ----

estimate_numerator_dependence_term_DD_single_t_cpp(
  data = data_far, t = t0[1], max_lag = 6, bw = 0.1,
  h = presmooth_bw, kernel_name = "epanechnikov", center = TRUE)

mat_num_DD <- estimate_numerator_dependence_term_DD_cpp(
  data = data_far, t = unique(t0), max_lag = 3, bw_vec = bw_grid,
  h = presmooth_bw, kernel_name = "epanechnikov", center = TRUE)
mat_num_DD
mat_num_DD[, 3] / mat_num_DD[, 4] ** 3
