source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
source("./R/07_estimate_constants.R")
Rcpp::sourceCpp("./src/07_estimate_constants_cpp.cpp")

# Import the data
dt <- readRDS("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=1_fts_model_2.RDS")
dt <- dt[ttag == "trandom"]

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 10)
presmooth_bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")

## Estimate autocovariance ----
# Estimation using current function
dt_autocov <- estimate_empirical_autocov(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, lag = c(0, 1, 2), h = presmooth_bw,
  smooth_ker = epanechnikov)

# Estimation using Rcpp
dt_autocov_rcpp <- estimate_empirical_autocov_cpp(
  data = data_prepared, t = t0, lag = c(0, 1, 2),
  h = presmooth_bw, kernel_name = "epanechnikov")

all.equal(target = dt_autocov[, autocov], current = dt_autocov_rcpp[, 3])

# Microbenchmark
res_benchmark <- microbenchmark::microbenchmark(
  Autocov_cur = dt_autocov <- estimate_empirical_autocov(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t0, lag = c(0, 1, 2), h = presmooth_bw,
    smooth_ker = epanechnikov),
  Autocov_rcpp = dt_autocov_rcpp <- estimate_empirical_autocov_cpp(
    data = data_prepared, t = t0, lag = c(0, 1, 2),
    h = presmooth_bw, kernel_name = "epanechnikov"),
  times = 50
)

print(res_benchmark, unit = "relative", order = "median")

## Estimate XsXt autocovariance ----
t0 <- seq(0.1, 0.7, len = 10)
s0 <- seq(0.3, 0.9, len = 10)
dt_XsXt_autocov <- estimate_empirical_XsXt_autocov(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, s = s0, cross_lag = 1, lag = 0:148, h = presmooth_bw,
  smooth_ker = epanechnikov, center = FALSE)

# Estimation using Rcpp
dt_XsXt_autocov_rcpp <- estimate_empirical_XsXt_autocov_cpp(
  data = data_prepared, t = t0, s = s0, lag = 0:148,
  cross_lag = 1, h = presmooth_bw, kernel_name = "epanechnikov", center = FALSE)

dt_XsXt_autocov_rcpp <- data.table::as.data.table(dt_XsXt_autocov_rcpp)
names(dt_XsXt_autocov_rcpp) <- c("s", "t", "cross_lag", "lag", "EXsXt_cross_lag_cpp", "XsXt_autocov_cpp")
dt_autocov_XsXt <- data.table::merge.data.table(
  x = dt_XsXt_autocov,
  y = dt_XsXt_autocov_rcpp,
  by = c("s", "t", "cross_lag", "lag"))

all.equal(target = dt_autocov_XsXt[, EXsXt_cross_lag], current = dt_autocov_XsXt[, EXsXt_cross_lag_cpp])
all.equal(target = dt_autocov_XsXt[, XsXt_autocov], current = dt_autocov_XsXt[, XsXt_autocov_cpp])

# Microbenchmark
res_benchmark_XsXt <- microbenchmark::microbenchmark(
  XsXt_autocov_cur = dt_XsXt_autocov <- estimate_empirical_XsXt_autocov(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t0, s = s0, cross_lag = 1, lag = c(0, 1, 2), h = presmooth_bw,
    smooth_ker = epanechnikov, center = TRUE),
  XsXt_autocov_rcpp = dt_XsXt_autocov_rcpp <- estimate_empirical_XsXt_autocov_cpp(
    data = data_prepared, t = t0, s = s0, lag = c(0, 1, 2),
    cross_lag = 1, h = presmooth_bw, kernel_name = "epanechnikov", center = TRUE),
  times = 50
)

print(res_benchmark_XsXt, unit = "relative", order = "median")

# Estimate sigma ----
dt_sig <- estimate_sigma(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X", t = t0)
dt_sig_cpp <- estimate_sigma_cpp(data = data_prepared, t = t0)

all.equal(target = dt_sig[, sig], current = dt_sig_cpp[, 2])

res_bm_sig <- microbenchmark::microbenchmark(
  sig_cur = dt_sig <- estimate_sigma(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X", t = t0),
  sig_cpp = dt_sig_cpp <- estimate_sigma_cpp(data = data_prepared, t = t0),
  times = 50
)

print(res_bm_sig, unit = "relative", order = "median")

# Estimate moment ----
mat_mom <- estimate_empirical_mom_cpp(
  data = data_prepared, t = t0,
  h = presmooth_bw,  mom_order = 1,
  center = FALSE, kernel_name = "epanechnikov")
