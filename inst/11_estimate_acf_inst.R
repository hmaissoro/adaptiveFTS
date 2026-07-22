# Adaptive functional autocorrelation function (FACF) demo.
#
# estimate_acf() computes the functional autocorrelation
#   rho_l = || Gamma_l || / integral Gamma_0(t, t) dt
# from the adaptive (auto)covariance estimator. data_far is a FAR(1) process, so
# the FACF is expected to decay with the lag.

# Load data
data("data_far")

# Estimate the FACF for lags 1..8 on a 20-point evaluation grid.
# The grid is chosen design-aware automatically; here data_far is an
# independent design, so a regular grid is used.
dt_acf <- estimate_acf(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  lag.max = 8, n_grid = 20, use_same_bw = FALSE,
  center = TRUE, kernel_name = "epanechnikov")

# The classed data.table: per-lag L2 norm and rho (see ?adaptiveFTS_est).
dt_acf

# Compact, design-aware text summary (design, grid, largest |rho|, ...).
summary(dt_acf)

# ACF-style lag plot of rho_l (requires ggplot2).
plot(dt_acf)
# equivalently: ggplot2::autoplot(dt_acf)

# A tabular view of the estimates.
DT::datatable(dt_acf[, lapply(.SD, function(x) round(x, 4))])
