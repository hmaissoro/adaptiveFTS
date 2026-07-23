# Adaptive functional autocorrelation function (FACF) demo.

# Load data
data("data_far")

# Estimate the FACF for lags 1..8 on a 20-point grid
dt_facf <- estimate_facf(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  lag.max = 8, n_grid = 20, use_same_bw = FALSE,
  center = TRUE, kernel_name = "epanechnikov")

# Summary and FACF lag plot
dt_facf
summary(dt_facf)
plot(dt_facf)

# Tabular view
dt_facf[, lapply(.SD, function(x) round(x, 4))]
