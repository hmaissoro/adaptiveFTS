source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
source("./R/07_estimate_constants.R")
source("./R/04_estimate_mean.R")
Rcpp::sourceCpp("./src/04_estimate_mean_cpp.cpp")

# Import the data
dt <- readRDS("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=1_fts_model_2.RDS")
dt <- dt[ttag == "trandom"]

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 25)
presmooth_bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")

## Estimate locreg
lambdahat <- mean(dt[, .N, by = "id_curve"][, N])
Deltahat <- 2 * exp(-log(lambdahat) ** 0.72)
dt_locreg <- estimate_locreg(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, Delta = Deltahat, h = presmooth_bw,
  smooth_ker = epanechnikov, center = TRUE)
Htvec <- dt_locreg[order(t), Ht]
Ltvec <- dt_locreg[order(t), Lt]

# Estimate the mean risk

## bw grid
N <- length(dt[, unique(id_curve)])
K <- 20
b0 <- (N * lambdahat) ** (- 0.9)
bK <- 2 * (N * lambdahat) ** (- 1 / 3)
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(K, b0, bK, a) ; gc()

## Current implementation
dt_mean_risk <- estimate_mean_risk(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, bw_grid = bw_grid,
  Ht = NULL, Lt = NULL, Delta = NULL, h = NULL,
  smooth_ker = epanechnikov)

## C++ Current implementation
res_mean_risk_cpp <- estimate_mean_risk_cpp(
  data = data_prepared,
  t = t0, bw_grid = bw_grid,
  kernel_name = "epanechnikov")

estimate_mean_risk_cpp(
  data = data_prepared,
  t = t0, bw_grid = NULL,
  kernel_name = "epanechnikov")

dt_mean_risk_cpp <- data.table::as.data.table(res_mean_risk_cpp)
names(dt_mean_risk_cpp) <- c("t", "h", "PN_cpp", "bias_term_cpp", "variance_term_cpp", "dependence_term_cpp", "mean_risk_cpp")
dt_mean_risk_all <- data.table::merge.data.table(x = dt_mean_risk, y = dt_mean_risk_cpp, by = c("t", "h"))

all.equal(target = dt_mean_risk_all[, bias_term_cpp], current = dt_mean_risk_all[, bias_term])
all.equal(target = dt_mean_risk_all[, variance_term_cpp], current = dt_mean_risk_all[, varriance_term])
all.equal(target = dt_mean_risk_all[, dependence_term_cpp], current = dt_mean_risk_all[, dependence_term])

## Graph
dygraphs::dygraph(dt_mean_risk_all[t == t0[1], .(h, mean_risk, mean_risk_cpp)])


## Benchmark
res_bm_mean_risk <- microbenchmark::microbenchmark(
  mean_risk_cur = dt_mean_risk <- estimate_mean_risk(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t0, bw_grid = bw_grid,
    Ht = Htvec, Lt = Ltvec, Delta = NULL, h = presmooth_bw,
    smooth_ker = epanechnikov),
  mean_risk_cpp = res_mean_risk_cpp <- estimate_mean_risk_cpp(
    data = data_prepared,
    t = t0, bw_grid = bw_grid,
    Ht = Htvec, Lt = Ltvec, h = rep(median(presmooth_bw), length(presmooth_bw)),
    kernel_name = "epanechnikov"),
  times = 10
)
print(res_bm_mean_risk, unit = "relative", order = "median")

# Estimate mean function ----
res_mean_estim <- estimate_mean_cpp(data = data_prepared, t = t0, optbw = NULL, bw_grid = NULL, kernel_name = "epanechnikov")
dt_mean_estim <- data.table::as.data.table(res_mean_estim)
names(dt_mean_estim) <- c("t", "optbw", "Ht", "Lt", "PN", "muhat")

dygraphs::dygraph(dt_mean_estim[, .(t, muhat)])
