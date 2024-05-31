source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
source("./R/07_estimate_constants.R")
source("./R/04_estimate_mean.R")
source("./R/05_estimate_autocovariance.R")
source("./R/07_estimate_autocovariance_2bw.R")
Rcpp::sourceCpp("./src/05_estimate_autocov_cpp.cpp")

# Import the data
dt <- readRDS("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=1_fts_model_2.RDS")
dt <- dt[ttag == "trandom"]

# Estimation parameters
t0 <- c(0.2, 0.4, 0.5, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]
rm(dt_st) ; gc()

presmooth_bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
data_prepared <- .format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")

## Estimate locreg
lambdahat <- mean(dt[, .N, by = "id_curve"][, N])
Deltahat <- 2 * exp(-log(lambdahat) ** 0.72)
dt_locreg_s <- estimate_locreg(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = s0, Delta = Deltahat, h = presmooth_bw,
  smooth_ker = epanechnikov, center = TRUE)
dt_locreg_t <- estimate_locreg(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, Delta = Deltahat, h = presmooth_bw,
  smooth_ker = epanechnikov, center = TRUE)
Hvec_s <- dt_locreg_s[, Ht]
Lvec_s <- dt_locreg_s[, Lt]
Hvec_t <- dt_locreg_t[, Ht]
Lvec_t <- dt_locreg_t[, Lt]

# Estimate the autocovariance risk

## bw grid
N <- length(presmooth_bw)
K <- 20
b0 <- 4 * max((N * lambdahat) ** (- 0.9), (N * (lambdahat ** 2)) ** (- 0.9))
bK <- 5 / 10
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(K, b0, bK, a) ; gc()

## One bandwidth
### Current implementation
dt_risk_gammahat <- estimate_autocov_risk(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = s0, t = t0, lag = 1, bw_grid = bw_grid,
  smooth_ker = epanechnikov, Hs = Hvec_s, Ls = Lvec_s,
  Ht = Hvec_t, Lt = Lvec_t, Delta = NULL, h = presmooth_bw)

### C++ implementation
res_gammahat_risk_cpp <- estimate_autocov_risk_cpp(
  data = data_prepared, s = s0, t = t0, lag = 1, bw_grid = bw_grid,
  use_same_bw = TRUE, center = FALSE,
  kernel_name = "epanechnikov")

dt_gmmahat_risk_cpp <- data.table::as.data.table(res_gammahat_risk_cpp)
names(dt_gmmahat_risk_cpp) <- c("s", "t", "h", "PN_cpp", "bias_term_cpp", "variance_term_cpp", "dependence_term_cpp", "autocov_risk_cpp")
dt_gammahat_risk_all <- data.table::merge.data.table(x = dt_risk_gammahat, y = dt_gmmahat_risk_cpp, by = c("s", "t", "h"))

all.equal(target = dt_gammahat_risk_all[, bias_term_cpp], current = dt_gammahat_risk_all[, bias_term])
all.equal(target = dt_gammahat_risk_all[, variance_term_cpp], current = dt_gammahat_risk_all[, varriance_term])
all.equal(target = dt_gammahat_risk_all[, dependence_term_cpp], current = dt_gammahat_risk_all[, dependence_term])

## Benchmark
res_bm_autocov_risk <- microbenchmark::microbenchmark(
  autocov_risk_cur = dt_risk_gammahat <- estimate_autocov_risk(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = s0, t = t0, lag = 1, bw_grid = bw_grid,
    smooth_ker = epanechnikov, Hs = Hvec_s, Ls = Lvec_s,
    Ht = Hvec_t, Lt = Lvec_t, Delta = NULL, h = presmooth_bw),
  autocov_risk_cpp = res_gammahat_risk_cpp <- estimate_autocov_risk_cpp(
    data = data_prepared, s = s0, t = t0, lag = 1, bw_grid = bw_grid,
    use_same_bw = TRUE, center = FALSE, kernel_name = "epanechnikov"),
  times = 60
)
print(res_bm_autocov_risk, unit = "relative", order = "median")

## two bandwidth
### Current implementation
dt_risk_gammahat_2bw <- estimate_autocov_risk_2bw(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = s0, t = t0, lag = 1, bw_grid = NULL,
  smooth_ker = epanechnikov, Hs = Hvec_s, Ls = Lvec_s,
  Ht = Hvec_t, Lt = Lvec_t, Delta = NULL, h = presmooth_bw)

### C++ implementation
res_gammahat_risk_2bw_cpp <- estimate_autocov_risk_cpp(
  data = data_prepared, s = s0, t = t0, lag = 1, bw_grid = NULL,
  use_same_bw = FALSE, center = TRUE, kernel_name = "epanechnikov")

dt_gmmahat_risk_2bw_cpp <- data.table::as.data.table(res_gammahat_risk_2bw_cpp)
names(dt_gmmahat_risk_2bw_cpp) <- c("s", "t", "hs", "ht", "PN_cpp", "presmooth_bw", "Hs", "Ls", "Ht", "Lt", "bias_term_cpp", "variance_term_cpp", "dependence_term_cpp", "autocov_risk_cpp")
dt_gammahat_risk_all_2bw <- data.table::merge.data.table(x = dt_risk_gammahat_2bw, y = dt_gmmahat_risk_2bw_cpp, by = c("s", "t", "hs", "ht"))

all.equal(target = dt_gammahat_risk_all_2bw[, bias_term_cpp], current = dt_gammahat_risk_all_2bw[, bias_term])
all.equal(target = dt_gammahat_risk_all_2bw[, variance_term_cpp], current = dt_gammahat_risk_all_2bw[, varriance_term])
all.equal(target = dt_gammahat_risk_all_2bw[, dependence_term_cpp], current = dt_gammahat_risk_all_2bw[, dependence_term])

## Graph
dygraphs::dygraph(dt_gammahat_risk_all[t == t0[1], .(h, autocov_risk, autocov_risk_cpp)])


## Benchmark
res_bm_mean_risk <- microbenchmark::microbenchmark(
  autocov_risk_2bw_cur = dt_risk_gammahat_2bw <- estimate_autocov_risk_2bw(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    s = s0, t = t0, lag = 1, bw_grid = bw_grid,
    smooth_ker = epanechnikov, Hs = Hvec_s, Ls = Lvec_s,
    Ht = Hvec_t, Lt = Lvec_t, Delta = NULL, h = presmooth_bw),
  autocov_risk_2bw_cpp = res_gammahat_risk_2bw_cpp <- estimate_autocov_risk_cpp(
    data = data_prepared, s = s0, t = t0, lag = 1, bw_grid = bw_grid,
    Hs = Hvec_s, Ls = Lvec_s, Ht = Hvec_t, Lt = Lvec_t,
    use_same_bw = FALSE, center = TRUE,
    h = presmooth_bw,
    kernel_name = "epanechnikov"),
  times = 50
)
print(res_bm_mean_risk, unit = "relative", order = "median")

## Estimate autocovariance function ----

estimate_autocov(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
                 s = s0, t = t0, lag = 1,
                 optbw = NULL,
                 bw_grid = seq(0.005, 0.15, len = 45),
                 Hs = NULL, Ls = NULL,
                 Ht = NULL, Lt = NULL,
                 Delta = NULL, h = NULL,
                 center = TRUE,
                 mean_estimates_s = NULL,
                 mean_estimates_t = NULL,
                 smooth_ker = epanechnikov)

Rcpp::sourceCpp("./src/05_estimate_autocov_cpp.cpp")
# t0 <- c(0.2, 0.4, 0.5, 0.7, 0.8)
t0 <- seq(0.1, 0.9, len = 200)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]
rm(dt_st) ; gc()

## Covariance function
microbenchmark::microbenchmark(
  "cov" = res_cov <- estimate_autocov_cpp(
    data = data_prepared, s = s0, t = t0, lag = 0,
    param_grid = c(0.2, 0.4, 0.5, 0.7, 0.8),
    use_same_bw = TRUE, center = TRUE, optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
    kernel_name = "epanechnikov"),
  times = 50
)
res_cov <- estimate_autocov_cpp(
  data = data_prepared, s = s0, t = t0, lag = 0,
  param_grid = c(0.2, 0.4, 0.5, 0.7, 0.8),
  use_same_bw = TRUE, center = FALSE, optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
  kernel_name = "epanechnikov")
dt_cov <- data.table::as.data.table(res_cov)
res_mat = reshape_matrix(A = res_cov, idx_col_s = 0, idx_col_t = 1, idx_col_value = 13)

res_gammahat_risk_2bw_cpp <- estimate_autocov_risk_cpp(
  data = data_prepared, s = s0, t = t0, lag = 1, bw_grid = NULL,
  use_same_bw = TRUE, center = TRUE, kernel_name = "epanechnikov")

res_gammahat_risk_2bw_cpp

library(plotly)
library(ggplot2)
library(reshape2)
df <- melt(volcano)

ggplot(dt_cov, aes(x = V1, y = V2, z = V14)) +
  geom_contour_filled(bins = 6) +
  xlab(label = "s") +
  ylab(label = "t") +
  labs(fill = expression(C(s,t))) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ggplot(dt_cov, aes(V1, V2, z= V10, colour=stat(level))) +
  geom_contour() +
  scale_colour_distiller(palette = "YlGn", direction = 1)

ggplotly(p)

## Plotly
dt_cov_mat <- data.table::dcast(data = dt_cov[order(V1,V2)], formula = V1 ~ V2, value.var = "V14")

plot_ly(x = dt_cov_mat[, 1], y = dt_cov_mat[, 1], z = as.matrix(dt_cov_mat[, -1])) %>% add_surface()

all.equal(target = as.matrix(dt_cov_mat[, -1])[, 20], current = res_mat[, 20])

