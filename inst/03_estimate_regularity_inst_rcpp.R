source("./R/02_smoothing.R")
source("./R/03_estimate_regularity.R")
Rcpp::sourceCpp("./src/03_estimate_locreg_rcpp.cpp")

# Import the data
dt <- readRDS("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=1_fts_model_2.RDS")
dt <- dt[ttag == "trandom"]

# Estimation parameters
t0 <- seq(0.1, 0.9, len = 10)
presmooth_bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
lambdahat <- mean(dt[, .N, by = "id_curve"][, N])
Deltahat <- 2 * exp(-log(lambdahat) ** 0.72)

# Estimation using current function
dt_locreg <- estimate_locreg(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, Delta = Deltahat, h = presmooth_bw,
  smooth_ker = epanechnikov, center = TRUE)

# Estimation using Rcpp
hsmooth <- rep(median(presmooth_bw), 150)
data_prepared <- .format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
dt_locreg_cpp <- estimate_locreg_cpp(
  data = data_prepared,
  t = t0, Delta = Deltahat, h = presmooth_bw,
  kernel_name = "epanechnikov", center = TRUE)

all.equal(target = dt_locreg[, Ht], current = dt_locreg_cpp[, 5])
all.equal(target = dt_locreg[, Lt], current = dt_locreg_cpp[, 6])

plot(x = dt_locreg[, t], y = dt_locreg[, Ht])
plot(x = dt_locreg_cpp[, 1], y = dt_locreg_cpp[, 4])

# Estimate local regularity wrap function

estimate_locreg_wrap <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                 kernel_name = "epanechnikov",
                                 t = 1/2, Delta = NULL, h = NULL, center = TRUE){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  # if (! methods::is(smooth_ker, "function"))
  #   stop("'smooth_ker' must be a function.")
  if (! methods::is(center, "logical"))
    stop("'center' must be a TRUE or FALSE.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Control and set arguments depending on the data
  if (! is.null(Delta)) {
    if (! (methods::is(Delta, "numeric") & data.table::between(Delta, 0, 1) & length(Delta) == 1))
      stop("'Delta' must be a numeric scalar value between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    Delta <- 2 * exp(-log(lambdahat) ** 0.72)
  }

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, smooth_ker = smooth_ker)
      h <- dt_optbw[, median(optbw)]
      rm(dt_optbw) ; gc()
    } else {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, smooth_ker = smooth_ker)
      h <- dt_optbw[, optbw]
      rm(dt_optbw) ; gc()
    }
  }
  dt_res_cpp <- estimate_locreg_cpp(data = data, t = t, Delta = Delta, h = h,
                                    kernel_name = kernel_name, center = center)
  dt_res_cpp <- data.table::as.data.table(dt_res_cpp)
  data.table::setnames(x = dt_res_cpp, new = c("t", "Delta", "Nused", "locreg_bw", "Ht", "Lt"))
  return(dt_res_cpp)
}

## Test the wrapped function
### Current function
dt_locreg_cur <- estimate_locreg(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, Delta = Deltahat, h = presmooth_bw,
  smooth_ker = epanechnikov, center = TRUE)

### cpp function
data_prepared <- .format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
dt_locreg_cpp <- estimate_locreg_cpp(
  data = data_prepared,
  t = t0, Delta = Deltahat, h = presmooth_bw,
  kernel_name = "epanechnikov", center = TRUE)

### Wrapped cpp function
dt_locreg_cpp_wrap <- estimate_locreg_wrap(
  data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = t0, Delta = Deltahat, h = presmooth_bw,
  kernel_name = "epanechnikov", center = TRUE)

### Microbenchmark
res_benchmark <- microbenchmark::microbenchmark(
  dt_locreg_cur = dt_locreg_cur <- estimate_locreg(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t0, Delta = Deltahat, h = presmooth_bw,
    smooth_ker = epanechnikov, center = TRUE),
  dt_locreg_cpp = dt_locreg_cpp <- estimate_locreg_cpp(
    data = data_prepared,
    t = t0, Delta = Deltahat, h = presmooth_bw,
    kernel_name = "epanechnikov", center = TRUE),
  dt_locreg_cpp_wrap = dt_locreg_cpp_wrap <- estimate_locreg_wrap(
    data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
    t = t0, Delta = Deltahat, h = presmooth_bw,
    kernel_name = "epanechnikov", center = TRUE),
  times = 50
)

print(res_benchmark, unit = "relative", order = "median")

# Permit Null values
estimate_locreg_cpp(
  data = data_prepared,
  t = t0, Delta = NULL, h = NULL,
  kernel_name = "epanechnikov", center = TRUE)
