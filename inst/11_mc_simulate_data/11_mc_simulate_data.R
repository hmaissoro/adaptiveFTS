# Mont2-Carlo simulation for the paper : data generation

library(data.table)

## Simulation global parameters----
sig <- 0.5
mc <- 75
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}
## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# Functions ----
zero_mean_func <- function(t) 0 * t

## Simulation function ----
## If center = TRUE, the mean function is added to a mean-zero FAR(1) process
sim_fun <- function(mc_i, Ni, lbda, t0, sig = 0.5,
                    kernel_fma = get_real_data_far_kenel,
                    mean_fma = get_real_data_mean,
                    hurst = Hlogistic, center = FALSE){

  # Generate FTS
  if (center) {
    # Simulate mean-zero FAR(1) process
    dt_fma <- simulate_fma(N = Ni, lambda = lbda,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = t0,
                           hurst_fun = hurst,
                           L = 4,
                           fma_kernel = kernel_fma,
                           fma_mean = zero_mean_func,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
  } else {
    dt_fma <- simulate_far(N = Ni, lambda = lbda,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = t0,
                           hurst_fun = hurst,
                           L = 4,
                           far_kernel = kernel_fma,
                           far_mean = zero_mean_func,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
    # add mean_fma
    dt_fma[, X := X + mean_fma(tobs)]
  }

  # Get pre-smoothing bandwidth
  ## Define and exponential bandwidth grid
  lambdahat <- mean(dt_fma[ttag == "trandom", .N, by = id_curve][, N])
  K <- 100
  b0 <- 1 / lambdahat
  bK <- lambdahat ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))

  ## Determine bw for each curve
  index <- dt_fma[, unique(id_curve)]
  dt <- data.table::rbindlist(lapply(index, function(id, dtt, bw_grid){
    # Filter data
    d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]

    # Add noise
    d[, X := X + rnorm(n = .N, mean = 0, sd = sig)]

    # Get be
    bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
    d[, presmooth_bw := bw]
    return(d)
  }, dtt = dt_fma, bw_grid = bw_grid))

  # Add fix design points
  dt_tcommon <- data.table::merge.data.table(
    x = dt_fma[ttag == "tcommon"],
    y = unique(dt[, .(id_curve, presmooth_bw)]),
    by = "id_curve")

  dt_res <- rbind(dt, dt_tcommon)

  # Clean
  rm(dt, dt_tcommon, index, K, dt_fma, b0, bK, a, bw_grid) ; gc()

  # Add MC index
  dt_res[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lbda)]
  data.table::setcolorder(
    x = dt_res,
    neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "fma_mean", "X", "presmooth_bw"))
  return(dt_res)
}

# Simulation - design 1 ----
## Mean function
far_mean_d1 <- function(t) 4 * sin(3 * pi * t / 2)

## Autoregressive kernel
far_ker_d1 <- function(s,t, operator_norm = 0.7){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d1 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_fma = far_ker_d1, mean_fma = far_mean_d1, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d1, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=400_lambda=300_d1.RDS")
rm(dt_mc_N400_lambda300_d1) ; gc()

# Centered version
dt_mc_N400_lambda300_d1 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_fma = far_ker_d1, mean_fma = far_mean_d1, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d1, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=400_lambda=300_centered_d1.RDS")
rm(dt_mc_N400_lambda300_d1) ; gc()

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d1 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_fma = far_ker_d1, mean_fma = far_mean_d1, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d1, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=1000_lambda=1000_d1.RDS")
rm(dt_mc_N1000_lambda1000_d1) ; gc()


# Simulation - design 2 ----
## Mean function from real data + Same kernel as design 1

### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d2 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_fma = far_ker_d1, mean_fma = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d2, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=400_lambda=300_d2.RDS")

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d2 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_fma = far_ker_d1, mean_fma = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d2, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=1000_lambda=1000_d2.RDS")
rm(dt_mc_N1000_lambda1000_d2) ; gc()

# Simulation - design 3 ----
## Mean and far kenel estimated from real data
### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d3 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_fma = get_real_data_far_kenel, mean_fma = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 77))
saveRDS(object = dt_mc_N400_lambda300_d3, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=400_lambda=300_d3.RDS")
rm(dt_mc_N400_lambda300_d3) ; gc()

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d3 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_fma = get_real_data_far_kenel, mean_fma = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d3, file = "./inst/11_mc_simulate_data/data/dt_mc_fma_N=1000_lambda=1000_d3.RDS")
rm(dt_mc_N1000_lambda1000_d3) ; gc()

