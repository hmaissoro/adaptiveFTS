# Mont2-Carlo simulation for the paper : data generation

## Simulation global parameters----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}
## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  round(runif(n = N, min = lambda * (1 -p), max = lambda * (1 + p)))
}

## Simulation function ----
sim_fun <- function(mc_i, Ni, lbda, t0, sig = 0.5,
                    kernel_far = get_real_data_far_kenel,
                    mean_far = get_real_data_mean, hurst = Hlogistic){

  # Generate FTS
  dt_far <- simulate_far(N = Ni, lambda = lbda,
                         tdesign = "random",
                         Mdistribution = bounded_uniform,
                         tdistribution = runif,
                         tcommon = t0,
                         hurst_fun = hurst,
                         L = 4,
                         far_kernel = kernel_far,
                         far_mean = mean_far,
                         int_grid = 100L,
                         burnin = 100L,
                         remove_burnin = TRUE)

  # Get pre-smoothing bandwidth
  ## Define and exponential bandwidth grid
  lambdahat <- mean(dt_far[ttag == "trandom", .N, by = id_curve][, N])
  K <- 100
  b0 <- 1 / lambdahat
  bK <- lambdahat ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))

  ## Determine bw for each curve
  index <- dt_far[, unique(id_curve)]
  dt <- data.table::rbindlist(lapply(index, function(id, dtt, bw_grid){
    # Filter data
    d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]

    # Add noise
    d[, X := X + rnorm(n = .N, mean = 0, sd = sig)]

    # Get be
    bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
    d[, presmooth_bw := bw]
    return(d)
  }, dtt = dt_far, bw_grid = bw_grid))

  # Add fix design points
  dt_tcommon <- data.table::merge.data.table(
    x = dt_far[ttag == "tcommon"],
    y = unique(dt[, .(id_curve, presmooth_bw)]),
    by = "id_curve")

  dt_res <- rbind(dt, dt_tcommon)

  # Clean
  rm(dt, dt_tcommon, index, K, dt_far, b0, bK, a, bw_grid) ; gc()

  # Add MC index
  dt_res[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lbda)]
  data.table::setcolorder(
    x = dt_res,
    neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "far_mean", "X", "presmooth_bw"))
  return(dt_res)
}

# Simulation - design 1 ----

## Mean function
far_mean_d1 <- function(t){
  res <- 4 * sin(3 * pi * t / 2)
  return(res)
}

## Autoregressive kernel
far_ker_d1 <- function(s,t, operator_norm = 0.7){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_nom / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

## Data generation
### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d1 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_far = far_ker_d1, mean_far = far_mean_d1, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d1, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=400_lambda=300_d1.RDS")

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d1 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_far = far_ker_d1, mean_far = far_mean_d1, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d1, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=1000_lambda=1000_d1.RDS")

# Simulation - design 2 ----
## Mean function from real data + Same kernel as design 1

## Data generation
### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d2 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_far = far_ker_d1, mean_far = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d2, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=400_lambda=300_d2.RDS")

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d2 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_far = far_ker_d1, mean_far = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d2, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=1000_lambda=1000_d2.RDS")

# Simulation - design 3 ----
## Mean and far kenel estimated from real data

## Data generation
### (N, lambda) = (400, 1000)
dt_mc_N400_lambda300_d3 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0, sig = sig,
                 kernel_far = get_real_data_far_kenel, mean_far = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N400_lambda300_d3, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=400_lambda=300_d3.RDS")

### (N, lambda) = (1000, 1000)
dt_mc_N1000_lambda1000_d3 <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_fun(mc_i = mc_i, Ni = 1000, lbda = 1000, t0 = t0, sig = sig,
                 kernel_far = get_real_data_far_kenel, mean_far = get_real_data_mean, hurst = Hlogistic)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_N1000_lambda1000_d3, file = "./inst/08_mc_simulate_data/data/dt_mc_far_N=1000_lambda=1000_d3.RDS")

