# Mont2-Carlo simulation for the paper : data generation

# Simulation for data in inst/data ----
## Simulation parameters----
N <- c(100L, 200L, 400L)
lambda <- c(25L, 50L, 100L, 200L, 300L)
sig <- 0.5
mc <- 500
t0 <- seq(0.2, 0.8, len = 6)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

## Simulation function ----
sim_fun <- function(mc_i, Ni, lbda, t0, sig = 0.5,
                    kernel_far = get_real_data_far_kenel,
                    mean_far = get_real_data_mean, hurst = Hlogistic){

  # Generate FTS
  dt_far <- simulate_far(N = Ni, lambda = lbda,
                         tdesign = "random",
                         Mdistribution = rpois,
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

## Simulattion loop ----
dt_N <- data.table::rbindlist(lapply(N, function(Ni){
  dt_lambda <- data.table::rbindlist(lapply(lambda, function(lambdai){
    dt_mc <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
      dt_ <- sim_fun(mc_i = mc_i, Ni = Ni, lbda = lambdai, t0 = t0, sig = sig,
                     kernel_far = get_real_data_far_kenel,
                     mean_far = get_real_data_mean, hurst = Hlogistic)
      return(dt_)
    }, mc.cores = 75))
    rds_name <- paste0("./inst/07_mc_simulate_data/data/dt_mc_fts_real_data_N=", Ni, "_lambda=", lambdai, "_Hlogistic_sig05.RDS")
    saveRDS(object = dt_mc, file = rds_name)
  }))
}))

