library(data.table)

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# Simulation function ----
simulate_data_fun <- function(mc_i, N, lambda, t0, sig = 0.25,
                              process_ker = get_real_data_far_kenel,
                              process_mean = get_real_data_mean,
                              Lt = 4, hurst = Hlogistic, design = "d5"){
  # Generate data according the choixe of the white noise

  ### Simulate mean-zero FAR(1) process
  dt_gen <- simulate_far(N = N, lambda = lambda,
                         tdesign = "random",
                         Mdistribution = bounded_uniform,
                         tdistribution = runif,
                         tcommon = t0,
                         hurst_fun = hurst,
                         L = Lt,
                         far_kernel = process_ker,
                         far_mean = process_mean,
                         int_grid = 100L,
                         burnin = 100L,
                         remove_burnin = TRUE)
  dt_gen[ttag == "tcommon", mean(X), by = "tobs"]
  ### Add mean function and noise
  dt_gen[, process_mean := far_mean]
  dt_gen[, far_mean := NULL]
  dt_gen[ttag == "trandom", X := X + rnorm(n = .N, mean = 0, sd = sig), by = "id_curve"]
  ### Get pre-smoothing bandwidth
  #### Define and exponential bandwidth grid
  lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
  K <- 30
  b0 <- 1 / lambdahat
  bK <- lambdahat ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(b0, bK, a, K) ; gc()

  #### Get optimal bw for each curve
  index <- dt_gen[, unique(id_curve)]
  best_bw <- median(unlist(lapply(tail(index, 20), function(id, dtt, bw_grid){
    # Filter data
    d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]
    # Get optimal bandwidth
    bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
    rm(d) ; gc()
    return(bw)
  }, dtt = dt_gen, bw_grid = bw_grid)))
  dt_gen[, presmooth_bw := best_bw]

  # Add MC index
  dt_gen[, c("id_mc", "N", "lambda") := .(mc_i, N, lambda)]
  data.table::setcolorder(
    x = dt_gen,
    neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "process_mean", "X", "presmooth_bw")
  )
  save_file_name <- paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda,
                           "/dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", mc_i, "_", design, ".RDS")
  saveRDS(dt_gen, save_file_name)
  rm(dt_gen) ; gc()

  # return the result
  return("Done !")
}

# Simulate all MC process

simulate_data <- function(index_mc = 1:400, N = 400, lambda = 300, t0, sig = 0.25,
                          process_ker = get_real_data_far_kenel,
                          process_mean = get_real_data_mean,
                          Lt = 4, hurst = Hlogistic, design = "d5"){
  # Not redo what is already done
  file_done <- list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda))
  id_mc_data <- gsub(pattern = paste0("dt_mc_FAR_mfBm_N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                     replacement = "", x = file_done)
  id_mc_data <- as.numeric(id_mc_data)
  index_mc <- setdiff(index_mc, id_mc_data)

  # Run in parallel
  parallel::mclapply(index_mc, function(mc_i, N, lambda, t0, sig,
                                        process_ker, process_mean, Lt, hurst, design){
    dt_ <- simulate_data_fun(mc_i = mc_i, N = N, lambda = lambda, t0 = t0, sig = sig,
                             process_ker = process_ker, process_mean = process_mean,
                             Lt = Lt, hurst = hurst, design = design)
    return(dt_)
  }, mc.cores = 30, N = N, lambda = lambda, t0 = t0, sig = sig,
  process_ker = process_ker, process_mean = process_mean, Lt = Lt, hurst = hurst, design)

  print(paste0("Done : dt_mc_FAR_mfBm_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Design 5 ----
## Mean function
dt_mean <- readRDS("./inst/12_mc_simulate_data/real_data/dt_emp_mean.RDS")
get_emp_mean_interpolate <- function(t, dt_mean = dt_mean){
  res <- approx(x = dt_mean[, tobs], y = dt_mean[, muhat], xout = t)$y
  return(res)
}
mean_d5 <- function(t) get_emp_mean_interpolate(t = t, dt_mean = dt_mean)

## Autoregressive kernel
ker_d1 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

# Clust 8
simulate_data(index_mc = 1:400, N = 150, lambda = 40, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")

simulate_data(index_mc = 1:400, N = 1000, lambda = 40, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")

simulate_data(index_mc = 1:400, N = 400, lambda = 300, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")

# Clust 8 mc = 100
simulate_data(index_mc = 1:100, N = 1000, lambda = 1000, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")
# Clust 6 mc = 200
simulate_data(index_mc = 101:300, N = 1000, lambda = 1000, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")
# Clust 6 mc = 100
simulate_data(index_mc = 301:400, N = 1000, lambda = 1000, t0, sig = 0.25,
              process_ker = ker_d1, process_mean = mean_d5,
              Lt = 4, hurst = Hlogistic, design = "d5")
