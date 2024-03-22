library(data.table)
library(foreach)
library(parallel)
library(doParallel)

# Simulation global parameters----


# Autocovariance optimal function ----
estim_autocov_bw_rp_by_mc <- function(mc_i = 1, N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  # Load raw data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
  optbw_file <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda,"/dt_mean_optbw_rp_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")

  dt <- readRDS(data_file_name)
  dt_optbw <- readRDS(optbw_file)
  if (white_noise == "mfBm") {
    # The bandiwdth R&P bandwidth for all mc
    K <- 20
    b0 <- N ** (- 5 / 4)
    bK <- N ** (- 1 / 4)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a) ; gc()

    start_time <- Sys.time()
    # Estimate the mean function at each tobs
    dt_mean_rp <- dt[order(tobs), list(id_curve, tobs)]
    optbw_rp <- dt_optbw[id_mc == mc_i, optbw_rp]
    dt_mean <- estimate_mean_rp(
      data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = dt_mean_rp[, tobs], h = optbw_rp, smooth_ker = epanechnikov)
    dt_mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
    rm(dt_mean) ; gc()

    system.time({
      dt_autocov_bw_rp <- estimate_autocov_bw_rp(
        data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
        Kfold = 10, bw_grid = bw_grid, optbw_mean = NULL,
        dt_mean_rp = dt_mean_rp, smooth_ker = epanechnikov)
    })
    end_time <- Sys.time()
    time_taken <- end_time - start_time

    # Return and clean
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda,
                                     "optbw_rp" = dt_autocov_bw_rp[, h[which.min(cv_error)]],
                                     "time_taken" = as.numeric(time_taken))
    rm(dt_autocov_bw_rp, dt) ; gc()

    file_save <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/autocov/", design, "_optbw/N", N, "lambda", lambda,"/dt_autocov_optbw_rp_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
    saveRDS(object = dt_res, file = file_save)
  } else if (white_noise == "fBm") {
    # Coming ...
  }

  return(paste0("Done : dt_autocov_optbw_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS at ", Sys.time()))
}

estim_autocov_bw_rp <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3){
  source("./R/03_estimate_regularity.R")
  source("./R/04_estimate_mean.R")
  source("./R/05_estimate_autocovariance.R")
  Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda)))
  Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/autocov/", design, "_optbw/N", N, "lambda", lambda)))

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  foreach::foreach(i = (Nmc_done + 1):Nmc, .inorder=FALSE,
                   .packages = c("adaptiveFTS", "data.table", "fastmatrix"),
                   .export = c("estim_autocov_bw_rp_by_mc", "estimate_mean_rp", "estimate_autocov_bw_rp", "estimate_autocov_rp", ".format_data"),
                   .verbose = TRUE) %dopar% {
                     estim_autocov_bw_rp_by_mc(mc_i = i, N = N, lambda = lambda, process = process, white_noise = white_noise, design = design)
                   }
  parallel::stopCluster(cl)
}
# Estimate rubin and panaretos bw ----
# Scenario 1 : ----

## clust-n11
estim_autocov_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)
estim_autocov_bw_rp(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)

## clust-n9
estim_autocov_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 6)

## clust-n13
estim_autocov_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 6)

## clust-n10
estim_autocov_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 7)

# Scenario 3 : ----
## clust-n11
estim_autocov_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 3)
estim_autocov_bw_rp(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 3)

## clust-n9
estim_autocov_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 6)

## clust-n13
estim_autocov_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 6)

## clust-n10
estim_autocov_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 7)

# Autocovariance RP function ----

estim_autocov_rp_fun <- function(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1){
  # Load raw data

  if (white_noise == "mfBm") {
    Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda)))
    Nmc_optbw <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/rubin_method/autocov/", design, "_optbw/N", N, "lambda", lambda)))
    if (Nmc_optbw < Nmc) Nmc <- Nmc_optbw
    # Estimate local regularity by mc
    dt_mean_rp_mc <- data.table::rbindlist(parallel::mclapply(1:Nmc, function(mc_i, N, lambda, t0){
      # Extract data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")

      dt_optbw_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/autocov/", design, "_optbw/N", N, "lambda", lambda,"/dt_autocov_optbw_rp_",
                                   process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt_optbw_mean_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda,"/dt_mean_optbw_rp_",
                                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_optbw_autocov <- readRDS(dt_optbw_file_name)
      dt_optbw_mean <- readRDS(dt_optbw_mean_file_name)

      optbw_rp_autocov <- dt_optbw_autocov[id_mc == mc_i, optbw_rp]
      optbw_rp_mean <- dt_optbw_mean[id_mc == mc_i, optbw_rp]

      start_time <- Sys.time()
      # Estimate the mean function at each tobs
      dt_mean_rp <- dt[order(tobs), list(id_curve, tobs)]
      dt_mean <- estimate_mean_rp(
        data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = dt_mean_rp[, tobs], h = optbw_rp_mean, smooth_ker = epanechnikov)
      dt_mean_rp[, muhat_RP := dt_mean[, muhat_RP]]
      rm(dt_mean) ; gc()

      # Estimate the autocovariance
      dt_autocov_rp <- estimate_autocov_rp(
        data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = s0, t = t0, lag = lag, h = optbw_rp_autocov, optbw_mean = NULL,
        dt_mean_rp = dt_mean_rp, smooth_ker = epanechnikov)
      end_time <- Sys.time()

      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + optbw_rp_autocov[id_mc == mc_i, time_taken]
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda, "time_taken" = time_taken, dt_autocov_rp)
      rm(dt_mean_rp, dt, optbw_rp_autocov, optbw_rp_mean, dt_optbw_mean, dt_optbw_autocov, data_file_name,
         dt_optbw_file_name, dt_optbw_mean_file_name) ; gc()
      return(dt_res)
    }, mc.cores = 75, N = N, lambda = lambda, t0 = t0))

  } else if (white_noise == "fBm") {
    # Comming ...
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/autocov/dt_autocov_rp_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_rp_mc, file = file_name)
  rm(file_name) ; gc() ; gc()

  return(paste0("Done : dt_mean_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

estim_autocov_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_rp_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_rp_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_rp_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)

# Scenario 3 :
estim_autocov_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1)
estim_autocov_rp_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1)
estim_autocov_rp_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1)



