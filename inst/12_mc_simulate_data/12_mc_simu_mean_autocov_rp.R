library(data.table)
library(foreach)
library(parallel)
library(doParallel)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Mean function estimation function ----

## Estimate mean function only for "plus_mean"
estim_mean_bw_rp_by_mc <- function(mc_i = 1, N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  # Load raw data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  if (white_noise == "mfBm") {
    # Estimate R&P bandwidth by mc
    K <- 20
    b0 <- N ** (- 5 / 4)
    bK <- N ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a) ; gc()

    # Estimate the mean function
    start_time <- Sys.time()
    dt_mean_bw_rp <- estimate_mean_bw_rp(
      data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
      Kfold = 10, bw_grid = bw_grid, smooth_ker = epanechnikov)
    end_time <- Sys.time()
    time_taken <- end_time - start_time

    # dt_mean_bw_rp[, h[which.min(cv_error)]]
    # Return and clean
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda,
                                     "optbw_rp" = dt_mean_bw_rp[, h[which.min(cv_error)]],
                                     "time_taken" = as.numeric(time_taken))
    rm(dt_mean_bw_rp, dt) ; gc()

    file_save <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda,"/dt_mean_optbw_rp_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
    saveRDS(object = dt_res, file = file_save)
  } else if (white_noise == "fBm") {
    # Coming ...
  }

  return(paste0("Done : dt_mean_optbw_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS at ", Sys.time()))
}

estim_mean_bw_rp <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3){
  source("./R/04_estimate_mean.R")
  source("./R/03_estimate_regularity.R")
  Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda)))
  Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda)))

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  foreach::foreach(i = (Nmc_done + 1):Nmc, .inorder=FALSE,
                   .packages = c("adaptiveFTS", "data.table", "fastmatrix"),
                   .export = c("estim_mean_bw_rp_by_mc", "estimate_mean_rp", "estimate_mean_bw_rp", ".format_data"),
                   .verbose = TRUE) %dopar% {
                     estim_mean_bw_rp_by_mc(mc_i = i, N = N, lambda = lambda, process = process, white_noise = white_noise, design = design)
                   }
  parallel::stopCluster(cl)
}
# Estimate rubin and panaretos bw ----
# Scenario 1 : ----

## clust-n11
estim_mean_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)
estim_mean_bw_rp(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)

## clust-n9
estim_mean_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 6)

## clust-n13
estim_mean_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 6)

## clust-n10
estim_mean_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 7)


# estim_mean_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)
# estim_mean_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)
# estim_mean_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)
# estim_mean_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 3)

# Scenario 3 : ----
## clust-n11
estim_mean_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 3)
estim_mean_bw_rp(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 3)

## clust-n9
estim_mean_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 6)

## clust-n13
estim_mean_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 6)

## clust-n10
estim_mean_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 7)

# Estimate mean RP function ----

estim_mean_rp_fun <- function(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0){
  # Load raw data

  if (white_noise == "mfBm") {
    Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda)))
    Nmc_optbw <- length(list.files(paste0("./inst/12_mc_simulate_data/FAR/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda)))
    if (Nmc_optbw < Nmc) Nmc <- Nmc_optbw
    # Estimate local regularity by mc
    dt_mean_rp_mc <- data.table::rbindlist(parallel::mclapply(1:Nmc, function(mc_i, N, lambda, t0){
      # Extract data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")

      dt_optbw_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/", design, "_optbw/N", N, "lambda", lambda,"/dt_mean_optbw_rp_",
                                   process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_optbw <- readRDS(dt_optbw_file_name)

      optbw_rp <- dt_optbw[id_mc == mc_i, optbw_rp]
      # Estimate the mean function
      start_time <- Sys.time()
      dt_mean_rp <- estimate_mean_rp(
        data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = t0, h = optbw_rp, smooth_ker = epanechnikov)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + dt_optbw[id_mc == mc_i, time_taken]
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda, "time_taken" = time_taken, dt_mean_rp)
      rm(dt_mean_rp, dt, dt_optbw, data_file_name, dt_optbw_file_name) ; gc()
      return(dt_res)
    }, mc.cores = 75, N = N, lambda = lambda, t0 = t0))

  } else if (white_noise == "fBm") {
    # Comming ...
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/rubin_method/mean/dt_mean_rp_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_rp_mc, file = file_name)
  rm(file_name) ; gc() ; gc()

  return(paste0("Done : dt_mean_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

estim_mean_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)

# Scenario 3 :
estim_mean_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)
estim_mean_rp_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)
estim_mean_rp_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)



