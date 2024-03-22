library(data.table)
library(parallel)
library(matrixStats)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Mean function estimation function ----

## Estimate mean function only for "plus_mean"
estim_mean_fun_RP_bw <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0, n_cores = 30){
  # Load mean risk data
  file_mean_rp <- paste0("../../matlab/rubin_mean_autocov_method/N", N, "lambda", lambda, design, ".csv")
  dt_mean_rp <- data.table::fread(file_mean_rp)

  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)

    # Get optimal bandwidth smoothing parameter
    dt_optbw <- dt_mean_rp[, .("optbw" = bdth), by = id_mc]
    # Estimate mean by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(id_mc_data, function(mc_i, dt_optbw, Ni, lambdai, t0){
      # Load data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)

      # Extract data
      dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][, optbw]
      optbw <- rep(optbw, 4)
      t0 <- dt[ttag == "tcommon" & id_mc == mc_i][order(tobs), unique(tobs)]
      t0 <- sort(t0)

      # Estimate the mean function
      start_time <- Sys.time()
      dt_mean <- estimate_mean(
        data = dt_random_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, optbw = optbw, smooth_ker = uniform)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      # time_taken <- time_taken + dt_optbw[id_mc == mc_i, unique(time_taken)]

      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_mean[, .(t, optbw, PN, muhat)])
      dt_res <- data.table::merge.data.table(
        x = dt_res,
        y = unique(dt[ttag == "tcommon", .("t" = tobs, "mutrue" = process_mean)]),
        by = "t")
      rm(optbw, dt_mean) ; gc()
      return(dt_res)
    }, mc.cores = n_cores, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    # Not need for the paper
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_rp_bw_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(file_name, dt_mean_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_mean_estimates_RP_bw_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean risk function ----
## design 1 ----

estim_mean_fun_RP_bw(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun_RP_bw(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun_RP_bw(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun_RP_bw(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)

## design 3 ----

estim_mean_fun_RP_bw(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun_RP_bw(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun_RP_bw(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun_RP_bw(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun_RP_bw(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)



