library(data.table)
library(parallel)
library(matrixStats)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Mean function estimation function ----
## Estimate mean functions risk
estim_mean_risk_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, n_cores = 75){

  # Load local regularity estimates
  locreg_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                             process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(locreg_file_name)

  # Define bandwidth grid
  K <- 40
  b0 <- 4 * (N * lambda) ** (- 0.9)
  # b0 <- max(b0, 1 / 1000)
  # bK <- 4 * (N * lambda) ** (- 1 / 3)
  bK <- 5 / 10
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()

  if (white_noise == "mfBm") {
    # Estimate mean risk by mc
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)
    index_mc <- id_mc_data

    # Get the already done MC repettion
    # Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda)))
    # if (Nmc_done == 0) {
    #   index_mc <- sort(id_mc_data)
    # } else if (Nmc_done < length(id_mc_data)) {
    #   data_file_done <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda))
    #   id_mc_done <- gsub(pattern = paste0("dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", "|", "_", design,".RDS"),
    #                      replacement = "", x = data_file_done)
    #   id_mc_done <- as.numeric(id_mc_done)
    #   index_mc <- setdiff(id_mc_data, id_mc_done)
    #   index_mc <- sort(index_mc)
    # }
    res <- parallel::mclapply(index_mc, function(mc_i, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Load data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_random_mc <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]

      # Estimate the risk of the mean function
      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg[id_mc == mc_i & order(t), Ht]
      Lt <- dt_locreg[id_mc == mc_i & order(t), Lt]

      start_time <- Sys.time()
      dt_risk_muhat <- estimate_mean_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = 'tobs', ycol = "X",
        t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
      )
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + dt_locreg[id_mc == mc_i, unique(time_taken)]

      # Save
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "time_taken" = time_taken, dt_risk_muhat)
      rm(dt_risk_muhat) ; gc()

      file_save <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda, "/dt_mean_risk_",
                          process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      saveRDS(object = dt_res, file = file_save)
      rm(dt_res) ; gc()
      return(paste0("Mean risk : done for MC = ", mc_i, ", N = ", N, " and lambda = ", lambda))
    }, mc.cores = n_cores, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0)

  } else if (white_noise == "fBm") {
    # Not need for the paper ...
  }

  rm(locreg_file_name, dt_locreg, res) ; gc() ; gc()

  return(paste0("Done : dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function only for "plus_mean"
estim_mean_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0, n_cores = 75){
  # Load mean risk data
  mean_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                                process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_bis.RDS")
  dt_mean_risk <- readRDS(mean_risk_file_name)

  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)

    # Get the already done MC repettion
    Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda)))
    if (Nmc_done < length(id_mc_data)) {
      data_file_done <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda))
      id_mc_done <- gsub(pattern = paste0("dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", "|", "_", design,".RDS"),
                         replacement = "", x = data_file_done)
      id_mc_done <- as.numeric(id_mc_done)
      index_mc <- id_mc_done
    } else {
      index_mc <- id_mc_data
    }

    # Get optimal bandwidth smoothing parameter
    dt_optbw <- dt_mean_risk[, .("optbw" = h[which.min(mean_risk)], "time_taken" = unique(time_taken)), by = c("id_mc", "t")]
    # Estimate mean by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_optbw, Ni, lambdai, t0){
      # Load data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)

      # Extract data
      dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
      t0 <- dt_optbw[id_mc == mc_i][order(t), t]

      # Estimate the mean function
      start_time <- Sys.time()
      dt_mean <- estimate_mean(
        data = dt_random_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, optbw = optbw)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + dt_optbw[id_mc == mc_i, unique(time_taken)]

      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai,
                                       "time_taken" = time_taken, dt_mean[, .(t, optbw, PN, muhat)])
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
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_bis.RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(mean_risk_file_name, file_name, dt_mean_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_mean_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean risk function ----
## design 1 ----

estim_mean_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, n_cores = 30)


estim_mean_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)

## design 3 ----
estim_mean_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, n_cores = 30)
estim_mean_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, n_cores = 30)


estim_mean_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)

## design 5 ----
estim_mean_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, n_cores = 35)
estim_mean_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, n_cores = 35)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, n_cores = 35)
estim_mean_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, n_cores = 35)


estim_mean_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", n_cores = 30)
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d5", n_cores = 30)
estim_mean_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d5", n_cores = 30)



