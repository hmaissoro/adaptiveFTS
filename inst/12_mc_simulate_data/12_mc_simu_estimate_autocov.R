library(data.table)
library(foreach)
library(parallel)
library(doParallel)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]

rm(dt_st) ; gc()

# Mean function estimation function ----
## Estimate mean functions risk
estim_autocov_risk_fun <- function(N = 400, lambda = 300, process = "FAR",
                                   white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1, n_cores = 75){
  # Load local regularity estimates
  locreg_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                             process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(locreg_file_name)

  # Define bandwidth grid
  K <- 40
  b0 <- 4 * max((N * lambda) ** (- 0.9), (N * (lambda ** 2)) ** (- 0.9))
  # b0 <- max(b0, 1 / 1000)
  # bK <- 4 * max((N * lambda) ** (- 1 / 3), (N * (lambda ** 2)) ** (- 1 / 3))
  # bK <- min(bK, 10 / 100)
  bK <- 5 / 10
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()

  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)
    index_mc <- sort(id_mc_data)
    # Get the already done MC repettion
    # Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda)))
    # if (Nmc_done == 0) {
    #   index_mc <- sort(id_mc_data)
    # } else if (Nmc_done < length(id_mc_data)) {
    #   data_file_done <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda))
    #   id_mc_done <- gsub(pattern = paste0("dt_autocov_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", "|", "_", design,".RDS"),
    #                      replacement = "", x = data_file_done)
    #   id_mc_done <- as.numeric(id_mc_done)
    #   index_mc <- setdiff(id_mc_data, id_mc_done)
    #   index_mc <- sort(index_mc)
    # }
    # Estimate the aucovariance risk by mc
    res <- parallel::mclapply(index_mc, function(mc_i, dt_locreg, N, lambda, process, white_noise, design){
      # Load data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt_random_mc <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]

      # Estimate the risk of the gamma function ----
      ## Extract local regularity parameters
      dt_locreg_mc <- dt_locreg[id_mc == mc_i]
      dt_locreg_mc <- data.table::merge.data.table(
        x = data.table::data.table("s" = s0, "t" = t0),
        y = dt_locreg[id_mc == mc_i, .("t" = t, Ht, Lt)],
        by = "t"
      )
      dt_locreg_mc <- data.table::merge.data.table(
        x = dt_locreg_mc,
        y = dt_locreg[id_mc == mc_i, .("s" = t, "Hs" = Ht, "Ls" = Lt)],
        by = "s"
      )

      bw <- unique(dt_random_mc[id_mc == mc_i, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg_mc[order(s, t), Ht]
      Lt <- dt_locreg_mc[order(s, t), Lt]
      Hs <- dt_locreg_mc[order(s, t), Hs]
      Ls <- dt_locreg_mc[order(s, t), Ls]

      start_time <- Sys.time()
      dt_risk_gammahat <- estimate_autocov_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = s0, t = t0, lag = 1, bw_grid = bw_grid,
        smooth_ker = epanechnikov, Hs = Hs, Ls = Ls,
        Ht = Ht, Lt = Lt, Delta = NULL, h = bw)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + dt_locreg[id_mc == mc_i, unique(time_taken)]

      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lambda" = lambda, "time_taken" = time_taken, dt_risk_gammahat)
      rm(dt_risk_gammahat, start_time, end_time) ; gc()

      file_save <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda, "/dt_autocov_risk_",
                          process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", mc_i, "_", design,".RDS")
      saveRDS(object = dt_res, file = file_save)
      rm(dt_res) ; gc()
      return(paste0("Autocov risk : done for MC = ", mc_i, ", N = ", N, " and lambda = ", lambda))
    }, mc.cores = n_cores, dt_locreg = dt_locreg, N = N, lambda = lambda, process = process, white_noise = white_noise, design = design)
  } else if (white_noise == "fBm") {
    # Not need for the paper.
  }
  rm(locreg_file_name, dt_locreg, res) ; gc()

  return(paste0("Done : dt_autocov_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function
estim_autocov_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", lag = 1, s0 = s0, t0 = t0, n_cores = 75){
  ## Load mean risk data
  autocov_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocov_risk_",
                                   process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_autocov_risk <- readRDS(autocov_risk_file_name)

  ## Load \widetilde{\gamma}_N(t)
  autocovtilde_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocovtilde_",
                                   process,"_", white_noise, "_", "N=5000_", design,".RDS")
  dt_autocovtilde <- readRDS(autocovtilde_file_name)

  ## Load \widehat \mu_N^*(t) data
  # mean_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
  #                          process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  # dt_mean <- readRDS(mean_file_name)

  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)

    # Get the already done MC repettion
    Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda)))
    if (Nmc_done < length(id_mc_data)) {
      data_file_done <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda))
      id_mc_done <- gsub(pattern = paste0("dt_autocov_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", "|", "_", design,".RDS"),
                         replacement = "", x = data_file_done)
      id_mc_done <- as.numeric(id_mc_done)
      index_mc <- id_mc_done
    } else {
      index_mc <- id_mc_data
    }

    # Estimate optimal bandwidth and autocovariance proxy
    dt_optbw <- unique(dt_autocov_risk[!is.nan(autocov_risk), .("optbw" = h[which.min(autocov_risk)], "taken_time" = unique(time_taken)), by = c("id_mc", "s", "t")])
    dt_optbw <- data.table::merge.data.table(x = dt_optbw, y = dt_autocovtilde, by = c("s", "t"))
    # Add optimal mean function
    ## Add \mu(t)
    dt_optbw[, c("muhat_opt_s", "mutrue_opt_s", "muhat_opt_t", "mutrue_opt_t") := .(0, 0, 0, 0)]
    # dt_optbw <- data.table::merge.data.table(
    #   x = dt_optbw,
    #   y = dt_mean[, .(id_mc, "t" = t, "muhat_opt_t" = muhat, "mutrue_opt_t" = mutrue)],
    #   by = c("id_mc", "t")
    # )
    # ## Add \mu(s)
    # dt_optbw <- data.table::merge.data.table(
    #   x = dt_optbw,
    #   y = dt_mean[, .(id_mc, "s" = t, "muhat_opt_s" = muhat, "mutrue_opt_s" = mutrue)],
    #   by = c("id_mc", "s"),
    #   all.x = TRUE
    # )
    # dt_time_mean <- unique(dt_mean[, .(id_mc, time_taken)])

    # Estimate autocovariance by mc
    dt_autocov_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_optbw, Ni, lambdai, lag, s0, t0){
      ## Load raw data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)

      # Extract data
      dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][order(s, t), optbw]
      muhat_opt_s <- dt_optbw[id_mc == mc_i][order(s, t), muhat_opt_s]
      muhat_opt_t <- dt_optbw[id_mc == mc_i][order(s, t), muhat_opt_t]

      # Estimate the mean function
      start_time <- Sys.time()
      dt_autocov <- estimate_autocov(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = s0, t = t0, lag = lag, optbw = optbw, bw_grid = NULL,
        Hs = NULL, Ls = NULL, Ht = NULL, Lt = NULL,
        Delta = NULL, h = NULL, center = TRUE,
        mean_estimates_s = muhat_opt_s, mean_estimates_t = muhat_opt_t,
        smooth_ker = epanechnikov)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      # time_taken <- time_taken + dt_optbw[id_mc == mc_i][order(s, t), taken_time] + dt_time_mean[id_mc == mc_i, unique(time_taken)]

      # Add estimate of the true gamma
      dt_autocov <- data.table::merge.data.table(
        x = dt_autocov,
        y = dt_optbw[id_mc == mc_i, .(s, t, mutilde_s, mutilde_t, gammatilde, autocovtilde)],
        by = c("s", "t")
      )
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "time_taken" = time_taken, dt_autocov)
      # rm(optbw) ; gc()
      return(dt_res)
    }, mc.cores = n_cores, dt_optbw = dt_optbw, Ni = N, lambdai = lambda, lag = lag, s0 = s0, t0 = t0))

  } else if (white_noise == "fBm") {
    # Not need for the paper ...
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocov_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_autocov_mc, file = file_name)
  rm(mean_file_name, dt_mean, autocov_risk_file_name, dt_autocov_risk, dt_autocovtilde, dt_autocov_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_autocov_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean risk function ----
## design 1
# estim_autocov_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
#
# estim_autocov_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
#
# ## Scenario 3:
# ### FAR
# estim_autocov_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
# estim_autocov_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
#
# estim_autocov_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
# estim_autocov_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", lag = 1, s0 = s0, t0 = t0, n_cores = 30)

# Scenario 4 : ----
estim_autocov_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
estim_autocov_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)
estim_autocov_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1, n_cores = 30)

estim_autocov_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
estim_autocov_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
estim_autocov_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)
estim_autocov_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d4", lag = 1, s0 = s0, t0 = t0, n_cores = 30)

