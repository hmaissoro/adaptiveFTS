library(data.table)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Mean function estimation function ----

## Estimate mean function only for "plus_mean"
estim_mean_bw_rp <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  # Load raw data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    # Estimate local regularity by mc
    dt_mean_rp_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai){
      # Extract data
      dt_random_mc <- data[ttag == "trandom"][id_mc == mc_i]
      K <- 20
      b0 <- 8 * (N * lambda) ** (- 1)
      bK <- 8 * (N * lambda) ** (- 1 / 3)
      a <- exp((log(bK) - log(b0)) / K)
      bw_grid <- b0 * a ** (seq_len(K))
      rm(K, b0, bK, a) ; gc()

      # Estimate the mean function
      dt_mean_bw_rp <- estimate_mean_bw_rp(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        Kfold = 10, bw_grid = bw_grid, smooth_ker = epanechnikov)
      dt_mean_bw_rp[, h[which.min(cv_error)]]
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai,
                                       "optbw_rp" = dt_mean_bw_rp[, h[which.min(cv_error)]])
      rm(dt_mean_bw_rp, dt_random_mc) ; gc()
      return(dt_res)
    }, mc.cores = 75, data = dt, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    if (N == 1000 & lambda == 1000)
      index_mc <- index_mc[-50]
    dt_optbw <- dt_mean_risk[, .("optbw" = h[which.min(mean_risk)]),
                             by = c("id_mc", "t", "Htrue")]
    # Estimate local regularity by mc
    dt_mean_rp_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai){
      # Extract and sort data
      data_mc <- data[id_mc == mc_i]
      Hvec <- data_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, data, dt_optbw, Ni, lambdai){
        # Extract data
        dt_random_mc_Hi <- data[ttag == "trandom"][id_mc == mc_i & Htrue == Hi]
        optbw <- dt_optbw[id_mc == mc_i & Htrue == Hi][order(t), optbw]
        t0 <- dt_optbw[id_mc == mc_i & Htrue == Hi][order(t), t]

        # Estimate the mean function
        dt_mean <- estimate_mean(
          data = dt_random_mc_Hi, idcol = "id_curve",
          tcol = "tobs", ycol = "X",
          t = t0, optbw = optbw)

        # Return and clean
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "Htrue" = Hi, dt_mean[, .(t, optbw, PN, muhat)])
        dt_res <- data.table::merge.data.table(
          x = dt_res,
          y = unique(data[ttag == "tcommon", .("t" = tobs, "Htrue" = Htrue, "mutrue" = process_mean)]),
          by = c("t", "Htrue"))
        rm(dt_mean, dt_random_mc_Hi, optbw) ; gc()
        return(dt_res)
      }, data = data_mc, dt_optbw = dt_optbw, Ni = Ni, lambdai = lambdai))

      return(dt_by_Hvec)
    }, mc.cores = 34, data = dt, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_optbw_rp_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_rp_mc, file = file_name)
  rm(data_file_name, file_name, dt_mean_rp_mc, dt, index_mc) ; gc() ; gc()

  return(paste0("Done : dt_mean_optbw_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate rubin and panaretos bw ----
# Scenario 1 :
estim_mean_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1")

# Scenario 3 :
estim_mean_bw_rp(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
estim_mean_bw_rp(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3")
estim_mean_bw_rp(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
estim_mean_bw_rp(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
estim_mean_bw_rp(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3")

# Estimate mean RP function ----

estim_mean_rp_fun <- function(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0){
  # Load raw data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_optbw_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_optbw_rp_",
                               process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  index_mc <- dt[, unique(id_mc)]

  dt_optbw <- readRDS(dt_optbw_file_name)

  if (white_noise == "mfBm") {
    # Estimate local regularity by mc
    dt_mean_rp_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai, t0){
      # Extract data
      dt_random_mc <- data[ttag == "trandom"][id_mc == mc_i]
      optbw_rp <- dt_optbw[id_mc == mc_i, optbw_rp]
      # Estimate the mean function
      dt_mean_rp <- estimate_mean_rp(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        t = t0, h = optbw_rp, smooth_ker = epanechnikov)
      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_mean_rp)
      rm(dt_mean_rp, dt_random_mc) ; gc()
      return(dt_res)
    }, mc.cores = 50, data = dt, dt_optbw = dt_optbw, Ni = N, lambdai = lambda, t0 = t0))

  } else if (white_noise == "fBm") {
    # Comming ...
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_rp_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_rp_mc, file = file_name)
  rm(data_file_name, file_name, dt_mean_rp_mc, dt_optbw_file_name, dt, dt_optbw) ; gc() ; gc()

  return(paste0("Done : dt_mean_rp_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

estim_mean_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_rp_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)

# Scenario 3 :
estim_mean_rp_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)
estim_mean_rp_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)



