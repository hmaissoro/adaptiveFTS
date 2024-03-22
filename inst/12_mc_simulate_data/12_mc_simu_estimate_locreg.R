library(data.table)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# local regularity estimation ----
estim_locreg_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0, nc_cores = 75){

  if (white_noise == "mfBm") {
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)

    # Nmc <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda)))
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(id_mc_data, function(mc_i, Ni, lambdai, t0){
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)
      dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]

      # Extract and sort data
      dt <- dt[order(id_curve, tobs)]
      bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

      # Estimate the local regularity
      ## Delta
      lambdahat <- mean(dt[, .N, by = id_curve][, N])
      delta <- min(exp(- log(lambdahat) ** ( 1 / 3)), 0.2)

      ## Centered process
      start_time <- Sys.time()
      dt_locreg <- estimate_locreg(
        data = dt, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = delta, h = bw,
        smooth_ker = epanechnikov, center = TRUE)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)

      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai,
                                       "lambdahat" = lambdahat, "time_taken" = time_taken, dt_locreg)
      rm(dt_locreg) ; gc()
      return(dt_res)
    }, mc.cores = nc_cores, Ni = N, lambdai = lambda, t0 = t0))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    # dt_reg_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai, t0){
    #   # Extract and sort data
    #   dt_mc <- dt[id_mc == mc_i]
    #   dt_mc <- dt_mc[order(id_curve, tobs)]
    #   Hvec <- dt_mc[, sort(unique(Htrue))]
    #
    #   dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_mc, t0){
    #     ## Set Delta
    #     lambdahat <- mean(dt_mc[Htrue == Hi, .N, by = id_curve][, N])
    #     # delta <- exp(- log(lambdahat) ** (1 / 4))
    #     delta <- min(exp(- log(lambdahat) ** (1 / 3)), 0.2)
    #
    #     ## Extract bandwidth
    #     bw <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
    #
    #     ## Centered process
    #     dt_locreg <- estimate_locreg(
    #       data = dt_mc[Htrue == Hi], idcol = "id_curve",
    #       tcol = "tobs", ycol = "X",
    #       t = t0, Delta = delta, h = bw,
    #       smooth_ker = epanechnikov, center = TRUE)
    #
    #     ## Merge local regularity parameters estimates and return the obtained result
    #     dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg)
    #     dt_res[, Htrue := Hi]
    #
    #     ## Clean
    #     rm(dt_locreg) ; gc()
    #     return(dt_res)
    #   }, dt_mc = dt_mc, t0 = t0))
    #
    #   return(dt_by_Hvec)
    # }, mc.cores = 75, dt = dt, Ni = N, lambdai = lambda, t0 = t0))
  }
  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_reg_mc, file = file_name)
  rm(file_name, dt_reg_mc) ; gc() ; gc()

  return(paste0("Done : dt_locreg_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate local regularity ----
## design 1 ----
### FAR ----
# mfBm
estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0, nc_cores = 30)

# fBm
# estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1")
# estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1")
# estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1")
# estim_locreg_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1")

# ### FMA ----
# estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1")
# estim_locreg_fun(N = 150, lambda = 40, process = "FMA", white_noise = "mfBm", design = "d1")
# estim_locreg_fun(N = 1000, lambda = 40, process = "FMA", white_noise = "mfBm", design = "d1")
# estim_locreg_fun(N = 1000, lambda = 1000, process = "FMA", white_noise = "mfBm", design = "d1")
#
# # fBm
# estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1")
#
# ## design 2 ----
# ### FAR
# estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2")
# estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2")
#
# ### FMA
# estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2")
# estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2")

## design 3 ----
### FAR
estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, nc_cores = 120)
estim_locreg_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0, nc_cores = 30)

## design 4 ----
### FAR
estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d4", t0 = t0, nc_cores = 30)

## design 5 ----
### FAR
estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, nc_cores = 30)
estim_locreg_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d5", t0 = t0, nc_cores = 30)

