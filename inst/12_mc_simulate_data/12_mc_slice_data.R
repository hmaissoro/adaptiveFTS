slice_data <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3"){
  # Load data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  # dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    # slice data
    parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai){
      # Extract data
      dt_mc <- dt[id_mc == mc_i]
      ## Save
      file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", Ni, "lambda", lambdai,"/dt_mc_",
                          process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      saveRDS(object = dt_mc, file = file_name)
      rm(dt_mc, file_name) ; gc()

    }, mc.cores = 35, dt = dt, Ni = N, lambdai = lambda)

  } else if (white_noise == "fBm") {
    # Coming ...
  }

  return(paste0("Done : Slice for dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Slice d3 data
slice_data(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
slice_data(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3")
slice_data(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
slice_data(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
slice_data(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3")

# Slice d1 data
slice_data(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
slice_data(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
slice_data(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
slice_data(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1")
