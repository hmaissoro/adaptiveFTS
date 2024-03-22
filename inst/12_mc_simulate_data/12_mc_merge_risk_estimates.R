merge_mean_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  data_path <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda)
  data_file_done <- list.files(data_path)
  dt_res <- data.table::rbindlist(lapply(data_file_done, function(ff){
    dt <- readRDS(file.path(data_path, ff))
  }))
  saveRDS(object = dt_res,
          file = paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}

merge_mean_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
merge_mean_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
merge_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
merge_mean_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1")

merge_mean_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
merge_mean_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
merge_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
merge_mean_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3")

merge_mean_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5")
merge_mean_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d5")
merge_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d5")
merge_mean_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d5")


merge_autocov_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  data_path <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/", design, "_autocov_risk/N", N, "lambda", lambda)
  data_file_done <- list.files(data_path)
  dt_res <- data.table::rbindlist(lapply(data_file_done, function(ff){
    dt <- readRDS(file.path(data_path, ff))
  }))
  saveRDS(object = dt_res,
          file = paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocov_risk_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}



merge_autocov_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
merge_autocov_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
merge_autocov_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
merge_autocov_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1")

merge_autocov_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
merge_autocov_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3")
merge_autocov_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
merge_autocov_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3")


merge_autocov_risk(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4")
merge_autocov_risk(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d4")
merge_autocov_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d4")
merge_autocov_risk(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d4")
