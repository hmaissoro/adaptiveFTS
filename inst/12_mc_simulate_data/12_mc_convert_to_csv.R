## Convert to csv
convert_to_csv_fun <- function(N = 150, lambda = 40, design = "d1"){
  file_list <- list.files(paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda, "/"))
  parallel::mclapply(file_list, function(ff){
    data_file <- paste0("./inst/12_mc_simulate_data/FAR/data/", design, "_slice/N", N, "lambda", lambda, "/", ff)
    dt <- readRDS(data_file)
    dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
    save_file <- paste0("../../matlab/rubin_mean_autocov_method/", design, "_slice/N", N, "lambda", lambda, "/",
                        gsub(".RDS", ".csv", ff))
    write.csv(x = dt, file = save_file)
  }, mc.cores = 35)
}

convert_to_csv_fun(N = 150, lambda = 40, design = "d1")
convert_to_csv_fun(N = 1000, lambda = 40, design = "d1")
convert_to_csv_fun(N = 400, lambda = 300, design = "d1")
convert_to_csv_fun(N = 1000, lambda = 1000, design = "d1")

convert_to_csv_fun(N = 150, lambda = 40, design = "d3")
convert_to_csv_fun(N = 1000, lambda = 40, design = "d3")
convert_to_csv_fun(N = 400, lambda = 300, design = "d3")
convert_to_csv_fun(N = 1000, lambda = 1000, design = "d3")

convert_to_csv_fun(N = 150, lambda = 40, design = "d5")
convert_to_csv_fun(N = 1000, lambda = 40, design = "d5")
convert_to_csv_fun(N = 400, lambda = 300, design = "d5")
convert_to_csv_fun(N = 1000, lambda = 1000, design = "d5")
