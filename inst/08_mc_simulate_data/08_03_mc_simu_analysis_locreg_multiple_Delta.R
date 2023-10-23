# Simulation global parameters----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

# Local regularity estimation function ----
estim_locreg_fun <- function(N = 400, lambda = 300, design = "d1", center = TRUE){
  dt <- readRDS(paste0("./inst/08_mc_simulate_data/data/dt_mc_far_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
  dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]

  # Chose a vector of Delta
  Delta_vec <- seq(0.01, 0.19, by = 0.005)

  # Estimate local regularity by mc
  dt_reg_mc <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i, Ni, lambdai, Delta_vec){
    # Extract and sort data
    dt_mc <- dt[id_mc == mc_i]
    dt_mc <- dt_mc[order(id_curve, tobs)]
    bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

    # Estimate the local regularity for each Delta in Delta_vec

    ## For exponential Delta
    dt_locreg <- data.table::rbindlist(lapply(Delta_vec, function(dd, t0, bw, dt_mc, center){
      dt_res <- estimate_locreg(
        data = dt_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = dd, h = bw,
        smooth_ker = epanechnikov, center = center)
      return(dt_res)
    }, t0 = t0, bw = bw, dt_mc = dt_mc, center = center))


    ## Return obtained result
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_locreg)
    return(dt_res)
  }, mc.cores = 75, Ni = N, lambdai = lambda, Delta_vec = Delta_vec))

  ## Save
  if (center) {
    file_name <- paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_Delta_Vec_N=", N, "_lambda=", lambda, "_", design, ".RDS")
  } else {
    file_name <- paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_Delta_Vec_N=", N, "_lambda=", lambda, "_not_centered_", design, ".RDS")
  }
  saveRDS(object = dt_reg_mc, file = file_name)
  return(dt_reg_mc)
}

## Estimate local regularity (Centered) ----
### d1
dt_N400_lambda300_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1")
rm(dt_N400_lambda300_d1) ; gc()
dt_N1000_lambda1000_d1 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d1")
rm(dt_N1000_lambda1000_d1) ; gc()

### d2
dt_N400_lambda300_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2")
rm(dt_N400_lambda300_d2) ; gc()
dt_N1000_lambda1000_d2 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d2")
rm(dt_N1000_lambda1000_d2) ; gc()

### d3
dt_N400_lambda300_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3")
rm(dt_N400_lambda300_d3) ; gc()
dt_N1000_lambda1000_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3")
rm(dt_N1000_lambda1000_d3) ; gc()

## Estimate local regularity (Not centered) ----
### d1
dt_N400_lambda300_not_centered_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1", center = FALSE)
rm(dt_N400_lambda300_not_centered_d1) ; gc()
dt_N1000_lambda1000_not_centered_d1 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d1", center = FALSE)
rm(dt_N1000_lambda1000_not_centered_d1) ; gc()

### d2
dt_N400_lambda300_not_centered_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2", center = FALSE)
rm(dt_N400_lambda300_not_centered_d2) ; gc()
dt_N1000_lambda1000_not_centered_d2 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d2", center = FALSE)
rm(dt_N1000_lambda1000_not_centered_d2) ; gc()

### d3
dt_N400_lambda300_not_centered_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3", center = FALSE)
rm(dt_N400_lambda300_not_centered_d3) ; gc()
dt_N1000_lambda1000_not_centered_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3", center = FALSE)
rm(dt_N1000_lambda1000_not_centered_d3) ; gc()

