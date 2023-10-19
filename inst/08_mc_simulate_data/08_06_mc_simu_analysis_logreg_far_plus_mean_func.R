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
  dt <- readRDS(paste0("./inst/08_mc_simulate_data/data/dt_mc_far_N=", N, "_lambda=", lambda, "_centered_", design, ".RDS"))
  dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]
  id_mc_vec <- dt[, unique(id_mc)]

  # Estimate local regularity by mc
  dt_reg_mc <- data.table::rbindlist(parallel::mclapply(id_mc_vec, function(mc_i, Ni, lambdai){
    # Extract and sort data
    dt_mc <- dt[id_mc == mc_i]
    dt_mc <- dt_mc[order(id_curve, tobs)]
    bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

    # Estimate the local regularity
    ## Delta
    lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
    delta <- exp(- log(lambdahat) ** 0.35)
    # delta <- (1 / lambdahat) ** 0.5

    ## For exponential Delta
    dt_locreg <- estimate_locreg(
      data = dt_mc, idcol = "id_curve",
      tcol = "tobs", ycol = "X",
      t = t0, Delta = delta, h = bw,
      smooth_ker = epanechnikov, center = center)

    ## Return obtained result
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg)
    return(dt_res)
  }, mc.cores = 75, Ni = N, lambdai = lambda))

  ## Save
  if (center) {
    file_name <- paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_far_centered_N=", N, "_lambda=", lambda, "_", design, ".RDS")
  } else {
    file_name <- paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_far_centered_N=", N, "_lambda=", lambda, "_not_centered", design, ".RDS")
  }
  saveRDS(object = dt_reg_mc, file = file_name)
  return(dt_reg_mc)
}

## Estimate local regularity (Centered) ----
### d1
dt_N400_lambda300_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1")
# dt_N1000_lambda1000_d1 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d1")

### d2
dt_N400_lambda300_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2")
# dt_N1000_lambda1000_d2 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d2")

### d3
dt_N400_lambda300_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3")
# dt_N1000_lambda1000_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3")

## Estimate local regularity (Not centered) ----
### d1
dt_N400_lambda300_not_centered_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1", center = FALSE)
# dt_N1000_lambda1000_not_centered_d1 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d1", center = FALSE)

### d2
dt_N400_lambda300_not_centered_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2", center = FALSE)
# dt_N1000_lambda1000_not_centered_d2 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d2", center = FALSE)

### d3
dt_N400_lambda300_not_centered_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3", center = FALSE)
# dt_N1000_lambda1000_not_centered_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3", center = FALSE)

