# Simulation global parameters----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

# Estimate mean function ----
## Estimation function
estim_mean_fun <- function(N = 400, lambda = 300, design = "d1"){
  dt <- readRDS(paste0("./inst/08_mc_simulate_data/data/dt_mc_far_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
  dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]
  dt_locreg <- readRDS(paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_N=", N, "_lambda=", lambda, "_", design, ".RDS"))

  # Estimate local regularity by mc
  dt_mean_mc <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i, Ni, lambdai, dt_locreg){
    # Extract and sort data
    dt_mc <- dt[id_mc == mc_i]
    dt_mc <- dt_mc[order(id_curve, tobs)]
    bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

    # Exponential grid
    lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
    K <- 100
    # b0 <- 1 / (N * lambdahat)
    b0 <- 0.005
    bK <- 5 * (N * lambdahat) ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))

    # Estimate mean function
    Ht <- dt_locreg[id_mc == mc_i][order(t), Ht]
    Lt <- dt_locreg[id_mc == mc_i][order(t), Lt]
    dt_mean <- estimate_mean(
      data = dt_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
      t = t0, bw_grid = bw_grid,
      Ht = Ht, Lt = Lt, Delta = NULL, h = bw,
      smooth_ker = epanechnikov)

    ## Return obtained result
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_mean)
    return(dt_res)
  }, mc.cores = 75, Ni = N, lambdai = lambda, dt_locreg = dt_locreg))

  ## Save
  saveRDS(
    object = dt_mean_mc,
    file = paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_mean_N=", N, "_lambda=", lambda, "_", design, ".RDS")
  )
  return(dt_mean_mc)
}

## Estimate mean function
### d1
dt_mean_N400_lambda300_d1 <- estim_mean_fun(N = 400, lambda = 300, design = "d1")
rm(dt_mean_N400_lambda300_d1) ; gc()
dt_mean_N1000_lambda1000_d1 <- estim_mean_fun(N = 1000, lambda = 1000, design = "d1")
rm(dt_mean_N1000_lambda1000_d1) ; gc()

### d2
dt_mean_N400_lambda300_d2 <- estim_mean_fun(N = 400, lambda = 300, design = "d2")
rm(dt_mean_N400_lambda300_d2) ; gc()
dt_mean_N1000_lambda1000_d2 <- estim_mean_fun(N = 1000, lambda = 1000, design = "d2")
rm(dt_mean_N1000_lambda1000_d2) ; gc()

### d3
dt_mean_N400_lambda300_d3 <- estim_mean_fun(N = 400, lambda = 300, design = "d3")
rm(dt_mean_N400_lambda300_d3) ; gc()
dt_mean_N1000_lambda1000_d3 <- estim_mean_fun(N = 1000, lambda = 1000, design = "d3")
rm(dt_mean_N1000_lambda1000_d3) ; gc()
