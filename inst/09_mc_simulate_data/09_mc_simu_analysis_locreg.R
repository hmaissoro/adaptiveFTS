# Simulation global parameters----
## Simulation global parameters----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hvec <- c(0.4, 0.5, 0.7)

# Local regularity estimation function ----
estim_locreg_fun <- function(N = 400, lambda = 300, design = "d1", center = FALSE){
  if (center) {
    dt <- readRDS(paste0("./inst/09_mc_simulate_data/data/dt_mc_far_fBm_N=", N, "_lambda=", lambda, "_centered_", design, ".RDS"))
    dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]
  } else {
    dt <- readRDS(paste0("./inst/09_mc_simulate_data/data/dt_mc_far_fBm_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
    dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]
  }

  # Estimate local regularity by mc
  dt_reg_mc <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i, Ni, lambdai){
    # Extract and sort data
    dt_mc <- dt[id_mc == mc_i]
    dt_mc <- dt_mc[order(id_curve, tobs)]

    # Estimate the local regularity
    ## Delta
    lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
    delta <- 2 * exp(- log(lambdahat) ** 0.72)
    # delta <- (1 / lambdahat) ** 0.5

    ## For exponential Delta
    dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_mc, t0, bw, delta, center){
      bw <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
      dt_locreg <- estimate_locreg(
        data = dt_mc[Htrue == Hi], idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = delta, h = bw,
        smooth_ker = epanechnikov, center = TRUE)
      dt_locreg[, Htrue := Hi]
    }, dt_mc = dt_mc, t0 = t0, bw = bw, delta = delta, center = center))


    ## Return obtained result
    dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_by_Hvec)
    return(dt_res)
  }, mc.cores = 25, Ni = N, lambdai = lambda))

  ## Save
  if (center) {
    file_name <- paste0("./inst/09_mc_simulate_data/locreg_estimates/dt_locreg_fBm_N=", N, "_lambda=", lambda, "_centered_", design, ".RDS")
  } else {
    file_name <- paste0("./inst/09_mc_simulate_data/locreg_estimates/dt_locreg_fBm_N=", N, "_lambda=", lambda, "_", design, ".RDS")
  }
  saveRDS(object = dt_reg_mc, file = file_name)
  return(dt_reg_mc)
}

## Estimate local regularity----
### d1
dt_N400_lambda300_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1")
rm(dt_N400_lambda300_d1) ; gc()

dt_N400_lambda300_centered_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1", center = TRUE)
rm(dt_N400_lambda300_centered_d1) ; gc()

### d1_bis
dt_N400_lambda300_d1_bis <- estim_locreg_fun(N = 400, lambda = 300, design = "d1_bis")
rm(dt_N400_lambda300_d1_bis) ; gc()

dt_N400_lambda300_centered_d1_bis <- estim_locreg_fun(N = 400, lambda = 300, design = "d1_bis", center = TRUE)
rm(dt_N400_lambda300_centered_d1_bis) ; gc()

### d2
dt_N400_lambda300_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2")
rm(dt_N400_lambda300_d2) ; gc()
dt_N400_lambda300_d2_bis <- estim_locreg_fun(N = 400, lambda = 300, design = "d2_bis")
rm(dt_N400_lambda300_d2_bis) ; gc()


### d3
dt_N400_lambda300_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3")
dt_N1000_lambda1000_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3")

## Estimate local regularity (Not centered) ----
### d1
dt_N400_lambda300_not_centered_d1 <- estim_locreg_fun(N = 400, lambda = 300, design = "d1", center = FALSE)
dt_N1000_lambda1000_not_centered_d1 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d1", center = FALSE)

### d2
dt_N400_lambda300_not_centered_d2 <- estim_locreg_fun(N = 400, lambda = 300, design = "d2", center = FALSE)
dt_N1000_lambda1000_not_centered_d2 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d2", center = FALSE)

### d3
dt_N400_lambda300_not_centered_d3 <- estim_locreg_fun(N = 400, lambda = 300, design = "d3", center = FALSE)
dt_N1000_lambda1000_not_centered_d3 <- estim_locreg_fun(N = 1000, lambda = 1000, design = "d3", center = FALSE)

