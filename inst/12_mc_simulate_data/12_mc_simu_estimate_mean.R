library(data.table)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Mean function estimation function ----
## Estimate mean functions risk
estim_mean_risk_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0){
  # Load data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  dt_random <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
  dt_common <- dt[ttag == "tcommon"][, .SD, .SDcols = ! "ttag"]
  index_mc <- dt_random[, unique(id_mc)]
  rm(dt) ; gc()

  # Load local regularity estimates
  locreg_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                             process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(locreg_file_name)

  if (white_noise == "mfBm") {
    # Define bandwidth grid
    K <- 20
    b0 <- 4 * (N * lambda) ** (- 0.9)
    bK <- 4 * (N * lambda) ** (- 1 / (2 * 0.6 + 1))
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(K, b0, bK, a) ; gc()

    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Extract data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]

      # Estimate pseudo-true mean function risk
      dt_risk_mutilde <- data.table::rbindlist(lapply(bw_grid, function(hi, data_common_mc, data_random_mc, bw_grid, t0){
        dt_Xhat <- data.table::rbindlist(lapply(data_random_mc[, unique(id_curve)], function(curve_index, t0, hi, data){
          Tn <- data[id_curve == curve_index, tobs]
          Yn <- data[id_curve == curve_index, X]
          Xhat <- estimate_nw(y = Yn, t = Tn, tnew = t0, h = hi, smooth_ker = epanechnikov)
          return(data.table("curve_index" = curve_index, "t" = t0, "h" = hi, "Xhat" = Xhat$yhat))
        }, data = data_random_mc, t0 = t0, hi = hi))

        # Estimate mean function and add true meanfunction
        dt_muhat <- dt_Xhat[!is.nan(Xhat), .("muhat" = mean(Xhat)), by = c("t", "h")]
        dt_muhat <- data.table::merge.data.table(
          x = dt_muhat,
          y = data_common_mc[, .("mutilde" = mean(X)), by = "tobs"],
          by.x = "t", by.y = "tobs")
        dt_muhat[, mutitle_mse := (muhat - mutilde) ** 2, by = t]

      }, data_common_mc = dt_common_mc, data_random_mc = dt_random_mc, bw_grid = bw_grid, t0 = t0))

      # Estimate the risk of the mean function
      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg[id_mc == mc_i & order(t), Ht]
      Lt <- dt_locreg[id_mc == mc_i & order(t), Lt]

      dt_risk_muhat <- estimate_mean_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = 'tobs', ycol = "X",
        t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
      )

      ## Add the MSE of mutilde
      dt_risk <- data.table::merge.data.table(
        x = dt_risk_muhat, y = dt_risk_mutilde, by = c("t", "h")
      )
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_risk)
      rm(dt_risk_muhat, dt_risk_mutilde, dt_risk) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, t0){
      # Extract and sort data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]
      Hvec <- dt_random_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_common_mc, dt_random_mc, dt_locreg, Ni, lambdai, t0){
        # Define bandwidth grid
        K <- 20
        b0 <- 4 * (N * lambda) ** (- 0.9)
        bK <- 4 * (N * lambda) ** (- 1 / (2 * (Hi + 0.05) + 1))
        a <- exp((log(bK) - log(b0)) / K)
        bw_grid <- b0 * a ** (seq_len(K))
        rm(K, b0, bK, a) ; gc()
        # Estimate pseudo-true mean function risk
        dt_risk_mutilde <- data.table::rbindlist(lapply(bw_grid, function(hi, data_common_mc, data_random_mc, t0){
          dt_Xhat <- data.table::rbindlist(lapply(data_random_mc[, unique(id_curve)], function(curve_index, t0, hi, data){
            Tn <- data[id_curve == curve_index, tobs]
            Yn <- data[id_curve == curve_index, X]
            Xhat <- estimate_nw(y = Yn, t = Tn, tnew = t0, h = hi, smooth_ker = epanechnikov)
            return(data.table("curve_index" = curve_index, "t" = t0, "h" = hi, "Xhat" = Xhat$yhat))
          }, data = data_random_mc, t0 = t0, hi = hi))

          # Estimate mean function and add true meanfunction
          dt_muhat <- dt_Xhat[!is.nan(Xhat), .("muhat" = mean(Xhat)), by = c("t", "h")]
          dt_muhat <- data.table::merge.data.table(
            x = dt_muhat,
            y = data_common_mc[, .("mutilde" = mean(X)), by = "tobs"],
            by.x = "t", by.y = "tobs")
          dt_muhat[, mutitle_mse := (muhat - mutilde) ** 2, by = t]

        }, data_common_mc = dt_common_mc[Htrue == Hi], data_random_mc = dt_random_mc[Htrue == Hi], t0 = t0))

        # Estimate the risk of the mean function
        bw <- unique(dt_random_mc[Htrue == Hi][, .(id_curve, presmooth_bw)])[, presmooth_bw]
        Ht <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Ht]
        Lt <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Lt]

        dt_risk_muhat <- estimate_mean_risk(
          data = dt_random_mc[Htrue == Hi], idcol = "id_curve", tcol = 'tobs', ycol = "X",
          t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
        )

        ## Add the MSE of mutilde and return
        dt_risk <- data.table::merge.data.table(
          x = dt_risk_muhat, y = dt_risk_mutilde, by = c("t", "h")
        )
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "Htrue" = Hi, dt_risk)
        rm(dt_risk_muhat, dt_risk_mutilde, dt_risk) ; gc()
        return(dt_res)
      }, dt_common_mc = dt_common_mc, dt_random_mc = dt_random_mc, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, t0 = t0))

      return(dt_by_Hvec)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, t0 = t0))
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_risk_mc, file = file_name)
  rm(data_file_name, locreg_file_name, file_name, dt_risk_mc) ; gc() ; gc()

  return(paste0("Done : dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function only for "plus_mean"
estim_mean_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  # Load mean risk data
  mean_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                                process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean_risk <- readRDS(mean_risk_file_name)
  # Load raw data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    dt_optbw <- dt_mean_risk[!is.nan(mean_risk), .("optbw" = h[which.min(mean_risk)]), by = c("id_mc", "t")]

    # Estimate local regularity by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai){
      # Extract data
      dt_random_mc <- data[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
      t0 <- dt_optbw[id_mc == mc_i][order(t), t]

      # Estimate the mean function
      dt_mean <- estimate_mean(
        data = dt_random_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, optbw = optbw)

      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_mean[, .(t, optbw, PN, muhat)])
      dt_res <- data.table::merge.data.table(
        x = dt_res,
        y = unique(data[ttag == "tcommon", .("t" = tobs, "mutrue" = process_mean)]),
        by = "t")
      rm(optbw, dt_mean) ; gc()
      return(dt_res)
    }, mc.cores = 50, data = dt, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    if (N == 1000 & lambda == 1000)
      index_mc <- index_mc[-50]
    dt_optbw <- dt_mean_risk[, .("optbw" = h[which.min(mean_risk)]),
                             by = c("id_mc", "t", "Htrue")]
    # Estimate local regularity by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai){
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
    }, mc.cores = 50, data = dt, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(data_file_name, mean_risk_file_name, file_name, dt_mean_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_mean_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean risk function ----
## design 1 ----
### FAR ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1")
estim_mean_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_fun(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1")
estim_mean_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1")
estim_mean_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1")
estim_mean_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1")

# estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
# estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1")
### FMA ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1")
estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1")

## design 2 ----
### FAR ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2")
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2")

### FMA ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2")
estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2")

## design 3 ----
### FAR ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
estim_mean_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3")

### FMA ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", t0 = t0)

estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3")
estim_mean_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3")
