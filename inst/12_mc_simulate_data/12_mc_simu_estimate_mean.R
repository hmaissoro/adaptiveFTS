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

  # Define bandwidth grid
  K <- 100
  b0 <- 2 / lambda
  bK <- lambda ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()

  if (white_noise == "mfBm") {

    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Extract data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]

      # Estimate pseudo-true mean function risk
      dt_risk_mutilde <- data.table::rbindlist(lapply(bw_grid, function(hi, data, bw_grid, t0){
        dt_Xhat <- data.table::rbindlist(lapply(dt_common_mc[, unique(id_curve)], function(curve_index, t0, hi, data){
          Tn <- data[id_curve == curve_index, tobs]
          Yn <- data[id_curve == curve_index, X_plus_mean]
          Xhat <- estimate_nw(y = Yn, t = Tn, tnew = t0, h = hi, smooth_ker = epanechnikov)
          return(data.table("curve_index" = curve_index, "t" = t0, "h" = hi, "Xhat" = Xhat$yhat))
        }, data = data, t0 = t0, hi = hi))

        # Estimate mean function and add true meanfunction
        dt_muhat <- dt_Xhat[!is.nan(Xhat), .("h" = unique(h), "mutilde" = mean(Xhat)), by = t]
        dt_muhat <- data.table::merge.data.table(
          x = dt_muhat,
          y = data[, .("mutrue" = unique(process_mean)), by = "tobs"],
          by.x = "t", by.y = "tobs")
        dt_muhat[, mutitle_mse := (mutilde - mutrue) ** 2, by = t]

      }, data = dt_common_mc, bw_grid = bw_grid, t0 = t0))

      # Estimate the risk of the mean function
      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw_plus_mean)])[, presmooth_bw_plus_mean]
      Ht <- dt_locreg[id_mc == mc_i & order(t), Ht]
      Lt <- dt_locreg[id_mc == mc_i & order(t), Lt]
      Ht_plus_mean <- dt_locreg[id_mc == mc_i & order(t), Ht_plus_mean]
      Lt_plus_mean <- dt_locreg[id_mc == mc_i & order(t), Lt_plus_mean]

      dt_risk_muhat <- estimate_mean_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = 'tobs', ycol = "X_plus_mean",
        t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
      )

      dt_risk_muhat_plus_mean <- estimate_mean_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = 'tobs', ycol = "X_plus_mean",
        t = t0, bw_grid = bw_grid, Ht = Ht_plus_mean, Lt = Lt_plus_mean, Delta = NULL, h = bw
      )
      data.table::setnames(x = dt_risk_muhat_plus_mean,
                           old = c("Ht", "Lt", "locreg_bw", "bias_term", "varriance_term", "dependence_term", "mean_risk"),
                           new = c("Ht_plus_mean", "Lt_plus_mean", "locreg_bw_plus_mean", "bias_term_plus_mean",
                                   "varriance_term_plus_mean", "dependence_term_plus_mean", "mean_risk_plus_mean"))

      ## Merge the mean risk estimates and return the obtained result
      dt_risk_muhat_res <- data.table::merge.data.table(
        x = dt_risk_muhat,
        y = dt_risk_muhat_plus_mean,
        by = c("t", "h")
      )
      rm(dt_risk_muhat, dt_risk_muhat_plus_mean) ; gc()

      ## Add the MSE of mutilde
      dt_risk <- data.table::merge.data.table(
        x = dt_risk_muhat_res, y = dt_risk_mutilde, by = c("t", "h")
      )
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_risk)
      rm(dt_risk_muhat_res, dt_risk_mutilde, dt_risk) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Extract and sort data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]
      Hvec <- dt_random_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_common_mc, dt_random_mc, dt_locreg, Ni, lambdai, bw_grid, t0){
        # Estimate pseudo-true mean function risk
        dt_risk_mutilde <- data.table::rbindlist(lapply(bw_grid, function(hi, data, bw_grid, t0){
          dt_Xhat <- data.table::rbindlist(lapply(dt_common_mc[, unique(id_curve)], function(curve_index, t0, hi, data){
            Tn <- data[id_curve == curve_index, tobs]
            Yn <- data[id_curve == curve_index, X_plus_mean]
            Xhat <- estimate_nw(y = Yn, t = Tn, tnew = t0, h = hi, smooth_ker = epanechnikov)
            return(data.table("curve_index" = curve_index, "t" = t0, "h" = hi, "Xhat" = Xhat$yhat))
          }, data = data, t0 = t0, hi = hi))

          # Estimate mean function and add true meanfunction
          dt_muhat <- dt_Xhat[!is.nan(Xhat), .("h" = unique(h), "mutilde" = mean(Xhat)), by = t]
          dt_muhat <- data.table::merge.data.table(
            x = dt_muhat,
            y = data[, .("mutrue" = unique(process_mean)), by = "tobs"],
            by.x = "t", by.y = "tobs")
          dt_muhat[, mutitle_mse := (mutilde - mutrue) ** 2, by = t]

        }, data = dt_common_mc[Htrue == Hi], bw_grid = bw_grid, t0 = t0))

        # Estimate the risk of the mean function
        bw <- unique(dt_random_mc[Htrue == Hi][, .(id_curve, presmooth_bw_plus_mean)])[, presmooth_bw_plus_mean]
        Ht <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Ht]
        Lt <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Lt]
        Ht_plus_mean <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Ht_plus_mean]
        Lt_plus_mean <- dt_locreg[Htrue == Hi][id_mc == mc_i & order(t), Lt_plus_mean]

        dt_risk_muhat <- estimate_mean_risk(
          data = dt_random_mc[Htrue == Hi], idcol = "id_curve", tcol = 'tobs', ycol = "X_plus_mean",
          t = t0, bw_grid = bw_grid, Ht = Ht, Lt = Lt, Delta = NULL, h = bw
        )

        dt_risk_muhat_plus_mean <- estimate_mean_risk(
          data = dt_random_mc[Htrue == Hi], idcol = "id_curve", tcol = 'tobs', ycol = "X_plus_mean",
          t = t0, bw_grid = bw_grid, Ht = Ht_plus_mean, Lt = Lt_plus_mean, Delta = NULL, h = bw
        )
        data.table::setnames(x = dt_risk_muhat_plus_mean,
                             old = c("Ht", "Lt", "locreg_bw", "bias_term", "varriance_term", "dependence_term", "mean_risk"),
                             new = c("Ht_plus_mean", "Lt_plus_mean", "locreg_bw_plus_mean", "bias_term_plus_mean",
                                     "varriance_term_plus_mean", "dependence_term_plus_mean", "mean_risk_plus_mean"))
        ## Merge the mean risk estimates and return the obtained result
        dt_risk_muhat_res <- data.table::merge.data.table(
          x = dt_risk_muhat,
          y = dt_risk_muhat_plus_mean,
          by = c("t", "h")
        )
        rm(dt_risk_muhat, dt_risk_muhat_plus_mean) ; gc()

        ## Add the MSE of mutilde and return
        dt_risk <- data.table::merge.data.table(
          x = dt_risk_muhat_res, y = dt_risk_mutilde, by = c("t", "h")
        )
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "Htrue" = Hi, dt_risk)
        rm(dt_risk_muhat_res, dt_risk_mutilde, dt_risk) ; gc()
        return(dt_res)
      }, dt_common_mc = dt_common_mc, dt_random_mc = dt_random_mc, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))

      return(dt_by_Hvec)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_risk_mc, file = file_name)
  rm(data_file_name, locreg_file_name, file_name, dt_risk_mc) ; gc() ; gc()

  return(paste0("Done : dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}

## Estimate mean function
estim_mean_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0){
  # comming ...
}

# Estimate local regularity ----
## FAR ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)

## FMA ----
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", t0 = t0)
estim_mean_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", t0 = t0)
