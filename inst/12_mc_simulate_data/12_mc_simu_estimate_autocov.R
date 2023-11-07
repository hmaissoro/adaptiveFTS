library(data.table)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]

# Mean function estimation function ----
## Estimate mean functions risk
estim_autocov_risk_fun <- function(N = 400, lambda = 300, process = "FAR",
                                   white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1){
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
  K <- 20
  b0 <- 4 * (N * lambda) ** (- 0.9)
  bK <- 4 * (N * lambda) ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()

  # Define

  if (white_noise == "mfBm") {

    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(
      parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, bw_grid, s0, t0, lag){
      # Extract data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]

      # Estimate the MSE ----
      dt_gammahat <- data.table::rbindlist(lapply(bw_grid, function(hi, data_random_mc, bw_grid, s0, t0, lag){
        dt_gammahat <- estimate_empirical_XsXt_autocov(
          data = data_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X_plus_mean",
          s = s0, t = t0, cross_lag = lag, lag = NULL, h = hi, smooth_ker = epanechnikov)
        return(dt_gammahat[, .("s" = s, "t" = t, "h" = hi, "gammahat" = EXsXt_cross_lag)])
      }, data_random_mc = dt_random_mc, bw_grid = bw_grid, s0 = s0, t0 = t0, lag = lag))

      ## Estimate \widetilde\gamma(s,t)
      dt_gammatilde_data <- data.table::rbindlist(lapply(dt_common_mc[, unique(id_curve)], function(curve_index, data, s0, t0){
        dt_true_X <- data.table::merge.data.table(
          x = data.table::data.table("s" = s0, "t" = t0),
          y = data[id_curve == curve_index, .("t" = tobs, "Xt" = X_plus_mean)],
          by = "t"
        )
        dt_true_X <- data.table::merge.data.table(
          x = dt_true_X,
          y = data[id_curve == curve_index, .("s" = tobs, "Xs" = X_plus_mean)],
          by = "s"
        )
        return(data.table::data.table("id_curve" = curve_index, dt_true_X))
      }, data = dt_common_mc, s0 = s0, t0 = t0))

      ### Take into account the cross_lag
      ### The argument s is associated to the curves n = 1,..., N - lag
      dt_gammatilde_data_s <- dt_gammatilde_data[, list(id_curve, s, t, Xs)]
      dt_gammatilde_data_s <- dt_gammatilde_data_s[id_curve %in% 1:(N - lag)]
      dt_id_lag_s <- data.table::data.table(
        "id_curve" = sort(unique(dt_gammatilde_data_s[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_gammatilde_data_s <- data.table::merge.data.table(
        x = dt_id_lag_s,
        y = dt_gammatilde_data_s,
        by = "id_curve")

      ### The argument t is associated to the curves n = 1 + lag,..., N
      dt_gammatilde_data_t <- dt_gammatilde_data[, list(id_curve, s, t, Xt)]
      dt_gammatilde_data_t <- dt_gammatilde_data_t[id_curve %in% (1 + lag):N]
      dt_id_lag_t <- data.table::data.table(
        "id_curve" = sort(unique(dt_gammatilde_data_t[, id_curve])),
        "id_lag" = paste0(1:(N - lag), "_", (1 + lag):N))
      dt_gammatilde_data_t <- data.table::merge.data.table(
        x = dt_id_lag_t,
        y = dt_gammatilde_data_t,
        by = "id_curve")

      ### Merge and clean
      dt_gammatilde_merge <- data.table::merge.data.table(
        x = dt_gammatilde_data_s,
        y = dt_gammatilde_data_t,
        by = c("id_lag", "s", "t"))
      dt_gammatilde_merge <- dt_gammatilde_merge[order(id_curve.x)]
      rm(dt_id_lag_s, dt_id_lag_t, dt_gammatilde_data_s, dt_gammatilde_data_t, dt_gammatilde_data) ; gc()

      ### The gamma function with true X's
      dt_gammatilde <- dt_gammatilde_merge[!(is.nan(Xs) | is.nan(Xt)), .("gammatilde" = mean(Xs * Xt)), by = c("s", "t")]
      rm(dt_gammatilde_merge) ; gc()

      ## Merge dt_gammahat and dt_gammatilde
      dt_gamma_mse <- data.table::merge.data.table(x = dt_gammahat, y = dt_gammatilde, by = c("s", "t"))
      dt_gamma_mse[, gammatilde_mse := (gammahat - gammatilde) ** 2]
      rm(dt_common_mc, dt_gammahat, dt_gammatilde) ; gc()

      # Estimate the risk of the gamma function ----
      ## Extract local regularity parameters
      dt_locreg_mc <- dt_locreg[id_mc == mc_i]
      dt_locreg_mc <- data.table::merge.data.table(
        x = data.table::data.table("s" = s0, "t" = t0),
        y = dt_locreg[id_mc == mc_i, .("t" = t, Ht, Lt, Ht_plus_mean, Lt_plus_mean)],
        by = "t"
      )
      dt_locreg_mc <- data.table::merge.data.table(
        x = dt_locreg_mc,
        y = dt_locreg[id_mc == mc_i, .("s" = t, "Hs" = Ht, "Ls" = Lt, "Hs_plus_mean" = Ht_plus_mean, "Ls_plus_mean" = Lt_plus_mean)],
        by = "s"
      )

      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw_plus_mean)])[, presmooth_bw_plus_mean]
      Ht <- dt_locreg_mc[order(s, t), Ht]
      Lt <- dt_locreg_mc[order(s, t), Lt]
      Hs <- dt_locreg_mc[order(s, t), Hs]
      Ls <- dt_locreg_mc[order(s, t), Ls]
      Ht_plus_mean <- dt_locreg_mc[order(s, t), Ht_plus_mean]
      Lt_plus_mean <- dt_locreg_mc[order(s, t), Lt_plus_mean]
      Hs_plus_mean <- dt_locreg_mc[order(s, t), Hs_plus_mean]
      Ls_plus_mean <- dt_locreg_mc[order(s, t), Ls_plus_mean]

      dt_risk_gammahat <- estimate_autocov_risk(
        data, idcol = "id_curve", tcol = "tobs", ycol = "X_plus_mean",
        s = s0, t = t0, lag = 1, bw_grid =bw_grid,
        smooth_ker = epanechnikov, Hs = Hs, Ls = Ls,
        Ht = Ht, Lt = Lt, Delta = NULL, h = bw)

      dt_risk_gammahat_plus_mean <- estimate_autocov_risk(
        data, idcol = "id_curve", tcol = "tobs", ycol = "X_plus_mean",
        s = s0, t = t0, lag = 1, bw_grid =bw_grid,
        smooth_ker = epanechnikov, Hs = Hs_plus_mean, Ls = Ls_plus_mean,
        Ht = Ht_plus_mean, Lt = Lt_plus_mean, Delta = NULL, h = bw)
      # data.table::setnames(x = dt_risk_muhat_plus_mean,
      #                      old = c("Ht", "Lt", "locreg_bw", "bias_term", "varriance_term", "dependence_term", "mean_risk"),
      #                      new = c("Ht_plus_mean", "Lt_plus_mean", "locreg_bw_plus_mean", "bias_term_plus_mean",
      #                              "varriance_term_plus_mean", "dependence_term_plus_mean", "mean_risk_plus_mean"))
      #
      # ## Merge the mean risk estimates and return the obtained result
      # dt_risk_muhat_res <- data.table::merge.data.table(
      #   x = dt_risk_muhat,
      #   y = dt_risk_muhat_plus_mean,
      #   by = c("t", "h")
      # )
      # rm(dt_risk_muhat, dt_risk_muhat_plus_mean) ; gc()
      #
      # ## Add the MSE of mutilde
      # dt_risk <- data.table::merge.data.table(
      #   x = dt_risk_muhat_res, y = dt_risk_mutilde, by = c("t", "h")
      # )
      # dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_risk)
      # rm(dt_risk_muhat_res, dt_risk_mutilde, dt_risk) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg,
    Ni = N, lambdai = lambda, bw_grid = bw_grid, s0 = s0, t0 = t0, lag = lag))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    dt_risk_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_random, dt_common, dt_locreg, Ni, lambdai, bw_grid, t0){
      # Extract and sort data
      dt_common_mc <- dt_common[id_mc == mc_i]
      dt_random_mc <- dt_random[id_mc == mc_i]
      Hvec <- dt_random_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_common_mc, dt_random_mc, dt_locreg, Ni, lambdai, bw_grid, t0){
        # Estimate pseudo-true mean function risk
        dt_risk_mutilde <- data.table::rbindlist(lapply(bw_grid, function(hi, data_common_mc, data_random_mc, bw_grid, t0){
          dt_Xhat <- data.table::rbindlist(lapply(data_random_mc[, unique(id_curve)], function(curve_index, t0, hi, data){
            Tn <- data[id_curve == curve_index, tobs]
            Yn <- data[id_curve == curve_index, X_plus_mean]
            Xhat <- estimate_nw(y = Yn, t = Tn, tnew = t0, h = hi, smooth_ker = epanechnikov)
            return(data.table("curve_index" = curve_index, "t" = t0, "h" = hi, "Xhat" = Xhat$yhat))
          }, data = data_random_mc, t0 = t0, hi = hi))

          # Estimate mean function and add true meanfunction
          dt_muhat <- dt_Xhat[!is.nan(Xhat), .("muhat" = mean(Xhat)), by = c("t", "h")]
          dt_muhat <- data.table::merge.data.table(
            x = dt_muhat,
            y = data_common_mc[, .("mutilde" = mean(X_plus_mean)), by = "tobs"],
            by.x = "t", by.y = "tobs")
          dt_muhat[, mutitle_mse := (muhat - mutilde) ** 2, by = t]

        }, data_common_mc = dt_common_mc[Htrue == Hi], data_random_mc = dt_random_mc[Htrue == Hi], bw_grid = bw_grid, t0 = t0))

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

  return(paste0("Done : dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function
estim_autocov_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0){
  # coming ...
}

# Estimate mean risk function ----
## design 1 ----
### FAR ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)

### FMA ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", t0 = t0)

## design 2 ----
### FAR ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", t0 = t0)

### FMA ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", t0 = t0)

## design 3 ----
### FAR ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", t0 = t0)

### FMA ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", t0 = t0)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", t0 = t0)

