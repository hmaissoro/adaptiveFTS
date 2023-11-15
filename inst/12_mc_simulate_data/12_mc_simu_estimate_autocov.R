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
  b0 <- 4 * max((N * lambda) ** (- 0.9), (N * (lambda ** 2)) ** (- 0.9))
  bK <- 4 * max((N * lambda) ** (- 1 / 3), (N * (lambda ** 2)) ** (- 1 / 3))
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))
  rm(K, b0, bK, a) ; gc()

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
          data = data_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
          s = s0, t = t0, cross_lag = lag, lag = NULL, h = hi, smooth_ker = epanechnikov)
        return(dt_gammahat[, .("s" = s, "t" = t, "h" = hi, "gammahat" = EXsXt_cross_lag)])
      }, data_random_mc = dt_random_mc, bw_grid = bw_grid, s0 = s0, t0 = t0, lag = lag))

      ## Estimate \widetilde\gamma(s,t)
      dt_gammatilde_data <- data.table::rbindlist(lapply(dt_common_mc[, unique(id_curve)], function(curve_index, data, s0, t0){
        dt_true_X <- data.table::merge.data.table(
          x = data.table::data.table("s" = s0, "t" = t0),
          y = data[id_curve == curve_index, .("t" = tobs, "Xt" = X)],
          by = "t"
        )
        dt_true_X <- data.table::merge.data.table(
          x = dt_true_X,
          y = data[id_curve == curve_index, .("s" = tobs, "Xs" = X)],
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
        y = dt_locreg[id_mc == mc_i, .("t" = t, Ht, Lt)],
        by = "t"
      )
      dt_locreg_mc <- data.table::merge.data.table(
        x = dt_locreg_mc,
        y = dt_locreg[id_mc == mc_i, .("s" = t, "Hs" = Ht, "Ls" = Lt)],
        by = "s"
      )

      bw <- unique(dt_random_mc[, .(id_curve, presmooth_bw)])[, presmooth_bw]
      Ht <- dt_locreg_mc[order(s, t), Ht]
      Lt <- dt_locreg_mc[order(s, t), Lt]
      Hs <- dt_locreg_mc[order(s, t), Hs]
      Ls <- dt_locreg_mc[order(s, t), Ls]

      dt_risk_gammahat <- estimate_autocov_risk(
        data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
        s = s0, t = t0, lag = 1, bw_grid = bw_grid,
        smooth_ker = epanechnikov, Hs = Hs, Ls = Ls,
        Ht = Ht, Lt = Lt, Delta = NULL, h = bw)

      ## Add the MSE of mutilde
      dt_risk <- data.table::merge.data.table(x = dt_risk_gammahat, y = dt_gamma_mse, by = c("s", "t", "h"))
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_risk)
      rm(dt_risk_gammahat, dt_gamma_mse, dt_risk) ; gc()
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
        # Estimate the MSE ----
        dt_gammahat <- data.table::rbindlist(lapply(bw_grid, function(hi, data_random_mc, bw_grid, s0, t0, lag){
          dt_gammahat <- estimate_empirical_XsXt_autocov(
            data = data_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
            s = s0, t = t0, cross_lag = lag, lag = NULL, h = hi, smooth_ker = epanechnikov)
          return(dt_gammahat[, .("s" = s, "t" = t, "h" = hi, "gammahat" = EXsXt_cross_lag)])
        }, data_random_mc = dt_random_mc[Htrue == Hi], bw_grid = bw_grid, s0 = s0, t0 = t0, lag = lag))

        ## Estimate \widetilde\gamma(s,t)
        dt_gammatilde_data <- data.table::rbindlist(lapply(dt_common_mc[Htrue == Hi][, unique(id_curve)], function(curve_index, data, s0, t0){
          dt_true_X <- data.table::merge.data.table(
            x = data.table::data.table("s" = s0, "t" = t0),
            y = data[id_curve == curve_index, .("t" = tobs, "Xt" = X)],
            by = "t"
          )
          dt_true_X <- data.table::merge.data.table(
            x = dt_true_X,
            y = data[id_curve == curve_index, .("s" = tobs, "Xs" = X)],
            by = "s"
          )
          return(data.table::data.table("id_curve" = curve_index, dt_true_X))
        }, data = dt_common_mc[Htrue == Hi], s0 = s0, t0 = t0))

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
        dt_locreg_mc <- dt_locreg[id_mc == mc_i & Htrue == Hi]
        dt_locreg_mc <- data.table::merge.data.table(
          x = data.table::data.table("s" = s0, "t" = t0),
          y = dt_locreg[id_mc == mc_i & Htrue == Hi, .("t" = t, Ht, Lt)],
          by = "t"
        )
        dt_locreg_mc <- data.table::merge.data.table(
          x = dt_locreg_mc,
          y = dt_locreg[id_mc == mc_i & Htrue == Hi, .("s" = t, "Hs" = Ht, "Ls" = Lt)],
          by = "s"
        )

        bw <- unique(dt_random_mc[Htrue == Hi, .(id_curve, presmooth_bw)])[, presmooth_bw]
        Ht <- dt_locreg_mc[order(s, t), Ht]
        Lt <- dt_locreg_mc[order(s, t), Lt]
        Hs <- dt_locreg_mc[order(s, t), Hs]
        Ls <- dt_locreg_mc[order(s, t), Ls]

        dt_risk_gammahat <- estimate_autocov_risk(
          data = dt_random_mc[Htrue == Hi], idcol = "id_curve", tcol = "tobs", ycol = "X",
          s = s0, t = t0, lag = 1, bw_grid = bw_grid,
          smooth_ker = epanechnikov, Hs = Hs, Ls = Ls,
          Ht = Ht, Lt = Lt, Delta = NULL, h = bw)

        ## Add the MSE of mutilde
        dt_risk <- data.table::merge.data.table(x = dt_risk_gammahat, y = dt_gamma_mse, by = c("s", "t", "h"))
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "Htrue" = Hi, dt_risk)
        rm(dt_risk_gammahat, dt_gamma_mse, dt_risk) ; gc()
        return(dt_res)
      }, dt_common_mc = dt_common_mc, dt_random_mc = dt_random_mc, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))

      return(dt_by_Hvec)
    }, mc.cores = 75, dt_random = dt_random, dt_common = dt_common, dt_locreg = dt_locreg, Ni = N, lambdai = lambda, bw_grid = bw_grid, t0 = t0))
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_auto_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_risk_mc, file = file_name)
  rm(data_file_name, locreg_file_name, file_name, dt_risk_mc) ; gc() ; gc()

  return(paste0("Done : dt_autocov_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

## Estimate mean function
estim_autocov_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1){
    ## Load mean risk data
    autocov_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_auto_risk_",
                                  process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
    dt_autocov_risk <- readRDS(autocov_risk_file_name)

    ## \widetilde{\mu}_N(t)
    mean_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                                  process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
    dt_mean_risk <- readRDS(mean_risk_file_name)

    ## Load raw data
    data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                             process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
    dt <- readRDS(data_file_name)
    index_mc <- dt[, unique(id_mc)]

    if (white_noise == "mfBm") {
      # Estimate optimal bandwidth and autocovariance proxy
      dt_optbw <- dt_autocov_risk[, .("optbw" = h[which.min(autocov_risk)], gammatilde), by = c("id_mc", "s", "t")]
      dt_optbw <- data.table::merge.data.table(
        x = dt_optbw,
        y = unique(dt_mean_risk[, .(id_mc, "s" = t, "mutilde_s" = mutilde)]),
        by = c("id_mc", "s")
      )
      dt_optbw <- data.table::merge.data.table(
        x = dt_optbw,
        y = unique(dt_mean_risk[, .(id_mc, "t" = t, "mutilde_t" = mutilde)]),
        by = c("id_mc", "t")
      )
      dt_optbw[, autocovtilde := gammatilde - mutilde_s * mutilde_t, by = c("s", "t")]

      # Estimate autocovariance by mc
      dt_autocov_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, data, dt_optbw, Ni, lambdai){
        # Extract data
        dt_random_mc <- data[ttag == "trandom"][id_mc == mc_i]
        optbw <- dt_optbw[id_mc == mc_i][order(s, t), optbw]
        s0 <- dt_optbw[id_mc == mc_i][order(s, t), s]
        t0 <- dt_optbw[id_mc == mc_i][order(s, t), t]

        # Estimate the mean function
        dt_autocov <- estimate_autocov(
          data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X",
          s = s0, t = t0, lag = 1, optbw = optbw, bw_grid = seq(0.005, 0.15, len = 45),
          Hs = NULL, Ls = NULL, Ht = NULL, Lt = NULL,
          Delta = NULL, h = NULL, center = TRUE,
          mean_estimates_s = NULL, mean_estimates_t = NULL,
          smooth_ker = epanechnikov)

        # Add estimate of the true gamma
        dt_autocov <- data.table::merge.data.table(
          x = dt_autocov,
          y = unique(dt_optbw[id_mc == mc_i, .(s, t, mutilde_s, mutilde_t, autocovtilde)]),
          by = c("s", "t")
        )
        # Return and clean
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, dt_autocov)
        rm(dt_optbw, optbw) ; gc()
        return(dt_res)
      }, mc.cores = 50, data = dt, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))

    } else if (white_noise == "fBm") {
      # Coming ...
    }

    ## Save
    file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocov_estimates_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
    saveRDS(object = dt_mean_mc, file = file_name)
    rm(data_file_name, autocov_risk_file_name, file_name, dt_mean_mc, dt_optbw) ; gc()

    return(paste0("Done : dt_autocov_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate mean risk function ----
## design 1 ----
### FAR ----
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)

estim_autocov_risk_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", s0 = s0, t0 = t0, lag = 1)

# estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0 = t0)
# estim_autocov_risk_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", t0 = t0)

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

