library(data.table)

# Simulation global parameters----
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

# (N, \lambda) \in {(400, 300), (1000, 1000), (150, 40), (1000, 40)}
# (N, \lambda) \in {(150, 40), (1000, 40)}
# (N, \lambda) = (150, 1000)

# Note that (N, \lambda) = (1000, 1000) for design = "d1" not computed

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

## Logistic constant hurst function
Hvec <- c(0.4, 0.5, 0.7)

# Simulation function ----
simate_data_common_fun <- function(mc_i, Ni, lambdai, t0, process = "FAR",
                                   process_ker = function(s,t) 9/4 * exp( - (t + 2 * s) ** 2),
                                   process_mean = function(t) 4 * sin(1.5 * pi * t),
                                   white_noise = "mfBm",
                                   hurst = Hlogistic, Hvec = Hvec){
  # Generate data according the choixe of the white noise
  if (white_noise == "mfBm") {

    ## Generate FAR(1) or FMA(1)
    if (process == "FAR") {
      ### Simulate mean-zero FAR(1) process
      dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                             tdesign = "common",
                             Mdistribution = NULL,
                             tdistribution = NULL,
                             tcommon = t0,
                             hurst_fun = hurst,
                             L = 4,
                             far_kernel = process_ker,
                             far_mean = process_mean,
                             int_grid = 100L,
                             burnin = 500L,
                             remove_burnin = TRUE)
      ### Add mean function
      dt_gen[, process_mean := far_mean]
      dt_gen[, far_mean := NULL]

    } else if (process == "FMA") {
      ### Simulate mean-zero FMA(1) process
      dt_gen <- simulate_fma(N = Ni, lambda = lambdai,
                             tdesign = "common",
                             Mdistribution = NULL,
                             tdistribution = NULL,
                             tcommon = t0,
                             hurst_fun = hurst,
                             L = 4,
                             fma_kernel = process_ker,
                             fma_mean = process_mean,
                             int_grid = 100L,
                             burnin = 500L,
                             remove_burnin = TRUE)
      ### Add mean function
      dt_gen[, process_mean := fma_mean]
      dt_gen[, fma_mean := NULL]
    }

    # Add MC index
    dt_gen[, c("id_mc", "N") := .(mc_i, Ni)]
    data.table::setcolorder(x = dt_gen, neworder = c("id_mc", "N", "id_curve", "tobs", "ttag", "process_mean", "X"))

  } else if (white_noise == "fBm") {
    dt_gen <- data.table::rbindlist(lapply(Hvec, function(Hi, process){
      ## Generate FAR(1) or FMA(1)
      Hfun <-  function(t) Hi + 0 * t
      if (process == "FAR") {
        ### Simulate mean-zero FAR(1) process
        dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                               tdesign = "common",
                               Mdistribution = NULL,
                               tdistribution = NULL,
                               tcommon = t0,
                               hurst_fun = Hfun,
                               L = 4,
                               far_kernel = process_ker,
                               far_mean = process_mean,
                               int_grid = 100L,
                               burnin = 500L,
                               remove_burnin = TRUE)
        ### Add mean function
        dt_gen[, process_mean := far_mean]
        dt_gen[, far_mean := NULL]

      } else if (process == "FMA") {
        ### Simulate mean-zero FMA(1) process
        dt_gen <- simulate_fma(N = Ni, lambda = lambdai,
                               tdesign = "common",
                               Mdistribution = NULL,
                               tdistribution = NULL,
                               tcommon = t0,
                               hurst_fun = Hfun,
                               L = 4,
                               fma_kernel = process_ker,
                               fma_mean = process_mean,
                               int_grid = 100L,
                               burnin = 500L,
                               remove_burnin = TRUE)
        ### Add mean function
        dt_gen[, process_mean := fma_mean]
        dt_gen[, fma_mean := NULL]
      }

      # Add MC index
      dt_gen[, c("id_mc", "N", "Htrue") := .(mc_i, Ni, Hi)]
      data.table::setcolorder(x = dt_gen, neworder = c("id_mc", "N", "Htrue", "id_curve", "tobs", "ttag", "process_mean", "X"))
      return(dt_gen)
    }, process = process))
  }
  # return the result
  return(dt_gen)
}

# Simulate all MC process

simate_data <- function(Nmc = mc, Ni = 5000, t0,
                        process = "FAR",
                        process_ker = get_real_data_far_kenel,
                        process_mean = get_real_data_mean,
                        white_noise = "mfBm",
                        hurst = Hlogistic, Hvec = Hvec, design = "d1", for_variance = FALSE){
  dt_res <- data.table::rbindlist(
    parallel::mclapply(seq_len(Nmc), function(mc_i, Ni, lambdai, t0, process,
                                              process_ker, process_mean, white_noise, hurst, Hvec){
      dt_ <- simate_data_common_fun(mc_i = mc_i, Ni = Ni, lambdai = 100, t0 = t0,
                                    process = process, process_ker = process_ker,
                                    process_mean = process_mean, white_noise = white_noise,
                                    hurst = hurst, Hvec = Hvec)
      return(dt_)
    }, Ni = Ni, lambdai = lambdai, t0 = t0,
    process = process, process_ker = process_ker, process_mean = process_mean,
    white_noise = white_noise, hurst = hurst, Hvec = Hvec, mc.cores = 30))

  ### Local Regularity
  if (for_variance) {
    file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_common_for_variance_",
                         process,"_", white_noise, "_", "N=", Ni, "_", design,".RDS")
  } else {
    file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_common",
                         process,"_", white_noise, "_", "N=", Ni, "_", design,".RDS")
  }

  saveRDS(object = dt_res, file = file_title)
  rm(dt_res, file_title) ; gc() ; gc()
}

# Data generation ----
## Simulation - design 1 ----

## Mean function
mean_d1 <- function(t) 4 * sin(3 * pi * t / 2)

## Autoregressive kernel
ker_d1 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

## For auto-covariance
simate_data(Nmc = 100, Ni = 5000, t0,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

## For variance
simate_data(Nmc = 30, Ni = 5000, t0 = sort(c(t0, seq(0.1, 0.9, len = 20))),
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1", for_variance = TRUE)

## For autocovariance
zero_mean <- function(t) 0
simate_data(Nmc = 1, Ni = 5000, t0 = seq(0.1, 0.9, len = 200),
            process = "FAR", process_ker = ker_d1,
            process_mean = zero_mean, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "zero_mean_d1", for_variance = TRUE)

## Simulation - design 3 ----
## For auto-covariance
ker_d3 <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.7)
simate_data(Nmc = 100, Ni = 5000, t0,
            process = "FAR", process_ker = ker_d3,
            process_mean = get_real_data_mean, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3")

## For variance
simate_data(Nmc = 30, Ni = 5000, t0 = sort(c(t0, seq(0.1, 0.9, len = 20))),
            process = "FAR", process_ker = ker_d3,
            process_mean = get_real_data_mean, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3", for_variance = TRUE)

## Simulation - design 4 ----

## Mean function
mean_d4 <- function(t) 0

## Autoregressive kernel
ker_d1 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

simate_data(Nmc = 30, Ni = 5000, t0,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d4, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d4")


# Estimation of true autocovariance ----
t0 <- c(0.2, 0.4, 0.7, 0.8)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]

estim_true_autocov_fun <- function(N = 5000, process = "FAR",
                                   white_noise = "mfBm", design = "d1",
                                   s0 = s0, t0 = t0, lag = 1){
  # Load data
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_common",
                           process,"_", white_noise, "_", "N=", N, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    # Estimate mean function
    dt[, mutilde := mean(X), by = c("id_mc", "tobs")]

    # Estimate local regularity by mc
    dt_autocovtilde_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, N, s0, t0, lag){
      dt_mc <- dt[id_mc == mc_i]
      ## Estimate \widetilde\gamma(s,t)
      dt_gammatilde_data <- data.table::rbindlist(lapply(dt[, unique(id_curve)], function(curve_index, data, s0, t0){
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
      }, data = dt_mc, s0 = s0, t0 = t0))

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

      ## Merge add mean function dt_gammatilde
      ### Add \mu(t)
      dt_autocovtilde <- data.table::merge.data.table(
        x = dt_gammatilde,
        y = unique(dt_mc[, .("t" = tobs, "mutilde_t" = mutilde)]),
        by = "t"
      )
      ### Add \mu(s)
      dt_autocovtilde <- data.table::merge.data.table(
        x = dt_autocovtilde,
        y = unique(dt_mc[, .("s" = tobs, "mutilde_s" = mutilde)]),
        by = "s"
      )
      dt_autocovtilde[, autocovtilde := gammatilde - mutilde_s * mutilde_t]

      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N, "lag" = lag, dt_autocovtilde)
      rm(dt_gammatilde, dt_autocovtilde) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt = dt, N = N, s0 = s0, t0 = t0, lag = lag))
    dt_autocovtilde_mc <- dt_autocovtilde_mc[
      ,  .("gammatilde" = mean(gammatilde),
           "mutilde_t" = mean(mutilde_t),
           "mutilde_s" = mean(mutilde_s),
           "autocovtilde" = mean(autocovtilde)),
      by = c("s", "t", "lag")]
  } else if (white_noise == "fBm") {
    # Estimate mean function
    dt[, mutilde := mean(X), by = c("id_mc", "tobs")]

    # Estimate local regularity by mc
    dt_autocovtilde_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, N, s0, t0, lag){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      Hvec <- dt_random_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt, N, s0, t0, lag){

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
        }, data = dt_mc[Htrue == Hi], s0 = s0, t0 = t0))

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

        ### Merge add mean function dt_gammatilde
        ### Add \mu(t)
        dt_autocovtilde <- data.table::merge.data.table(
          x = dt_gammatilde,
          y = unique(dt_mc[Htrue == Hi, .("t" = tobs, "mutilde_t" = mutilde)]),
          by = "t"
        )
        ### Add \mu(s)
        dt_autocovtilde <- data.table::merge.data.table(
          x = dt_autocovtilde,
          y = unique(dt_mc[Htrue == Hi, .("s" = tobs, "mutilde_s" = mutilde)]),
          by = "s"
        )
        dt_autocovtilde[, autocovtilde := gammatilde - mutilde_s * mutilde_t]

        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = N,  "lag" = lag, "Htrue" = Htrue, dt_autocovtilde)
        rm(dt_gammatilde, dt_autocovtilde) ; gc()
        return(dt_res)
      }, dt = dt, N = N, s0 = s0, t0 = t0, lag = lag))

      return(dt_by_Hvec)
    }, mc.cores = 75, dt = dt, N = N, s0 = s0, t0 = t0, lag = lag))
    dt_autocovtilde_mc <- dt_autocovtilde_mc[
      ,  .("gammatilde" = mean(gammatilde),
           "mutilde_t" = mean(mutilde_t),
           "mutilde_s" = mean(mutilde_s),
           "autocovtilde" = mean(autocovtilde)),
      by = c("s", "t", "lag", "Htrue")]
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_autocovtilde_",
                      process,"_", white_noise, "_", "N=", N, "_", design,".RDS")
  saveRDS(object = dt_autocovtilde_mc, file = file_name)
  rm(data_file_name, file_name, dt_autocovtilde_mc) ; gc() ; gc()

  return(paste0("Done : dt_autocovtilde_", process,"_", white_noise, "_", "N=", N, "_", design,".RDS at ", Sys.time()))
}

estim_true_autocov_fun(N = 5000, process = "FAR", white_noise = "mfBm", design = "d1", s0 = s0, t0 = t0, lag = 1)
estim_true_autocov_fun(N = 5000, process = "FAR", white_noise = "mfBm", design = "d3", s0 = s0, t0 = t0, lag = 1)
estim_true_autocov_fun(N = 5000, process = "FAR", white_noise = "mfBm", design = "d4", s0 = s0, t0 = t0, lag = 1)

# Estimate zero_mean_d1 autocovariance ----
t0 <- seq(0.1, 0.9, len = 200)
dt_st <- data.table::as.data.table(expand.grid(s = t0, t = t0))
dt_st <- dt_st[order(s,t)]
s0 <- dt_st[, s]
t0 <- dt_st[, t]

# Load data
dt <- readRDS("./inst/12_mc_simulate_data/FAR/data/dt_mc_common_for_variance_FAR_mfBm_N=5000_zero_mean_d1.RDS")
N <- 5000
## Estimate \widetilde\gamma(s,t)
dt_gammatilde_data <- data.table::rbindlist(lapply(dt[, unique(id_curve)], function(curve_index, data, s0, t0){
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
}, data = dt, s0 = s0, t0 = t0))

### Take into account the lag
lag <- 1
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

saveRDS(dt_gammatilde, "./inst/12_mc_simulate_data/FAR/data/dt_true_gammatilde_FAR_mfBm_N=5000_zero_mean_d1.RDS")




