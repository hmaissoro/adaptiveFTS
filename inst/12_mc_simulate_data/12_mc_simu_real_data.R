library(data.table)

t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hvec <- c(0.4, 0.5, 0.7)

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}


Hlogistic_d3 <- function(t){
  hurst_logistic(t, h_left = 0.35, h_right = 0.65,
                 change_point_position = 0.6, slope = 5)
}

get_real_data_mean_bis <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    242.466591816367, 0.20197443026174, 0.140309556534677, 0.354532601148011,
    0.106485439213739, 0.209292561333274, 0.649941621023466, 1.4463552202352,
    -0.308858146955855, -0.402541770481494, 0.145343656745755
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

get_real_data_far_kenel_bis <- function(s = 0.2, t = 0.3, operator_norm = 0.5){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    0.887265486496153, -0.158284777828367, -0.433123270896265, -0.383368407909871,
    0.145655492369033, -0.00932791858596785, 0.25405721049976, 0.0360507006945973,
    0.0389539855934984, 0, -0.0133553863644848, 0.0177582032888235, 0.189761421268642,
    0.0195864450427664, 0.0887495150023169, 0, 0.0347257788913602, 0, 0.298938773778208,
    0.360062724244617, 0.00694075505838772, 0.0383993219719295, 0.0889742879270508,
    0.108124616829882, 0.597015339786177
  )
  # Transform to (K, L) matrix
  basis_coef_mat <- t(matrix(data = basis_coef, ncol = 5))

  ker_values <- mapply(function(s,t, coef_mat){
    # \eta(s)
    etas <- c(1, sqrt(2) * cos(2 * pi * 1:2 * s), sqrt(2) * sin(2 * pi * 1:2 * s))

    # \theta(t)
    thetat <- c(1, sqrt(2) * cos(2 * pi * 1:2 * t), sqrt(2) * sin(2 * pi * 1:2 * t))

    # Basis function
    ker_val <- matrix(etas, nrow = 1) %*% coef_mat %*% matrix(thetat, ncol = 1)
    return(c(ker_val))
  }, s = s, t = t, MoreArgs = list(coef_mat = basis_coef_mat))

  # Normalize values using operator norm
  op_norm <- 0.9128311
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef)
  gc()

  return(ker_values)
}


### FAR process ----
## mfBm
ker_d3_bis <- function(s,t) get_real_data_far_kenel_bis(s = s, t = t, operator_norm = 0.2)

# simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3,
#               process_mean = get_real_data_mean, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")
#
# simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3,
#               process_mean = get_real_data_mean, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")
#
# dt <- simulate_data_fun(mc_i = 1, Ni = 400, lambdai = 300, t0 = t0, sig = 0.25,
#                         process = "FAR",
#                         process_ker = ker_d3_bis,
#                         process_mean = get_real_data_mean,
#                         white_noise = "mfBm",
#                         hurst = Hlogistic_d3, Hvec = Hvec)
# dt <- dt[ttag == "trandom"]
# dt_mc <- copy(dt)
# dt_mc[, id_curve := as.factor(id_curve)]
# library(ggplot2)
#
# ggplot(data = dt_mc, mapping = aes(x = tobs, y = X, group = id_curve)) +
#   geom_line() +
#   theme_minimal()
#
# dt_locreg <- estimate_locreg(
#   data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
#   t = c(0.2, 0.3, 0.4, 0.6, 0.7, 0.8), Delta = exp(-log(300) ** (1/3)),
#   h = unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw], smooth_ker = epanechnikov, center = TRUE)
#
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")

simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")

simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")

simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")

simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")

simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")

simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")

simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")

simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")

simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")


## Merge data simulated data
merge_data_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211"){
  # design = "d3_new_2211"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/to_delete/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  id_mc_max <- max(dt[, unique(id_mc)])

  # design = "d3_new_2211_2"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/to_delete/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_2.RDS")
  dt_2 <- readRDS(file_title)

  # Merge
  dt_2[, id_mc := id_mc + id_mc_max]
  dt_merge <- rbind(dt, dt_2)

  # Save
  saveRDS(object = dt_merge,
          file = paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_d3.RDS"))
}

merge_data_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")


# Simulate data for large (N, lambda)
simulate_data_large_N_fun <- function(mc_i, Ni, lambdai, t0, sig = 0.25,
                                      process = "FAR",
                                      process_ker = get_real_data_far_kenel,
                                      process_mean = get_real_data_mean,
                                      white_noise = "mfBm",
                                      hurst = Hlogistic, Hvec = Hvec){
  # Generate data according the choixe of the white noise
  if (white_noise == "mfBm") {

    ## Generate FAR(1) or FMA(1)
    if (process == "FAR") {
      ### Simulate mean-zero FAR(1) process
      dt_gen <- simulate_far(
        N = Ni, lambda = lambdai, tdesign = "random", Mdistribution = bounded_uniform,
        tdistribution = runif, tcommon = t0, hurst_fun = hurst, L = 4, far_kernel = process_ker,
        far_mean = process_mean, int_grid = 100L, burnin = 100L, remove_burnin = TRUE)
      ### Add mean function
      dt_gen[, process_mean := far_mean]
      dt_gen[, far_mean := NULL]

    } else if (process == "FMA") {
      ### Simulate mean-zero FMA(1) process
      dt_gen <- simulate_fma(
        N = Ni, lambda = lambdai, tdesign = "random", Mdistribution = bounded_uniform,
        tdistribution = runif, tcommon = t0, hurst_fun = hurst, L = 4, fma_kernel = process_ker,
        fma_mean = process_mean, int_grid = 100L, burnin = 100L, remove_burnin = TRUE)
      ### Add mean function
      dt_gen[, process_mean := fma_mean]
      dt_gen[, far_mean := NULL]
    }

    ### Get pre-smoothing bandwidth
    #### Define and exponential bandwidth grid
    lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
    K <- 30
    b0 <- 1 / lambdahat
    bK <- lambdahat ** (- 1 / 3)
    a <- exp((log(bK) - log(b0)) / K)
    bw_grid <- b0 * a ** (seq_len(K))
    rm(b0, bK, a, K) ; gc()

    #### Add noise and get optimal bw for each curve
    index <- dt_gen[, unique(id_curve)]
    dt <- data.table::rbindlist(parallel::mclapply(index, function(id, dtt, bw_grid){
      # Filter data
      d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]

      # Add noise
      d[, X := X + rnorm(n = .N, mean = 0, sd = sig)]

      # Get optimal bandwidth
      bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
      d[, presmooth_bw := bw]
      return(d)
    }, mc.cores = 75, dtt = dt_gen, bw_grid = bw_grid))

    ### Add fix design points
    dt_tcommon <- data.table::merge.data.table(
      x = dt_gen[ttag == "tcommon"],
      y = unique(dt[, .(id_curve, presmooth_bw)]),
      by = "id_curve")
    dt_res <- rbind(dt, dt_tcommon)
    rm(dt, dt_tcommon, index, bw_grid) ; gc() ; gc()

    # Add MC index
    dt_res[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lambdai)]
    data.table::setcolorder(
      x = dt_res,
      neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "process_mean", "X", "presmooth_bw"))

  } else if (white_noise == "fBm") {
    dt_res <- data.table::rbindlist(lapply(Hvec, function(Hi, process){
      ## Generate FAR(1) or FMA(1)
      Hfun <-  function(t) Hi + 0 * t
      if (process == "FAR") {
        ### Simulate mean-zero FAR(1) process
        dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                               tdesign = "random",
                               Mdistribution = bounded_uniform,
                               tdistribution = runif,
                               tcommon = t0,
                               hurst_fun = Hfun,
                               L = 4,
                               far_kernel = process_ker,
                               far_mean = process_mean,
                               int_grid = 100L,
                               burnin = 100L,
                               remove_burnin = TRUE)
        ### Add mean function
        dt_gen[, process_mean := far_mean]
        dt_gen[, far_mean := NULL]

      } else if (process == "FMA") {
        ### Simulate mean-zero FMA(1) process
        dt_gen <- simulate_fma(N = Ni, lambda = lambdai,
                               tdesign = "random",
                               Mdistribution = bounded_uniform,
                               tdistribution = runif,
                               tcommon = t0,
                               hurst_fun = Hfun,
                               L = 4,
                               fma_kernel = process_ker,
                               fma_mean = process_mean,
                               int_grid = 100L,
                               burnin = 100L,
                               remove_burnin = TRUE)
        ### Add mean function
        dt_gen[, process_mean := fma_mean]
        dt_gen[, far_mean := NULL]
      }

      ### Get pre-smoothing bandwidth
      #### Define an exponential bandwidth grid
      lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
      K <- 50
      b0 <- 1 / lambdahat
      bK <- lambdahat ** (- 1 / 3)
      a <- exp((log(bK) - log(b0)) / K)
      bw_grid <- b0 * a ** (seq_len(K))
      rm(b0, bK, a, K) ; gc()

      #### Add noise and get optimal bw for each curve
      index <- dt_gen[, unique(id_curve)]
      dt <- data.table::rbindlist(lapply(index, function(id, dtt, bw_grid){
        # Filter data
        d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]

        # Add noise
        d[, X := X + rnorm(n = .N, mean = 0, sd = sig)]

        # Get optimal bandwidth
        bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
        d[, presmooth_bw := bw]
        return(d)
      }, dtt = dt_gen, bw_grid = bw_grid))

      ### Add fix design points
      dt_tcommon <- data.table::merge.data.table(
        x = dt_gen[ttag == "tcommon"],
        y = unique(dt[, .(id_curve, presmooth_bw)]),
        by = "id_curve")
      dt_by_H <- rbind(dt, dt_tcommon)
      rm(dt, dt_tcommon, index, bw_grid) ; gc() ; gc()

      # Add MC index
      dt_by_H[, c("id_mc", "N", "lambda", "Htrue") := .(mc_i, Ni, lambdai, Hi)]
      data.table::setcolorder(
        x = dt_by_H,
        neworder = c("id_mc", "N", "lambda", "Htrue", "id_curve", "tobs", "ttag", "process_mean", "X", "presmooth_bw"))
      return(dt_by_H)
    }, process = process))
  }

  # return the result
  return(dt_res)
}

system.time(
  dt_mc_i <- simulate_data_large_N_fun(
    mc_i = 1, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
    process = "FAR", process_ker = ker_d3_bis,
    process_mean = get_real_data_mean_bis, white_noise = "mfBm",
    hurst = Hlogistic_d3, Hvec = Hvec)
)
# Old version ----
Hlogistic_d3_old <- function(t){
  hurst_logistic(t, h_left = 0.25, h_right = 0.45,
                 change_point_position = 0.6, slope = 5)
}
ker_d3_bis_old <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.2)

#
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_old,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_old, Hvec = Hvec, design = "d3_old_mean_old_ker")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_old,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_old, Hvec = Hvec, design = "d3_old_mean_old_ker")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_old,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_old, Hvec = Hvec, design = "d3_old_mean_old_ker")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_old,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_old, Hvec = Hvec, design = "d3_old_mean_old_ker")


estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")

estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker")
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker")
estim_locreg_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker")

Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.35, h_right = 0.65,
                 change_point_position = 0.6, slope = 5)
}
gridExtra::grid.arrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Ht"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Lt"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211", param = "Lt"),
  nrow = 2, ncol = 4
)

Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.25, h_right = 0.45,
                 change_point_position = 0.6, slope = 5)
}

gridExtra::grid.arrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Ht"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Lt"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_old_mean_old_ker", param = "Lt"),
  nrow = 2, ncol = 4
)

# If we use All real data ----
Hlogistic_d3_all <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

ker_d3_all <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.7)

#
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
#
gc() ; gc()
simulate_data(Nmc = 25, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 25, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")

merge_data_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all"){
  # design = "d3_new_2211"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  id_mc_max <- dt[, max(unique(id_mc))]

  for(i in 1:n_sup_basis){
    file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_", i, ".RDS")
    dt_2 <- readRDS(file_name)

    # Merge
    dt_2[, id_mc := id_mc + id_mc_max]
    dt <- rbind(dt, dt_2)
    id_mc_max <- dt[, max(unique(id_mc))]
  }

  # Save
  saveRDS(object = dt,
          file = paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_d3.RDS"))
}

merge_data_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", n_sup_basis = 2, design = "d3_all")




estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all")
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_all")

Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

gridExtra::grid.arrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Ht", Hfun = Hlogistic_d3_all),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Ht", Hfun = Hlogistic_d3_all),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Ht", Hfun = Hlogistic_d3_all),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Lt", Hfun = Hlogistic_d3_all),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Lt", Hfun = Hlogistic_d3_all),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Lt", Hfun = Hlogistic_d3_all),
  nrow = 2, ncol = 3
)

file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                     process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
dt <- readRDS(file_title)
dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
dt_mc <- copy(dt[id_mc == 1])
dt_mc[, id_curve := as.factor(id_curve)]
library(ggplot2)

ggplot(data = dt_mc, mapping = aes(x = tobs, y = X, group = id_curve)) +
  geom_line() +
  theme_minimal()

## Estimate localisation
lambda <- 150
K <- 50
b0 <- 1 / lambda
bK <- lambda ** (- 1 / 3)
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(b0, bK, a, K) ; gc()

###
epanechnikov_bis <- function(u){
  ifelse(abs(u) <= 1, (3/4) * (1 - u ** 2), 0)
}
###


optbw <- estimate_nw_bw(y = dt_mc[id_curve == 1, X],
                        t = dt_mc[id_curve == 1, tobs],
                        bw_grid = bw_grid,
                        smooth_ker = epanechnikov)
dt_nw <- estimate_nw(y = dt_mc[id_curve == 1, X],
                     t = dt_mc[id_curve == 1, tobs],
                     tnew = dt_mc[id_curve == 1, tobs],
                     h = optbw,
                     smooth_ker = epanechnikov)
dygraphs::dygraph(data = dt_nw[, .("t" = tnew, "Xhat" = yhat)])
dygraphs::dygraph(data = dt_mc[id_curve == 1, .(tobs, X)])
## Load data
file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                    process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
dt_locreg <- readRDS(file_name)

dt_locreg[, .("Ht" = median(Ht), "Lt" = median(Lt)), by = t]
