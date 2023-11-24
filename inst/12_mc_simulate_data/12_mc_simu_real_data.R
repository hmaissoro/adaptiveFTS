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
get_real_data_mean_all <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    240.851203112216, 0.511305802240387, 0.0686061390530684, 0.404870711185904,
    0.163861146404104, 0.114790691976136, 0.104308592101677, 0.142469296029503,
    0.130499379756486, -0.0361260524404399, -0.0678188776098126, 0.00380999048511739,
    -0.0503325208862268, -0.0224294395747781, -0.0577702125832923, -0.0392160482687889,
    0.0287914974650662, -0.0088367407619057, 0.0342607237900621, 0.0273805020027171,
    0.0316142323237759, -0.0041035647844982, -0.00674473387407733, -0.013464230018548,
    -0.0534193098888712, 0, -0.0314081783398105, -0.00635579389377677, -0.0109138999402073,
    -0.00258921655373403, 0.0210582431037979, 0.00270039632612643, 0.0220059327400361,
    0.00498016602040829, 0.00612809662966225, -0.000468228505337897, -0.015169843922895,
    -0.0128459609384212, -0.0106467952350991, -0.0127090790046628, -0.00165413277187792,
    0.00462608335044012, 0.00161552211498331, 0.00387471684545869, 0.0168292422759318,
    0.00724332519473612, -0.00489353583983273, -0.0127769082071871, -0.0178778867020446,
    -0.0194110273298402, -0.0188371035035502, 0.422453270228039, 1.00539572630755,
    -0.24482290499743, -0.261068571761477, 0.00307387023311774, -0.116983058709441,
    -0.232838350406328, -0.0329785119290798, -0.0816897661769816, -0.0647624678180593,
    -0.0220324199101271, 0.0885506433535071, 0.0146671099049825, -0.037714370452663,
    0.0436172347334241, 0.0274110323188496, 0.00692447566082562, -0.0230790920539398, 0,
    -0.0055879729086357, -0.024539268711286, -0.00612357817533665, 0.0330311138598858,
    -0.0349369137221206, 0.00811058206224394, 0.0086712004621395, 0.00485477680861141,
    0.000146907988786585, 0.017293445884018, 0.00678226738589246, -0.0237728236554015,
    -0.00555310365060365, 0.00790951023830742, -0.0134462843070566, -0.0194592780694928,
    -0.0110838300364163, 0.00437116487075195, 0, 0.0165232688198793, 0.000983634631283174,
    -0.00550366281650836, 0.00681777170946617, 0, -0.00151544500544163, -0.0129957746034038,
    0.0164014779142208, 0, -0.0210844330910173, 0.00363142721126929, -0.00200062166486549
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

get_real_data_far_kenel_all <- function(s = 0.2, t = 0.3, operator_norm = 0.5){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    0.898834000399914, -0.201020928946389, -0.236306422293472, 0.157836273520224,
    1.01946351212488, -0.732328012035491, -0.148542720086944, -0.146138621170237,
    1.2498942752218, -0.262845805101749, 1.97752361499692, -0.0109391342372169,
    0.262264781394272, 0.0607664153340225, 0, 0, 0.0533916320226869, 0.0433328320283003,
    0, 0, 0.00761337862710542, 0, -0.0125748460527178, 0.0234006978007677, 0.208242319966059,
    0.103124708111194, 0.00174591039927232, 0, 0.00465195260850987, 0.0753961965397498, 0, 0,
    0, 0, 0.0375740200989248, 0.112959396051784, 0.238095512534355, 0.00607729321428306, 0, 0,
    0.143011279530699, 0.0318162577318275, -0.0244062233040078, 0, 0, 0.0453436880346667,
    0.0569139676816576, -8.87068855687694e-05, 0.225992788863933, 0.0145639760309315,
    -0.0637891015214431, -0.0330855529173304, 0.0565341880432872, 0.0128435034029214, 0,
    0.020467626032159, 0.0163418561957687, 0, 0, 0.0640521639428075, 0.201895467411645,
    -0.0189203783109461, 0.000103927869561494, 0.0163412381713805, 0, 0.0683672469636619,
    -0.0048884113709868, 0.0262713119122406, -0.0389549091547215, 0.154478808389718,
    -0.0784238035098326, 0, 0.258766115097621, 0.264930155593137, -0.151898071985805,
    -0.245727860836431, 0.060028283121522, -0.0024269422889436, 0.0395698627288913,
    0.0553485446969482, 0.027356336342131, -0.128196819234174, 0.0316074694721634,
    0.0946611667970468, 0.528595114162558, -0.0283715450353531, -0.322289325862917, 0,
    0.0333992601855177, -0.003502762254535, 0.0123171004670106, -0.0184507793031618,
    0.0350729045700525, 0, -0.096999393039024, -0.0470379375757071, 0.159039956211493,
    0.0868772924626463, -0.147372681109599, -0.00829179510305117, -0.0543763300163787,
    -0.0456116471106727, 0, 0, -0.0158253729218884, 0.00185249667365958, -0.133180555340607,
    0.0426019180545379, 0.212961857294148, -0.0477126915276283, -0.0262173009754391,
    -0.0254349902338133, -0.0277831179697749, 0, -0.0800833395190534, 0.0468971951127422,
    0.0391750128621857, 0.0186352741893125, -0.0848502585970805, -0.0647570634965357,
    0.21055148878229
  )
  # Transform to (K, L) matrix
  basis_coef_mat <- t(matrix(data = basis_coef, ncol = 11))

  ker_values <- mapply(function(s,t, coef_mat){
    # \eta(s)
    etas <- c(1, sqrt(2) * cos(2 * pi * 1:5 * s), sqrt(2) * sin(2 * pi * 1:5 * s))

    # \theta(t)
    thetat <- c(1, sqrt(2) * cos(2 * pi * 1:5 * t), sqrt(2) * sin(2 * pi * 1:5 * t))

    # Basis function
    ker_val <- matrix(etas, nrow = 1) %*% coef_mat %*% matrix(thetat, ncol = 1)
    return(c(ker_val))
  }, s = s, t = t, MoreArgs = list(coef_mat = basis_coef_mat))

  # Normalize values using operator norm
  op_norm <- 1.009763
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef) ; gc() ; gc()

  return(ker_values)
}

Hlogistic_d3_all <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

ker_d3_bis_all <- function(s,t) get_real_data_far_kenel_all(s = s, t = t, operator_norm = 0.2)

#
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_all,
              process_mean = get_real_data_mean_all, white_noise = "mfBm",
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_all,
              process_mean = get_real_data_mean_all, white_noise = "mfBm",
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_all,
              process_mean = get_real_data_mean_all, white_noise = "mfBm",
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis_all,
              process_mean = get_real_data_mean_all, white_noise = "mfBm",
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")


estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all")
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all")

Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

gridExtra::grid.arrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Ht"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Lt"),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_all", param = "Lt"),
  nrow = 2, ncol = 2
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

