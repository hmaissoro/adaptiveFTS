library(data.table)

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

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

# zero-mean function
zero_mean_func <- function(t) 0 * t

# Simulation function ----
simate_data_fun <- function(mc_i, Ni, lambdai, t0, sig = 0.25,
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
      dt_gen <- simulate_far(N = Ni, lambda = lambdai,
                             tdesign = "random",
                             Mdistribution = bounded_uniform,
                             tdistribution = runif,
                             tcommon = t0,
                             hurst_fun = hurst,
                             L = 4,
                             far_kernel = process_ker,
                             far_mean = zero_mean_func,
                             int_grid = 100L,
                             burnin = 100L,
                             remove_burnin = TRUE)
      dt_gen[, far_mean := NULL]

    } else if (process == "FMA") {
      ### Simulate mean-zero FMA(1) process
      dt_gen <- simulate_fma(N = Ni, lambda = lambdai,
                             tdesign = "random",
                             Mdistribution = bounded_uniform,
                             tdistribution = runif,
                             tcommon = t0,
                             hurst_fun = hurst,
                             L = 4,
                             fma_kernel = process_ker,
                             fma_mean = zero_mean_func,
                             int_grid = 100L,
                             burnin = 100L,
                             remove_burnin = TRUE)
      dt_gen[, fma_mean := NULL]
    }

    ### Add mean function
    dt_gen[, process_mean := process_mean(tobs)]
    dt_gen[, X_plus_mean := X + process_mean]

    ### Get pre-smoothing bandwidth
    #### Define and exponential bandwidth grid
    lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
    K <- 100
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
      d[, X_plus_mean := X_plus_mean + rnorm(n = .N, mean = 0, sd = sig)]

      # Get optimal bandwidth
      bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
      bw_plus_mean <- estimate_nw_bw(y = d[, X_plus_mean], t = d[, tobs], bw_grid = bw_grid)
      d[, presmooth_bw := bw]
      d[, presmooth_bw_plus_mean := bw_plus_mean]
      return(d)
    }, dtt = dt_gen, bw_grid = bw_grid))

    ### Add fix design points
    dt_tcommon <- data.table::merge.data.table(
      x = dt_gen[ttag == "tcommon"],
      y = unique(dt[, .(id_curve, presmooth_bw, presmooth_bw_plus_mean)]),
      by = "id_curve")
    dt_res <- rbind(dt, dt_tcommon)
    rm(dt, dt_tcommon, index, bw_grid) ; gc() ; gc()

    # Add MC index
    dt_res[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lambdai)]
    data.table::setcolorder(
      x = dt_res,
      neworder = c("id_mc", "N", "lambda", "id_curve", "tobs", "ttag", "process_mean",
                   "X", "X_plus_mean", "presmooth_bw", "presmooth_bw_plus_mean"))

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
                               far_mean = zero_mean_func,
                               int_grid = 100L,
                               burnin = 100L,
                               remove_burnin = TRUE)
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
                               fma_mean = zero_mean_func,
                               int_grid = 100L,
                               burnin = 100L,
                               remove_burnin = TRUE)
        dt_gen[, fma_mean := NULL]
      }

      ### Add mean function
      dt_gen[, process_mean := process_mean(tobs)]
      dt_gen[, X_plus_mean := X + process_mean]

      ### Get pre-smoothing bandwidth
      #### Define an exponential bandwidth grid
      lambdahat <- mean(dt_gen[ttag == "trandom", .N, by = id_curve][, N])
      K <- 100
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
        d[, X_plus_mean := X_plus_mean + rnorm(n = .N, mean = 0, sd = sig)]

        # Get optimal bandwidth
        bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
        bw_plus_mean <- estimate_nw_bw(y = d[, X_plus_mean], t = d[, tobs], bw_grid = bw_grid)
        d[, presmooth_bw := bw]
        d[, presmooth_bw_plus_mean := bw_plus_mean]
        return(d)
      }, dtt = dt_gen, bw_grid = bw_grid))

      ### Add fix design points
      dt_tcommon <- data.table::merge.data.table(
        x = dt_gen[ttag == "tcommon"],
        y = unique(dt[, .(id_curve, presmooth_bw, presmooth_bw_plus_mean)]),
        by = "id_curve")
      dt_by_H <- rbind(dt, dt_tcommon)
      rm(dt, dt_tcommon, index, bw_grid) ; gc() ; gc()

      # Add MC index
      dt_by_H[, c("id_mc", "N", "lambda", "Htrue") := .(mc_i, Ni, lambdai, Hi)]
      data.table::setcolorder(
        x = dt_by_H,
        neworder = c("id_mc", "N", "lambda", "Htrue", "id_curve", "tobs", "ttag",
                     "process_mean", "X", "X_plus_mean", "presmooth_bw", "presmooth_bw_plus_mean"))
      return(dt_by_H)
    }, process = process))
  }

  # return the result
  return(dt_res)
}

# Simulate all MC process

simate_data <- function(Nmc = mc, Ni = 400, lambdai = 300, t0, sig = 0.25,
                        process = "FAR",
                        process_ker = get_real_data_far_kenel,
                        process_mean = get_real_data_mean,
                        white_noise = "mfBm",
                        hurst = Hlogistic, Hvec = Hvec, design = "d1"){
  dt_res <- data.table::rbindlist(
    parallel::mclapply(seq_len(Nmc), function(mc_i, Ni, lambdai, t0, sig, process,
                                              process_ker, process_mean, white_noise, hurst, Hvec){
      dt_ <- simate_data_fun(mc_i = mc_i, Ni = Ni, lambdai = lambdai, t0 = t0, sig = sig,
                             process = process, process_ker = process_ker,
                             process_mean = process_mean, white_noise = white_noise,
                             hurst = hurst, Hvec = Hvec)
      return(dt_)
    }, Ni = Ni, lambdai = lambdai, t0 = t0, sig = sig,
    process = process, process_ker = process_ker, process_mean = process_mean,
    white_noise = white_noise, hurst = hurst, Hvec = Hvec, mc.cores = 75))

  ### Local Regularity
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                       process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", design,".RDS")

  saveRDS(object = dt_res, file = file_title)
  rm(dt_res, file_title) ; gc() ; gc()
}

# Data generation ----
## Simulation - design 1 ----

## Mean function
mean_d1 <- function(t) 4 * sin(3 * pi * t / 2)

## Autoregressive kernel
ker_d1 <- function(s,t, operator_norm = 0.7){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

### FAR process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

simate_data(Nmc = 100, Ni = 400, lambdai = 50, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

simate_data(Nmc = 100, Ni = 400, lambdai = 50, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

### FMA process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

simate_data(Nmc = 100, Ni = 400, lambdai = 50, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")

simate_data(Nmc = 100, Ni = 400, lambdai = 50, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d1, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d1")


## Simulation - design 2 ----
## Mean function
mean_d2 <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  eta <- splines2::nsp(x = t, df = 12, intercept = TRUE)

  # Basis coeffients
  basis_coef <- matrix(
    data = c(731.389571901671, -0.52006982827689, 246.098166747946,
             243.059070562643, 239.632261351878, 242.246665268361,
             241.476817297158, 246.281122481508, 241.410779068207,
             238.415397856336, -5.9613322334532, 734.400332002359),
    ncol = 1)

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(eta, basis_coef) ; gc()

  return(muhat[, 1])
}

## Autoregressive kernel
ker_d1 <- function(s,t, operator_norm = 0.7){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

### FAR process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d2, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d2")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d1,
            process_mean = mean_d2, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d2")

### FMA process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d2, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d2")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d1,
            process_mean = mean_d2, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d2")

## Simulation - design 3 ----
## Mean function
mean_d2 <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  eta <- splines2::nsp(x = t, df = 12, intercept = TRUE)

  # Basis coeffients
  basis_coef <- matrix(
    data = c(731.389571901671, -0.52006982827689, 246.098166747946,
             243.059070562643, 239.632261351878, 242.246665268361,
             241.476817297158, 246.281122481508, 241.410779068207,
             238.415397856336, -5.9613322334532, 734.400332002359),
    ncol = 1)

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(eta, basis_coef) ; gc()

  return(muhat[, 1])
}

## Autoregressive kernel
ker_d3 <- function(s = 0.2, t = 0.3, operator_norm = 0.7){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    -152.339167113767, 109.891801325509, 8.01218483618224, 40.7773902997773,
    9.29538508447265, -130.869348227434, 114.349788382373, -86.6452258918943,
    66.6620550406495, 2.85516854777547, 5.4726154735248, 8.16996388562305,
    -56.2974864050636, 82.3309005212671, -104.821148135886, 70.806466322378,
    3.20308498330693, 6.34208585976186, 12.8240898172477, -68.4244432976903,
    57.5035523478083, -79.009963744684, 55.5058015942419, 7.42288493639332,
    20.1537864363113, 2.35376692860042, -61.7830117509254, 77.2298006322408,
    -80.679220971239, 59.4135944681524, 4.1358055271817, 18.1150446807183,
    5.8204767292337, -60.5952551453644, 60.9830028643384, -137.595821369872,
    83.2753892033171, -7.64204301659673, -4.06678434607622, 11.8936539033712,
    -66.1910800695624, 85.8780261168699, -142.247945431835, 121.819129805756,
    9.45741883242566, 39.8753096159282, 6.85395396362438, -132.236594114549,
    128.212492255768)

  # Transform to (K, L) matrix
  basis_coef_mat <- matrix(data = basis_coef, ncol = 1)

  # \eta(s)
  etas <- splines2::nsp(x = s, df = 5 + 1 + 1, intercept = TRUE)

  # \theta(t)
  thetat <- splines2::nsp(x = t, df = 5 + 1 + 1, intercept = TRUE)

  mat_basis <- lapply(1:ncol(etas), function(id_eta_col, eta_mat, theta_mat){
    eta_mat[, id_eta_col] * theta_mat
  }, eta_mat = etas, theta_mat = thetat)
  mat_basis <- do.call(cbind, mat_basis)

  # Basis function
  ker_values <-  mat_basis %*% basis_coef_mat
  ker_values <-  c(ker_values)

  # Normalize values using operator norm
  op_norm <- 3.446037
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef, etas, thetat, mat_basis) ; gc()

  return(ker_values)
}

### FAR process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d3,
            process_mean = mean_d2, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FAR", process_ker = ker_d3,
            process_mean = mean_d2, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3")

### FMA process ----
## mfBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d3,
            process_mean = mean_d2, white_noise = "mfBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3")

## fBm
simate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
            process = "FMA", process_ker = ker_d3,
            process_mean = mean_d2, white_noise = "fBm",
            hurst = Hlogistic, Hvec = Hvec, design = "d3")


