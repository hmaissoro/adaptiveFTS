library(data.table)

# Simulation global parameters----
sig <- 0.25
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

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# Simulation function ----
simulate_data_fun <- function(mc_i, Ni, lambdai, t0, sig = 0.25,
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
                             hurst_fun = hurst,
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

# Simulate all MC process

simulate_data <- function(Nmc = mc, Ni = 400, lambdai = 300, t0, sig = 0.25,
                          process = "FAR",
                          process_ker = get_real_data_far_kenel,
                          process_mean = get_real_data_mean,
                          white_noise = "mfBm",
                          hurst = Hlogistic, Hvec = Hvec, design = "d1"){
  dt_res <- data.table::rbindlist(
    parallel::mclapply(seq_len(Nmc), function(mc_i, Ni, lambdai, t0, sig, process,
                                              process_ker, process_mean, white_noise, hurst, Hvec){
      dt_ <- simulate_data_fun(mc_i = mc_i, Ni = Ni, lambdai = lambdai, t0 = t0, sig = sig,
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
  print(paste0("Done : dt_mc_", process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", design,".RDS at ", Sys.time()))
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

### FAR process ----
## mfBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

## fBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

### FMA process ----
## mfBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

## fBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d1")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d1, white_noise = "mfBm",
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
ker_d1 <- function(s,t, operator_norm = 0.5){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

### FAR process ----
## mfBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d2, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d2")

## fBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d1,
              process_mean = mean_d2, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d2")

### FMA process ----
## mfBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d2, white_noise = "mfBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d2")

## fBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d1,
              process_mean = mean_d2, white_noise = "fBm",
              hurst = Hlogistic, Hvec = Hvec, design = "d2")

## Simulation - design 3 ----
##
## Logistic constant hurst function
# Hlogistic_d3 <- function(t){
#   hurst_logistic(t, h_left = 0.25, h_right = 0.3,
#                  change_point_position = 0.5, slope = 50)
# }

Hlogistic_d3 <- function(t){
  hurst_logistic(t, h_left = 0.3, h_right = 0.4,
                 change_point_position = 0.5, slope = 50)
}

### FAR process ----
## mfBm
ker_d3 <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.2)

simulate_data(Nmc = 100, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")

simulate_data(Nmc = 100, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")

simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")
parallel::mckill()
gc() ; gc() ; gc() ; gc()

simulate_data(Nmc = 100, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")

## d3_bis
Hlogistic_d3_bis <- function(t){
  0.25 + 0 * t
}
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_bis")

simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_bis")

simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_bis")
gc() ; gc() ; gc() ; gc()

simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3,
              process_mean = get_real_data_mean, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_bis")

### FMA process ----
## mfBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d3,
              process_mean = mean_d2, white_noise = "mfBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")

## fBm
simulate_data(Nmc = 100, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FMA", process_ker = ker_d3,
              process_mean = mean_d2, white_noise = "fBm",
              hurst = Hlogistic_d3, Hvec = Hvec, design = "d3")


