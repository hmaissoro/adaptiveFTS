# Mont2-Carlo simulation for the paper : data generation
#
# Objectif : estimation of \widetilde H and its comparaison to \widehat H and the true H.
#

library(data.table)

# Simulation global parameters----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# Simulation functions ----
## Mean functions
zero_mean_func <- function(t){
  0 * t
}
far_mean_d1 <- function(t){
  res <- 4 * sin(3 * pi * t / 2)
  return(res)
}

## Autoregressive kernel
far_ker_d1 <- function(s,t, operator_norm = 0.7){
  # Note that : \kappa_c * k = operator_norm
  k <- sqrt(pi) / 2 * (
    pnorm(q = 2, mean = 0, sd = sqrt(1/2)) - pnorm(q = 0, mean = 0, sd = sqrt(1/2))
  )
  kappa_c <- operator_norm / k
  res <- kappa_c * exp(- (s - 2 * t) ** 2)
}

## Simulation function ----
## If center = TRUE, the mean function is added to a mean-zero FAR(1) process
sim_estim_locreg_proxy <- function(mc_i, Ni, lbda, t0,
                                   kernel_far = far_ker_d1,
                                   mean_far = far_mean_d1,
                                   hurst = Hlogistic, center = FALSE){
  # Chose a vector of Delta
  Delta_vec <- seq(0.01, 0.19, by = 0.005)

  # Calculte t1, t_2, and t_3
  time_grid <- expand.grid("t0" = t0, "Delta" = Delta_vec)
  dt_time <- data.table::as.data.table(time_grid)
  dt_time[, c("t1", "t2", "t3") := .(t0 - Delta /2, t0, t0 + Delta /2)]

  dt_time <- data.table::melt(
    data = dt_time, id.vars = c("t0", "Delta"),
    measure.var = c("t1", "t2", "t3"),
    variable.name ="time_tag", value.name = "time"
  )
  dt_time[, id_time := paste0(time_tag, "_", t0, "_", Delta)]
  dt_time <- dt_time[order(time)]

  # Generate FTS
  ## If center = TRUE, the mean function is added to a mean-zero FAR(1) process
  if (center) {
    dt_far <- simulate_far(N = Ni, lambda = lbda,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = dt_time[, time],
                           hurst_fun = hurst,
                           L = 4,
                           far_kernel = kernel_far,
                           far_mean = zero_mean_func,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
  } else {
    dt_far <- simulate_far(N = Ni, lambda = lbda,
                           tdesign = "random",
                           Mdistribution = bounded_uniform,
                           tdistribution = runif,
                           tcommon = dt_time[, time],
                           hurst_fun = hurst,
                           L = 4,
                           far_kernel = kernel_far,
                           far_mean = zero_mean_func,
                           int_grid = 100L,
                           burnin = 100L,
                           remove_burnin = TRUE)
    dt_far[, X := X + mean_far(tobs)]
  }
  # Estimate local regularity parameters ----
  # Get pre-smoothing bandwidth
  ## Define and exponential bandwidth grid
  lambdahat <- mean(dt_far[ttag == "trandom", .N, by = id_curve][, N])
  K <- 100
  b0 <- 1 / lambdahat
  bK <- lambdahat ** (- 1 / 3)
  a <- exp((log(bK) - log(b0)) / K)
  bw_grid <- b0 * a ** (seq_len(K))

  ## Determine bw for each curve
  index <- dt_far[, unique(id_curve)]
  dt <- data.table::rbindlist(lapply(index, function(id, dtt, bw_grid){
    # Filter data
    d <- dtt[id_curve == id & ttag == "trandom"][order(tobs)]

    # Add noise
    d[, X := X + rnorm(n = .N, mean = 0, sd = sig)]

    # Get be
    bw <- estimate_nw_bw(y = d[, X], t = d[, tobs], bw_grid = bw_grid)
    d[, presmooth_bw := bw]
    return(d)
  }, dtt = dt_far, bw_grid = bw_grid))

  bw <- unique(dt[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

  # Estimate the local regularity for each Delta in Delta_vec
  ## For exponential Delta
  dt_locreg <- data.table::rbindlist(lapply(Delta_vec, function(dd, t0, bw, dt){
    dt_res <- estimate_locreg(
      data = dt, idcol = "id_curve",
      tcol = "tobs", ycol = "X",
      t = t0, Delta = dd, h = bw,
      center = FALSE, smooth_ker = epanechnikov)
    return(dt_res)
  }, t0 = t0, bw = bw, dt = dt))
  data.table::setnames(x = dt_locreg, old = c("Ht", "Lt"), new = c("hatH", "hatL"))

  ## clean
  rm(dt, bw) ; gc()

  # Estimate proxies ----
  ## Extract data and Add dt_time
  dt_far_tcommon <- dt_far[ttag == "tcommon"]
  dt_far_tcommon[, id_time := unique(dt_time[, id_time]), by = "id_curve"]
  dt_far_tcommon <- data.table::merge.data.table(x = dt_far_tcommon, y = dt_time, by = "id_time")

  ## Reshape merged data
  dt_far_tcommon_dcast <- data.table::dcast(
    data = dt_far_tcommon,
    formula = id_curve + Delta + t0 ~ time_tag,
    value.var = "X")
  data.table::setnames(
    x = dt_far_tcommon_dcast,
    old = c("t1", "t2", "t3"),
    new = c("xt1", "xt2", "xt3"))

  ## Clean
  rm(dt_far_tcommon, dt_time, dt_far) ; gc()

  ## Estimate theta
  dt_locreg_proxies <- dt_far_tcommon_dcast[
    between(xt1, quantile(xt1, 0.025), quantile(xt1, 0.975)) &
      between(xt1, quantile(xt1, 0.025), quantile(xt1, 0.975)) &
      between(xt1, quantile(xt1, 0.025), quantile(xt1, 0.975)),
    .("theta_t1_t2" = mean((xt1 - xt2) ** 2),
      "theta_t1_t3" = mean((xt1 - xt3) ** 2),
      "theta_t2_t3" = mean((xt2 - xt3) ** 2)),
    by = c("t0", "Delta")
  ]
  ## Estimate local regularity parameters proxies
  dt_locreg_proxies[, tildeH := (log(theta_t1_t3) - log(theta_t2_t3)) / (2 * log(2))]
  dt_locreg_proxies[, tildeL := theta_t2_t3 / ((Delta / 2) ** (2 * tildeH))]
  data.table::setnames(x = dt_locreg_proxies, old = "t0", new = "t")

  # Merge estimates and proxies and Add MC index ----
  ## Merge
  dt_locreg_res <- data.table::merge.data.table(
    x = dt_locreg,
    y = dt_locreg_proxies,
    by = c("t", "Delta")
  )
  ## Add true local regularity parameters
  dt_locreg_res[, Htrue := hurst(t)]
  dt_locreg_res[, Ltrue := 4]

  ## Add MC index
  dt_locreg_res[, c("id_mc", "N", "lambda") := .(mc_i, Ni, lbda)]
  data.table::setcolorder(
    x = dt_locreg_res,
    neworder = c("id_mc", "N", "lambda", "t", "Delta",
                 "theta_t1_t2", "theta_t1_t3", "theta_t2_t3",
                 "Htrue", "Ltrue", "tildeH", "tildeL",
                 "Nused",  "locreg_bw","hatH", "hatL"))
  return(dt_locreg)
}

# Data generation ----
## Not centered FAR - (N, lambda) = (400, 300)
dt_mc_locreg_proxy <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_estim_locreg_proxy(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0,
                                kernel_far = far_ker_d1, mean_far = far_mean_d1,
                                hurst = Hlogistic, center = FALSE)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_locreg_proxy, file = "./inst/10_mc_simulate_data/data/dt_mc_far_proxy_N=400_lambda=300.RDS")
rm(dt_mc_locreg_proxy) ; gc()

## Centered FAR - (N, lambda) = (400, 300)
dt_mc_locreg_proxy_centered <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i){
  dt_ <- sim_estim_locreg_proxy(mc_i = mc_i, Ni = 400, lbda = 300, t0 = t0,
                                kernel_far = far_ker_d1, mean_far = far_mean_d1,
                                hurst = Hlogistic, center = TRUE)
  return(dt_)
}, mc.cores = 75))
saveRDS(object = dt_mc_locreg_proxy_centered, file = "./inst/10_mc_simulate_data/data/dt_mc_far_proxy_centered_N=400_lambda=300.RDS")
rm(dt_mc_locreg_proxy_centered) ; gc()


