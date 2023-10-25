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

# local regularity estimation ----
estim_locreg_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      dt_mc <- dt_mc[order(id_curve, tobs)]
      bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
      bw_plus_mean <- unique(dt_mc[, .(id_curve, presmooth_bw_plus_mean)])[order(id_curve), presmooth_bw_plus_mean]

      # Estimate the local regularity
      ## Delta
      lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
      delta <- exp(- log(lambdahat) ** 0.35)

      ## Centered process
      dt_locreg <- estimate_locreg(
        data = dt_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = delta, h = bw,
        smooth_ker = epanechnikov, center = FALSE)

      ## process plus mean function
      dt_locreg_plus_mean <- estimate_locreg(
        data = dt_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X_plus_mean",
        t = t0, Delta = delta, h = bw_plus_mean,
        smooth_ker = epanechnikov, center = TRUE)
      data.table::setnames(x = dt_locreg_plus_mean,
                           old = c("locreg_bw", "Delta", "Nused", "Ht", "Lt"),
                           new = c("locreg_bw_plus_mean", "Delta_plus_mean", "Nused_plus_mean", "Ht_plus_mean", "Lt_plus_mean"))

      ## Merge local regularity parameters estimates and return the obtained result
      dt_locreg_res <- data.table::merge.data.table(x = dt_locreg, y = dt_locreg_plus_mean, by = "t")
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg_res)
      rm(dt_locreg_res, dt_locreg, dt_locreg_plus_mean) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt = dt, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai, t0){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      dt_mc <- dt_mc[order(id_curve, tobs)]
      Hvec <- dt_mc[, sort(unqiue(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_mc, t0){
        ## Set Delta
        lambdahat <- mean(dt_mc[Htrue == Hi, .N, by = id_curve][, N])
        delta <- exp(- log(lambdahat) ** 0.25)

        ## Extract bandwidth
        bw <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
        bw_plus_mean <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw_plus_mean)])[order(id_curve), presmooth_bw_plus_mean]

        ## Centered process
        dt_locreg <- estimate_locreg(
          data = dt_mc[Htrue == Hi], idcol = "id_curve",
          tcol = "tobs", ycol = "X",
          t = t0, Delta = delta, h = bw,
          smooth_ker = epanechnikov, center = FALSE)

        ## process plus mean function
        dt_locreg_plus_mean <- estimate_locreg(
          data = dt_mc[Htrue == Hi], idcol = "id_curve",
          tcol = "tobs", ycol = "X_plus_mean",
          t = t0, Delta = delta, h = bw_plus_mean,
          smooth_ker = epanechnikov, center = TRUE)
        data.table::setnames(x = dt_locreg_plus_mean,
                             old = c("locreg_bw", "Delta", "Nused", "Ht", "Lt"),
                             new = c("locreg_bw_plus_mean", "Delta_plus_mean", "Nused_plus_mean", "Ht_plus_mean", "Lt_plus_mean"))

        ## Merge local regularity parameters estimates and return the obtained result
        dt_locreg_res <- data.table::merge.data.table(x = dt_locreg, y = dt_locreg_plus_mean, by = "t")
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg_res)
        dt_res[, Htrue := Hi]

        ## Clean
        rm(dt_locreg_res, dt_locreg, dt_locreg_plus_mean) ; gc()
      }, dt_mc = dt_mc, t0 = t0))

      return(dt_res)
    }, mc.cores = 75, dt = dt, Ni = N, lambdai = lambda, t0 = t0))
  }


  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_reg_mc, file = file_name)
  rm(dt, file_title, file_name, dt_reg_mc) ; gc() ; gc()

  return(paste0("Done : dt_locreg_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS"))
}

# Estimate local regularity ----
## FAR ----
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1")

## FMA ----
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1")
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1")
