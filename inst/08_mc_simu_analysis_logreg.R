library(data.table)
library(ggplot2)
library(ggpubr)
library(magrittr)

# Simulation parameters ----
N <- c(100L, 200L, 400L)
lambda <- c(25L, 50L, 100L, 200L, 300L)
sig <- 0.5
mc <- 500
t0 <- seq(0.2, 0.8, len = 6)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

# Estimate local regularity ----
dt_locreg <- data.table::rbindlist(lapply(N, function(Ni){
  dt_lbda <- data.table::rbindlist(lapply(lambda, function(lambdai, Ni){
    # Import data
    dt <- readRDS(paste0("./inst/data/dt_mc_fts_real_data_N=", Ni, "_mu=", lambdai, "_Hlogistic_sig05.RDS"))
    dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "far_mean")]

    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(seq_len(mc), function(mc_i, Ni, lambdai){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      dt_mc <- dt_mc[order(id_curve, tobs)]
      bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]

      # Estimate the local regularity
      ## Take two a grid of Delta
      ## \Delta = exp(-log(\lambda)^\gamma), \gamma \in (0, 1) and
      ## \Delta = (1 / \lambda)^\gamma, \gamma \in (0, 1)
      lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
      g <- seq(0, 0.99, len = 100)
      g <- g[-1]
      dgrid_expo <- exp(- log(lambdahat) ** g)
      dgrid_poly <- (1 / lambdahat) ** g

      ## For exponential Delta
      dt_delta_expo <- data.table::rbindlist(lapply(1:length(dgrid_expo), function(idd, dgrid, gg, dt_mc_i, t0){
        dt_locreg <- estimate_locreg(
          data = dt_mc_i, idcol = "id_curve",
          tcol = "tobs", ycol = "X",
          t = t0, Delta = dgrid[idd], h = bw,
          smooth_ker = epanechnikov)
        dt_locreg[, g := gg[idd]]
      }, dgrid = dgrid_expo, gg = g, dt_mc_i = dt_mc, t0 = t0))

      ### Get best local regularity parameters and the Deltas which allow to obtain them
      dt_delta_expo[, delta_formula := "exponential"]
      dt_delta_expo[, c("Htrue", "Ltrue") := .(Hlogistic(t), 4)]
      dt_delta_expo[, c("H_mc", "DeltaH_mc", "gH_mc") := .(Ht[which.min(abs(Ht - Htrue))], Delta[which.min(abs(Ht - Htrue))], g[which.min(abs(Ht - Htrue))]), by = "t"]
      dt_delta_expo[, c("L_mc", "DeltaL_mc", "gL_mc") := .(Lt[which.min(abs(Lt - Ltrue))], Delta[which.min(abs(Lt - Ltrue))], g[which.min(abs(Lt - Ltrue))]), by = "t"]
      dt_delta_expo <- unique(dt_delta_expo[, .SD, .SDcols = c("t", "locreg_bw", "delta_formula", "DeltaH_mc", "DeltaL_mc", "gH_mc", "gL_mc", "Htrue", "H_mc", "Ltrue", "L_mc")])

      ## For polynomial Delta
      dt_delta_poly <- data.table::rbindlist(lapply(1:length(dgrid_poly), function(idd, dgrid, gg, dt_mc_i, t0){
        dt_locreg <- estimate_locreg(
          data = dt_mc_i, idcol = "id_curve",
          tcol = "tobs", ycol = "X",
          t = t0, Delta = dgrid[idd], h = bw,
          smooth_ker = epanechnikov)
        dt_locreg[, g := gg[idd]]
      }, dgrid = dgrid_poly, gg = g, dt_mc_i = dt_mc, t0 = t0))

      ### Get best local regularity parameters and the Deltas which allow to obtain them
      dt_delta_poly[, delta_formula := "polynomial"]
      dt_delta_poly[, c("Htrue", "Ltrue") := .(Hlogistic(t), 4)]
      dt_delta_poly[, c("H_mc", "DeltaH_mc", "gH_mc") := .(Ht[which.min(abs(Ht - Htrue))], Delta[which.min(abs(Ht - Htrue))], g[which.min(abs(Ht - Htrue))]), by = "t"]
      dt_delta_poly[, c("L_mc", "DeltaL_mc", "gL_mc") := .(Lt[which.min(abs(Lt - Ltrue))], Delta[which.min(abs(Lt - Ltrue))], g[which.min(abs(Lt - Ltrue))]), by = "t"]
      dt_delta_poly <- unique(dt_delta_poly[, .SD, .SDcols = c("t", "locreg_bw", "delta_formula", "DeltaH_mc", "DeltaL_mc", "gH_mc", "gL_mc", "Htrue", "H_mc", "Ltrue", "L_mc")])

      ## Merge the two estimates fof the local regularity parameters
      ## That dt_delta_expo and dt_delta_poly
      dt_delta <- rbind(dt_delta_expo, dt_delta_poly)
      rm(dt_delta_expo, dt_delta_poly) ; gc()

      ## Return obtained result
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_delta)
      return(dt_res)
    }, mc.cores = 75, Ni = Ni, lambdai = lambdai))
    return(dt_reg_mc)
  }, Ni = Ni))
  return(dt_lbda)
}))

saveRDS(dt_locreg, "./inst/locreg_estimates/dt_Hlogistic_bw_sig05.RDS")
