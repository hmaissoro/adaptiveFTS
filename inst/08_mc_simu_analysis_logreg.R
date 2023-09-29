library(data.table)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(latex2exp)

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

# Plot true local regularity parameters ----
## Generate data
dt_Hlogistic <- data.table::data.table(
  "t" = seq(0.01, 0.99, len = 100),
  "Ht" = Hlogistic(t = seq(0.01, 0.99, len = 100))
)
dt_Lconstant <- data.table::data.table(
  "t" = seq(0.01, 0.99, len = 100),
  "Lt" = rep(4, 100)
)

## Plot and save
figures_path <- "../../../report/learning-smmoothness/Learning-smoothness/figures/"
theme_set(theme_minimal())

### Local exponent
ggplot(dt_Hlogistic, aes(x = t, y = Ht)) +
  geom_line(size = 1.5) +
  ylim(0.4, 0.8) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$H_t$")) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "local_exponent.png"), units = "px", dpi = 300)

### Hölder function
ggplot(dt_Lconstant, aes(x = t, y = Lt)) +
  geom_line(size = 1.5) +
  ylim(3, 5) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$L_t^2$")) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "holder_constant.png"), units = "px", dpi = 300)

# Estimates of the local regularity parameters ----
dt_locreg[, c("N", "lambda") := .(as.factor(N), as.factor(lambda))]

ggplot(dt_locreg[t == t0[6] & delta_formula == "exponential"], aes(x = lambda, y = H_mc, fill = N)) +
  geom_boxplot() +
  ggtitle(paste0("Delta = exponential - t = ", t0[6] , " - H_t = ", round(Hlogistic(t0[6]), 3))) +
  geom_hline(aes(yintercept = Hlogistic(t0[6])), colour = "black") +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(filename = file.path("./inst/locreg_estimates/", "H_locreg_estimate_t6.png"), units = "px", dpi = 300)

## Plot Local exponent H
for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "exponential"], aes(x = lambda, y = H_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("Delta = exponential - t = ", t0[i] , " - H_t = ", round(Hlogistic(t0[i]), 3))) +
    geom_hline(aes(yintercept = Hlogistic(t0[i])), colour = "black") +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "H_locreg_Delta_exponential_estimate_t=",t0[i], ".png")), units = "px", dpi = 300)
}

for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "polynomial"], aes(x = lambda, y = H_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("Delta = polynomial - t = ", t0[i] , " - H_t = ", round(Hlogistic(t0[i]), 3))) +
    geom_hline(aes(yintercept = Hlogistic(t0[i])), colour = "black") +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "H_locregDelta_polynomial_estimate_t=",t0[i], ".png")), units = "px", dpi = 300)
}

## Hölder constant
for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "exponential"], aes(x = lambda, y = L_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("Delta = exponential - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    geom_hline(aes(yintercept = 4), colour = "black") +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "L_locreg_Delta_exponential_estimate_t=",t0[i], "_Lt.png")), units = "px", dpi = 300)
}

for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "polynomial"], aes(x = lambda, y = L_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("Delta = polynomial - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    geom_hline(aes(yintercept = 4), colour = "black") +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "L_locregDelta_polynomial_estimate_t=",t0[i], "_Lt.png")), units = "px", dpi = 300)
}

## Gamma exponent for Delta exponential

### Local exponent
for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "exponential"], aes(x = lambda, y = gH_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("gH - Delta = exponential - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "gH_Delta_exponential_estimate_t=",t0[i], "_Ht.png")), units = "px", dpi = 300)
}

for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "polynomial"], aes(x = lambda, y = gH_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("gH - Delta = polynomial - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "gH_Delta_polynomial_estimate_t=",t0[i], "_Ht.png")), units = "px", dpi = 300)
}

### Holder constant
for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "exponential"], aes(x = lambda, y = gL_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("gL - Delta = exponential - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "gL_Delta_exponential_estimate_t=",t0[i], "_Lt.png")), units = "px", dpi = 300)
}

for(i in 1:6){
  ggplot(dt_locreg[t == t0[i] & delta_formula == "polynomial"], aes(x = lambda, y = gL_mc, fill = N)) +
    geom_boxplot() +
    ggtitle(paste0("gL - Delta = polynomial - t = ", t0[i] , " - L_t = 4 - H_t = ", round(Hlogistic(t0[i]), 3))) +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top")
  ggsave(filename = file.path(paste0("./inst/locreg_estimates/", "gL_Delta_polynomial_estimate_t=",t0[i], "_Lt.png")), units = "px", dpi = 300)
}

