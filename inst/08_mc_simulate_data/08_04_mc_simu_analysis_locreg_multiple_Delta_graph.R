library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Simulation parameters ----
sig <- 0.5
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

Hlogistic(t0)

# Estimates of the local regularity parameters ----

## Plot function
gplot <- function(N = 400, lambda = 300, ti = 0.2, design = "d1", param = "Ht"){
  dt_locreg <- readRDS(paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_Delta_Vec_N=", N, "_lambda=", lambda, "_not_centered_", design, ".RDS"))
  dt_locreg_proxy <- readRDS(paste0("./inst/08_mc_simulate_data/data/dt_mc_locreg_proxy_far_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
  dt_locreg <- data.table::merge.data.table(
    x = dt_locreg, y = dt_locreg_proxy,
    by.x = c("id_mc", "N", "lambda", "t", "Delta"),
    by.y = c("id_mc", "N", "lambda", "t0", "Delta")
  )
  dt_locreg[, Htrue := Hlogistic(t)]
  dt_locreg[, ecartH := H - Htrue]
  dt_locreg[, ecartL := L - 4]
  dt_locreg <- dt_locreg[t == ti]
  dt_locreg[, t := as.factor(t)]
  dt_locreg[, Delta := as.factor(Delta)]
  # title_exp <- paste0("$\\Delta=e^{-\\log(\\lambda)^\\gamma}$, N = ", N , ", $\\lambda$=", lambda, ", $L_t^2 =4$")
  title_exp <- paste0("N = ", N , ", $\\lambda$=", lambda, " - t = ", ti)
  if (param == "Ht") {
    gplt <- ggplot(dt_locreg, aes(x = Delta, y = ecartH)) +
      geom_boxplot() +
      ylim(-0.2, 1.5) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("Delta") +
      ylab(latex2exp::TeX("$\\widetilde{H}_t - H_t$")) +
      geom_hline(yintercept = 0, color = "black", linetype = "longdash") +
      theme_minimal() +
      scale_fill_grey() +
      theme(axis.title = element_text(size = 9),
            axis.title.x = element_text(size = 9, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 9, margin = margin(t = , r = 10, b = 0, l = 0)),
            axis.text.x =  element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
            axis.text.y =  element_text(size = 9),
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 9),
            legend.key.width= unit(0.8, 'cm'), legend.position = "top")
  } else if (param == "Lt") {
    gplt <- ggplot(dt_locreg, aes(x = Delta, y = ecartL)) +
      geom_boxplot() +
      ggtitle(latex2exp::TeX(title_exp)) +
      ylim(-30, 500) +
      xlab("Delta") +
      ylab(latex2exp::TeX("$\\widetilde{L}_t^2 - L_t^2$")) +
      geom_hline(yintercept = 0, color = "black", linetype = "longdash") +
      theme_minimal() +
      scale_fill_grey() +
      theme(axis.title = element_text(size = 9),
            axis.title.x = element_text(size = 9, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 9, margin = margin(t = , r = 10, b = 0, l = 0)),
            axis.text.x =  element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
            axis.text.y =  element_text(size = 9),
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 9),
            legend.key.width= unit(0.8, 'cm'), legend.position = "top")
  }
  return(gplt)
}

## d1
for(i in 1:4){
  g_d1 <- gridExtra::grid.arrange(
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d1", param = "Ht"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d1", param = "Ht"),
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d1", param = "Lt"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d1", param = "Lt"),
    ncol = 2, nrow = 2
  )
  ggsave(
    filename = file.path(paste0("./inst/08_mc_simulate_data/locreg_estimates/locreg_delta_Htrue_t",t0[i] ,"_d1.png")),
    plot = g_d1, units = "px", dpi = 300, bg = "white")
}

## d2
for(i in 1:4){
  g_d2 <- gridExtra::grid.arrange(
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d2", param = "Ht"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d2", param = "Ht"),
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d2", param = "Lt"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d2", param = "Lt"),
    ncol = 2, nrow = 2
  )
  ggsave(
    filename = file.path(paste0("./inst/08_mc_simulate_data/locreg_estimates/locreg_delta_Htrue_t",t0[i] ,"_d2.png")),
    plot = g_d2, units = "px", dpi = 300, bg = "white")
}

## d3
for(i in 1:4){
  g_d3 <- gridExtra::grid.arrange(
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d3", param = "Ht"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d3", param = "Ht"),
    gplot(N = 400, lambda = 300, ti = t0[i], design = "d3", param = "Lt"),
    gplot(N = 1000, lambda = 1000, ti = t0[i], design = "d3", param = "Lt"),
    ncol = 2, nrow = 2
  )
  ggsave(
    filename = file.path(paste0("./inst/08_mc_simulate_data/locreg_estimates/locreg_delta_Htrue_t",t0[i] ,"_d3.png")),
    plot = g_d3, units = "px", dpi = 300, bg = "white")
}

