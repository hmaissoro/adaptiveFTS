library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Simulation parameters ----
## Simulation global parameters----
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hvec <- c(0.4, 0.5, 0.7)

# Estimates of the local regularity parameters ----

## Plot function
gplot <- function(N = 400, lambda = 300, design = "d1", param = "Ht", center = TRUE){
  if (center) {
    dt_locreg <- readRDS(paste0("./inst/09_mc_simulate_data/locreg_estimates/dt_locreg_fBm_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
    title_exp <- paste0("N = ", N , ", $\\lambda$=", lambda, " - Centered")
    ylim_Ht <- c(0.1, 0.9)
    ylim_Lt <- c(2, 16)
  } else {
    dt_locreg <- readRDS(paste0("./inst/09_mc_simulate_data/locreg_estimates/dt_locreg_fBm_N=", N, "_lambda=", lambda, "_not_centered_", design, ".RDS"))
    title_exp <- paste0("N = ", N , ", $\\lambda$=", lambda, " - Not centered")
    ylim_Ht <- c(0.1, 1.5)
    ylim_Lt <- c(2, 500)
  }
  dt_locreg[, t := as.factor(t)]
  dt_locreg[, Htrue := as.factor(Htrue)]

  if (param == "Ht") {
    gplt <- ggplot(dt_locreg, aes(x = Htrue, y = Ht, fill = t)) +
      geom_boxplot() +
      ylim(ylim_Ht) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("H true") +
      ylab(latex2exp::TeX("$H_t$")) +
      geom_hline(
        yintercept = Hvec,
        color = c("#839192", "#717D7E", "#616A6B"),
        linetype = c("dashed", "dotted", "dotdash")) +
      scale_fill_grey() +
      theme_minimal() +
      theme(axis.title = element_text(size = 12),
            axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, margin = margin(t = , r = 10, b = 0, l = 0)),
            axis.text.x =  element_text(size = 12),
            axis.text.y =  element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.key.width= unit(0.8, 'cm'))
  } else if (param == "Lt") {
    gplt <- ggplot(dt_locreg, aes(x = Htrue, y = Lt, fill = t)) +
      geom_boxplot() +
      ylim(2, max(dt_locreg[, Lt])) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("H true") +
      ylab(latex2exp::TeX("$L_t^2$")) +
      geom_hline(yintercept = 4, color = "#283747", linetype = "longdash") +
      scale_fill_grey() +
      theme_minimal() +
      theme(axis.title = element_text(size = 12),
            axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, margin = margin(t = , r = 10, b = 0, l = 0)),
            axis.text.x =  element_text(size = 12),
            axis.text.y =  element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.key.width= unit(0.8, 'cm'))
  }
  return(gplt)
}

## d1
g_d1 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d1", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d1", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d1", param = "Lt"),
  gplot(N = 400, lambda = 300, design = "d1", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/09_mc_simulate_data/locreg_estimates/locreg_fBm_d1.png"),
  plot = g_d1, bg = "white")

## d2
g_d2 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d2", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d2", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d2", param = "Lt"),
  gplot(N = 400, lambda = 300, design = "d2", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/09_mc_simulate_data/locreg_estimates/locreg_fBm_d2.png"),
  plot = g_d2, bg = "white")

## d3
g_d3 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d3", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d3", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d3", param = "Lt"),
  gplot(N = 400, lambda = 300, design = "d3", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/09_mc_simulate_data/locreg_estimates/locreg_fBm_d3.png"),
  plot = g_d2, bg = "white")
