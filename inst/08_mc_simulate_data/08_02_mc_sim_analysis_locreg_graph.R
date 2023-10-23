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


# Estimates of the local regularity parameters ----
## Plot function
gplot <- function(N = 400, lambda = 300, design = "d1", param = "Ht", center = TRUE){
  if (center) {
    dt_locreg <- readRDS(paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
    ylim_H <- c(0.1, 0.9)
    ylim_L <- c(2, 15)
  } else {
    dt_locreg <- readRDS(paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_locreg_N=", N, "_lambda=", lambda, "_not_centered", design, ".RDS"))
    ylim_H <- c(0.1, 1.5)
    ylim_L <- c(0, 500)
  }
  dt_locreg[, t := as.factor(t)]
  # title_exp <- paste0("$\\Delta=e^{-\\log(\\lambda)^\\gamma}$, N = ", N , ", $\\lambda$=", lambda, ", $L_t^2 =4$")
  title_exp <- paste0("N = ", N , ", $\\lambda$=", lambda)
  if (param == "Ht") {
    gplt <- ggplot(dt_locreg, aes(x = t, y = get(param))) +
      geom_boxplot() +
      ylim(ylim_H) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("t") +
      ylab(latex2exp::TeX("$H_t$")) +
      geom_hline(
        yintercept = Hlogistic(t0),
        color = c("#839192", "#717D7E", "#616A6B", "#283747"),
        linetype = c("dashed", "dotted", "dotdash", "longdash")) +
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
    gplt <- ggplot(dt_locreg, aes(x = t, y = get(param))) +
      geom_boxplot() +
      ylim(ylim_L) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("t") +
      ylab(latex2exp::TeX("$L_t^2$")) +
      geom_hline(yintercept = 4, color = "#283747", linetype = "longdash") +
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

## Local regularity estimates on centered ----
### d1
g_d1 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d1", param = "Ht"),
  gplot(N = 1000, lambda = 1000, design = "d1", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d1", param = "Lt"),
  gplot(N = 1000, lambda = 1000, design = "d1", param = "Lt"),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_d1.png"),
  plot = g_d1, units = "px", dpi = 300, bg = "white")

### d2
g_d2 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d2", param = "Ht"),
  gplot(N = 1000, lambda = 1000, design = "d2", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d2", param = "Lt"),
  gplot(N = 1000, lambda = 1000, design = "d2", param = "Lt"),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_d2.png"),
  plot = g_d2, units = "px", dpi = 300, bg = "white")

### d3
g_d3 <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d3", param = "Ht"),
  gplot(N = 1000, lambda = 1000, design = "d3", param = "Ht"),
  gplot(N = 400, lambda = 300, design = "d3", param = "Lt"),
  gplot(N = 1000, lambda = 1000, design = "d3", param = "Lt"),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_d3.png"),
  plot = g_d3, units = "px", dpi = 300, bg = "white")

## Local regularity estimates on not centered ----
### d1
g_d1_not_centered <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d1", param = "Ht", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d1", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d1", param = "Lt", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d1", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_not_centered_d1.png"),
  plot = g_d1_not_centered, units = "px", dpi = 300, bg = "white")

### d2
g_d2_not_centered <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d2", param = "Ht", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d2", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d2", param = "Lt", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d2", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_not_centered_d2.png"),
  plot = g_d2_not_centered, units = "px", dpi = 300, bg = "white")

### d3
g_d3_not_centered <- gridExtra::grid.arrange(
  gplot(N = 400, lambda = 300, design = "d3", param = "Ht", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d3", param = "Ht", center = FALSE),
  gplot(N = 400, lambda = 300, design = "d3", param = "Lt", center = FALSE),
  gplot(N = 1000, lambda = 1000, design = "d3", param = "Lt", center = FALSE),
  ncol = 2, nrow = 2
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/locreg_not_centered_d3.png"),
  plot = g_d3_not_centered, units = "px", dpi = 300, bg = "white")

# Mean function estimation ----
far_mean_d1 <- function(t){
  res <- 4 * sin(3 * pi * t / 2)
  return(res)
}
## Plot function
gplot_mean <- function(N = 400, lambda = 300, design = "d1"){
  dt_mean <- readRDS(paste0("./inst/08_mc_simulate_data/locreg_estimates/dt_mean_N=", N, "_lambda=", lambda, "_", design, ".RDS"))
  dt_mean[, t := as.factor(t)]
  # title_exp <- paste0("$\\Delta=e^{-\\log(\\lambda)^\\gamma}$, N = ", N , ", $\\lambda$=", lambda, ", $L_t^2 =4$")
  title_exp <- paste0("Mean function - N = ", N , ", $\\lambda$=", lambda)
  if (design == "d1") {

    gplt <- ggplot(dt_mean, aes(x = t, y = muhat)) +
      geom_boxplot() +
      ylim(-4, 5) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("t") +
      ylab(latex2exp::TeX("$\\mu_N(t)$")) +
      geom_hline(
        yintercept = far_mean_d1(t0),
        color = c("#839192", "#717D7E", "#616A6B", "#283747"),
        linetype = c("dashed", "dotted", "dotdash", "longdash")) +
      theme_minimal() +
      theme(axis.title = element_text(size = 12),
            axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, margin = margin(t = , r = 10, b = 0, l = 0)),
            axis.text.x =  element_text(size = 12),
            axis.text.y =  element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.key.width= unit(0.8, 'cm'))

  } else {
    gplt <- ggplot(dt_mean, aes(x = t, y = muhat)) +
      geom_boxplot() +
      ylim(237, 247) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab("t") +
      ylab(latex2exp::TeX("$\\mu_N(t)$")) +
      geom_hline(
        yintercept = get_real_data_mean(t0),
        color = c("#839192", "#717D7E", "#616A6B", "#283747"),
        linetype = c("dashed", "dotted", "dotdash", "longdash")) +
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
  gplot_mean(N = 400, lambda = 300, design = "d1"),
  gplot_mean(N = 1000, lambda = 1000, design = "d1"),
  ncol = 2, nrow = 1
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/mean_d1.png"),
  plot = g_d1, units = "px", dpi = 300, bg = "white")

## d2
g_d2 <- gridExtra::grid.arrange(
  gplot_mean(N = 400, lambda = 300, design = "d2"),
  gplot_mean(N = 1000, lambda = 1000, design = "d2"),
  ncol = 2, nrow = 1
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/mean_d2.png"),
  plot = g_d2, units = "px", dpi = 300, bg = "white")

## d3
g_d3 <- gridExtra::grid.arrange(
  gplot_mean(N = 400, lambda = 300, design = "d3"),
  gplot_mean(N = 1000, lambda = 1000, design = "d3"),
  ncol = 2, nrow = 1
)
ggsave(
  filename = file.path("./inst/08_mc_simulate_data/locreg_estimates/mean_d3.png"),
  plot = g_d3, units = "px", dpi = 300, bg = "white")
