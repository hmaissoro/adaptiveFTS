library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
# Hlogistic <- function(t){
#   hurst_logistic(t, h_left = 0.4, h_right = 0.6,
#                  change_point_position = 0.5, slope = 50)
# }

Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

# Local regularity graph function ----

ggplot_locreg <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht", Hfun = Hlogistic){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(file_name)

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'),
          legend.position = "top")

  if (white_noise == "mfBm") {
    ## define segment and set scale label
    dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(t - 0.05),
                                  "xend" = as.factor(t + 0.049), "Htrue" =  Hfun(t))])
    scale_label <- c(dt_pr[, t], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(scale_label))
    scale_label[- which(scale_label %in% as.character(dt_pr[, t]))] <- as.character(" ")

    ## set t as factor
    dt_locreg[, t := as.factor(t)]
    if (param == "Ht") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{H}_t$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.1, 0.9)
      x_lab <- "t"
      y_lab <- ""
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = "twodash", linewidth = 0.9)
    } else if (param == "Lt") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)*
      title_exp <- paste0("$\\widehat{L}_t^2 - $", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 12)
      x_lab <- "t"
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = "twodash", linewidth = 0.9)
      y_lab <- ""
      scale_label <- scale_label[! scale_label == " "]
    }

    ggplt <- ggplot(data = dt_locreg, mapping = aes(x = t, y = get(param))) +
      geom_boxplot() +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      geom_theme

  } else if (white_noise == "fBm") {
    ## define segment and set scale label
    dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(Htrue - 0.05),
                                  "xend" = as.factor(Htrue + 0.05), Htrue)])
    scale_label <- c(dt_pr[, as.factor(Htrue)], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(unique(scale_label)))
    scale_label[- which(scale_label %in% as.character(dt_pr[, as.factor(Htrue)]))] <- ""

    ## set t and Htrue as factor
    dt_locreg[, t := as.factor(t)]
    dt_locreg[, Htrue := as.factor(Htrue)]

    if (param == "Ht") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{H}_t$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.1, 0.9)
      x_lab <- latex2exp::TeX("True $H_t$")
      y_lab <- ""
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = "twodash", linewidth = 0.9)
    } else if (param == "Lt") {
      # title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      title_exp <- paste0("$\\widehat{L}_t^2$ - ", "N=", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 12)
      x_lab <- latex2exp::TeX("True $H_t$")
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = "twodash", linewidth = 0.9)
      y_lab <- ""
      scale_label <- scale_label[! scale_label == ""]
    }

    ggplt <- ggplot(data = dt_locreg, mapping = aes(x = Htrue, y = get(param), fill = t)) +
      geom_boxplot() +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      scale_fill_grey() +
      geom_theme
  }
  # Save plots
  # plot_name <- paste0("./inst/12_mc_simulate_data/graphs/locreg/locreg_estimates_",
  #                     process,"_", white_noise, "_", param, "_N=", N, "_lambda=", lambda, "_", design,".png")
  # ggsave(filename = plot_name, plot = ggplt,
  #        width = 7.5 / 1.5, height = 5.97 / 1.5, units = "in", dpi = 300, bg = "white")
  return(ggplt)
}

# Plot local regularity parameters ----
# ## design 1 ----
# g_locreg_far_fma_mfBm_d1 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_mfBm_d1.png", plot = g_locreg_far_fma_mfBm_d1,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)
#
# g_locreg_far_fma_fBm_d1 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_fBm_d1.png", plot = g_locreg_far_fma_fBm_d1,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)
#
# ## design 2 ----
# g_locreg_far_fma_mfBm_d2 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_mfBm_d2.png", plot = g_locreg_far_fma_mfBm_d2,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)
#
# g_locreg_far_fma_fBm_d2 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_fBm_d2.png", plot = g_locreg_far_fma_fBm_d2,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)
#
# ## design 3 ----
# g_locreg_far_fma_mfBm_d3 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_mfBm_d3.png", plot = g_locreg_far_fma_mfBm_d3,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)
#
# g_locreg_far_fma_fBm_d3 <- gridExtra::grid.arrange(
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Ht"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Lt"),
#   ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Lt"),
#   nrow = 2, ncol = 2
# )
# ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg/locreg_far_fma_fBm_d2.png", plot = g_locreg_far_fma_fBm_d3,
#        width = 9.98, height = 8.5, units = "in", dpi = 300)


