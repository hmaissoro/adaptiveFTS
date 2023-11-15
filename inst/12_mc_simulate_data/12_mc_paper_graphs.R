library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

source("./inst/12_mc_simulate_data/12_mc_simu_estimate_locreg_graph.R")

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

# Local regularity parameters ----

g_locreg_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  nrow = 2, ncol = 4)

ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_far_mfBm_d1.png", plot = g_locreg_far_mfBm_d1,
  width = 9, height = 6, units = "in", bg = "white")

g_locreg_far_fBm_d1  <- ggpubr::ggarrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")

ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/locreg_far_fBm_d1.png", plot = g_locreg_far_fBm_d1,
  width = 9, height = 6, units = "in", bg = "white")

# Estimate mean function ----

## Risk function of the mean function
ggplot_mean_risk_by_t <- function(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean_risk <- readRDS(file_name)
  dt_mean_risk <- dt_mean_risk[! (is.nan(mutitle_mse) | is.nan(mean_risk))]

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 9, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    dt_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk),
                                "mutitle_mse" = mean(mutitle_mse)),
                            by = c("t", "h")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("mean_risk", "mutitle_mse"),
                                     id.vars = c("t", "h"))
      title_exp <- paste0("t=", ti, " - ", "N=", N , ", $\\lambda$=", lambda)
      ggplt <- ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", "mean_risk") & t == ti],
                      mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
        geom_line(linewidth = 0.8) +
        ylim(0, 0.4) +
        ggtitle(latex2exp::TeX(title_exp)) +
        labs(y = "", x = "h") +
        scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk" = "twodash"),
                              labels = c("mutitle_mse" = "MSE", "mean_risk" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                              name = "Risk") +
        geom_theme

    return(ggplt)
  }
}
theme_bank_margin <- theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"))
theme_last_plot <- theme(
  plot.margin = unit(c(-8, 0, -30, 0), "pt"),
  axis.title.x = element_text(size = 10, margin = margin(t = 1, r = 0, b = 18, l = 0, unit = "pt"))
  )


g_mean_risk_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_mean_risk_by_t(N = 150, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 150, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 150, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 150, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 1000, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 1000, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 1000, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 1000, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 400, lambda = 300, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 400, lambda = 300, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 400, lambda = 300, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_risk_by_t(N = 1000, lambda = 1000, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_risk_by_t(N = 1000, lambda = 1000, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_risk_by_t(N = 1000, lambda = 1000, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_risk_by_t(N = 1000, lambda = 1000, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  nrow = 4, ncol = 4, common.legend = TRUE, legend = "bottom")
g_mean_risk_far_mfBm_d1

ggsave(filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_risk_far_mfBm_d1.png",
       plot = g_mean_risk_far_mfBm_d1, width = 9, height = 6, units = "in", bg = "white")

## Mean estimate of the risk function
g_mean_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_mean(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_far_mfBm_d1.png", plot = g_mean_far_mfBm_d1,
  width = 9, height = 6, units = "in", bg = "white")

g_mean_far_fBm_d1  <- ggpubr::ggarrange(
  ggplot_mean(N = 150, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1"),
  ggplot_mean(N = 1000, lambda = 40, process = "FAR", white_noise = "fBm", design = "d1"),
  ggplot_mean(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1"),
  ggplot_mean(N = 1000, lambda = 1000, process = "FAR", white_noise = "fBm", design = "d1"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_far_fBm_d1.png", plot = g_mean_far_fBm_d1,
  width = 9, height = 6, units = "in", bg = "white")

# Estimate autocovariance function ----

## Risk function of the autocovariance function
ggplot_autocov_risk_by_st <- function(N = 400, lambda = 300, si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_auto_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_autocov_risk <- readRDS(file_name)
  dt_autocov_risk <- dt_autocov_risk[! (is.nan(gammatilde_mse) | is.nan(autocov_risk))]

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 7, hjust = 0.7, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    dt_risk <- dt_autocov_risk[, .("autocov_risk" = mean(autocov_risk),
                                "gammatilde_mse" = mean(gammatilde_mse)),
                            by = c("s", "t", "h")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("autocov_risk", "gammatilde_mse"),
                                     id.vars = c("s", "t", "h"))
    title_exp <- paste0("(s,t)=(", si, ",", ti, ") - ", "(N,$\\lambda$)=(", N , ",", lambda, ")")
    # title_exp <- paste0("(N,$\\lambda$)=(", N , ",", lambda, ") $\\newline$", "(s,t)=(", si, ",", ti, ")")
    ggplt <- ggplot(data = dt_risk_melt[variable %in% c("gammatilde_mse", "autocov_risk") & s == si & t == ti],
                    mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
      geom_line(linewidth = 0.8) +
      ylim(0, 65) +
      ggtitle(latex2exp::TeX(title_exp)) +
      labs(y = "", x = "h") +
      scale_linetype_manual(values = c("gammatilde_mse" = "solid", "autocov_risk" = "twodash"),
                            labels = c("gammatilde_mse" = "MSE", "autocov_risk" = latex2exp::TeX("  $R_\\gamma(h, \\widehat{H}_s, \\widehat{H}_t, \\widehat{L_t^2}, \\widehat{L_s^2})$")),
                            name = "Risk") +
      geom_theme

    return(ggplt)
  }
}
theme_bank_margin <- theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           plot.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"))
theme_last_plot <- theme(
  plot.margin = unit(c(-8, 0, -30, 0), "pt"),
  axis.title.x = element_text(size = 10, margin = margin(t = 1, r = 0, b = 18, l = 0, unit = "pt"))
)


g_autocov_risk_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_autocov_risk_by_st(N = 150, lambda = 40, si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 150, lambda = 40, si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 150, lambda = 40, si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 150, lambda = 40, si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 40, si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 40, si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 40, si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 40, si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 400, lambda = 300, si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 400, lambda = 300, si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 400, lambda = 300, si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 400, lambda = 300, si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 1000, si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 1000, si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 1000, si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_autocov_risk_by_st(N = 1000, lambda = 1000, si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  nrow = 4, ncol = 4, common.legend = TRUE, legend = "bottom")
g_autocov_risk_far_mfBm_d1

ggsave(filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/autocov_risk_far_mfBm_d1.png",
       plot = g_autocov_risk_far_mfBm_d1, width = 9, height = 6, units = "in", bg = "white")


