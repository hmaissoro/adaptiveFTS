# Mean function estimates graph ----
ggplot_mean_qq_by_t <- function(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean <- readRDS(file_name)
  dt_mean <- unique(dt_mean)

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 9, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 10, margin = margin(t = , r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'))

  title_exp <- paste0("t=", ti, " - ", "(N, $\\lambda$)=(", N , ",", lambda, ")")
  dt_mean[, "mu" := sqrt(PN) * (muhat - mutrue), by = t]
  ggplt <-  ggplot(data = dt_mean[t == ti], aes(sample = mu)) +
    stat_qq() +
    stat_qq_line() +
    labs(y = "", x = "N(0,1)") +
    ggtitle(latex2exp::TeX(title_exp)) +
    geom_theme
  return(ggplt)
}

theme_bank_margin <- theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           plot.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"))
theme_last_plot <- theme(
  plot.margin = unit(c(-8, 0, -30, 0), "pt"),
  axis.title.x = element_text(size = 10, margin = margin(t = 1, r = 0, b = 18, l = 0, unit = "pt"))
)

g_mean_qq_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_bank_margin,
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1") + theme_last_plot,
  nrow = 4, ncol = 4, common.legend = TRUE, legend = "bottom")
g_mean_qq_far_mfBm_d1

ggsave(filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_d1.png",
       plot = g_mean_qq_far_mfBm_d1, width = 9, height = 6, units = "in", bg = "white")
