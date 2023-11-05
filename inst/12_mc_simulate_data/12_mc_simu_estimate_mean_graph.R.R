library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Simulation global parameters----
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

# Local regularity graph function ----

ggplot_mean_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "mean_risk_plus_mean"){
  ## Load data, remove NaN values and reshape data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean_risk <- readRDS(file_name)
  dt_mean_risk <- dt_mean_risk[! (is.nan(mutitle_mse) | is.nan(mean_risk) | is.nan(mean_risk_plus_mean))]

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 10, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 10, margin = margin(t = , r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    dt_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk),
                                "mean_risk_plus_mean" = mean(mean_risk_plus_mean),
                                "mutitle_mse" = mean(mutitle_mse)),
                            by = c("t", "h")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("mean_risk", "mean_risk_plus_mean", "mutitle_mse"),
                                     id.vars = c("t", "h"))

    fplot <- function(ti, param = param){
      if (param == "mean_risk_plus_mean") {
        ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", param) & t == ti],
               mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
          geom_line(linewidth = 0.8) +
          ylim(0, 0.9) +
          labs(title = paste0("t = ", ti, " - H = ", round(Hlogistic(ti), 2)), y = "", x = "h") +
          scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk_plus_mean" = "twodash"),
                                labels = c("mutitle_mse" = "MSE", "mean_risk_plus_mean" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                                name = "Risk") +
          geom_theme
      } else if (param == "mean_risk") {
        ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", param) & t == ti],
               mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
          geom_line(linewidth = 0.8) +
          ylim(0, 0.9) +
          labs(title = paste0("t = ", ti, " - H = ", round(Hlogistic(ti), 2)), y = "", x = "h") +
          scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk" = "twodash"),
                                labels = c("mutitle_mse" = "MSE", "mean_risk" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                                name = "Risk") +
          geom_theme
      }
    }

    # Builb plots
    ggt1 <- fplot(ti = 0.2, param = param)
    ggt2 <- fplot(ti = 0.4, param = param)
    ggt3 <- fplot(ti = 0.7, param = param)
    ggt4 <- fplot(ti = 0.8, param = param)
    gplot <- ggpubr::ggarrange(ggt1, ggt2, ggt3, ggt4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

    # Save plots
    plot_name <- paste0("./inst/12_mc_simulate_data/graphs/mean/", param,"_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".png")
    ggsave(filename = plot_name, plot = gplot,
           width = 6.81, height = 5.97, units = "in", dpi = 300, bg = "white")
    rm(ggt1, ggt2, ggt3, ggt4, gplot)

  } else if (white_noise == "fBm") {
    dt_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk),
                                "mean_risk_plus_mean" = mean(mean_risk_plus_mean),
                                "mutitle_mse" = mean(mutitle_mse)),
                            by = c("t", "h", "Htrue")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("mean_risk", "mean_risk_plus_mean", "mutitle_mse"),
                                     id.vars = c("t", "h", "Htrue"))

    fplot <- function(ti, Hi = 0.4, param = param){
      if (param == "mean_risk_plus_mean") {
        ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", param) & t == ti & Htrue == Hi],
               mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
          geom_line(linewidth = 0.8) +
          ylim(0, 0.9) +
          labs(title = paste0("t = ", ti, " - H = ", Hi), y = "", x = "h") +
          scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk_plus_mean" = "twodash"),
                                labels = c("mutitle_mse" = "MSE", "mean_risk_plus_mean" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                                name = "Risk") +
          geom_theme
      } else if (param == "mean_risk") {
        ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", param) & t == ti & Htrue == Hi],
               mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
          geom_line(linewidth = 0.8) +
          ylim(0, 0.9) +
          labs(title = paste0("t = ", ti, " - H = ", Hi), y = "", x = "h") +
          scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk" = "twodash"),
                                labels = c("mutitle_mse" = "MSE", "mean_risk" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                                name = "Risk") +
          geom_theme
      }
    }
    for(Hi in dt_risk_melt[, unique(Htrue)]){
      # build plots
      ggt1 <- fplot(ti = 0.2, Hi = Hi, param = param)
      ggt2 <- fplot(ti = 0.4, Hi = Hi, param = param)
      ggt3 <- fplot(ti = 0.7, Hi = Hi, param = param)
      ggt4 <- fplot(ti = 0.8, Hi = Hi, param = param)
      gplot <- ggpubr::ggarrange(ggt1, ggt2, ggt3, ggt4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

      # Save plots
      plot_name <- paste0("./inst/12_mc_simulate_data/graphs/mean/", param,"_",
                          process,"_", white_noise, "_H=", Hi, "_N=", N, "_lambda=", lambda, "_", design,".png")
      ggsave(filename = plot_name, plot = gplot,
             width = 6.81, height = 5.97, units = "in", dpi = 300, bg = "white")
      rm(ggt1, ggt2, ggt3, ggt4, gplot)
    }
  }

  return(paste0("Done : plot for ", process," with ", white_noise, " for ", "N=", N, " and lambda=", lambda, " for ", design))
}

# Plot local regularity parameters ----
## design 1 ----
### FAR ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "mean_risk_plus_mean")

### FMA ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "mean_risk_plus_mean")

## design 2 ----
### FAR ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "mean_risk_plus_mean")

### FMA ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "mean_risk_plus_mean")

## design 3 ----
### FAR ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "mean_risk_plus_mean")

### FMA ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "mean_risk_plus_mean")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "mean_risk")
ggplot_mean_risk(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "mean_risk_plus_mean")



