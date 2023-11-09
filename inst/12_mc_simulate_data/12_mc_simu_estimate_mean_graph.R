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

ggplot_mean_risk <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  ## Load data, remove NaN values and reshape data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean_risk <- readRDS(file_name)
  dt_mean_risk <- dt_mean_risk[! (is.nan(mutitle_mse) | is.nan(mean_risk))]

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 10, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 9, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 9, margin = margin(t = , r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 9),
          axis.text.y =  element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    dt_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk),
                                "mutitle_mse" = mean(mutitle_mse)),
                            by = c("t", "h")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("mean_risk", "mutitle_mse"),
                                     id.vars = c("t", "h"))

    fplot <- function(ti){
        ggplot(data = dt_risk_melt[variable %in% c("mutitle_mse", "mean_risk") & t == ti],
               mapping = aes(x = h, y = risk, group = variable, linetype = variable)) +
          geom_line(linewidth = 0.8) +
          ylim(0, 0.3) +
          labs(title = paste0("t = ", ti, " - H = ", round(Hlogistic(ti), 2)), y = "", x = "h") +
          scale_linetype_manual(values = c("mutitle_mse" = "solid", "mean_risk" = "twodash"),
                                labels = c("mutitle_mse" = "MSE", "mean_risk" = latex2exp::TeX("  $R_\\mu(h, \\widehat{H}_t, \\widehat{L_t^2})$")),
                                name = "Risk") +
          geom_theme
    }

    # Builb plots
    ggt1 <- fplot(ti = 0.2)
    ggt2 <- fplot(ti = 0.4)
    ggt3 <- fplot(ti = 0.7)
    ggt4 <- fplot(ti = 0.8)
    gplot <- ggpubr::ggarrange(ggt1, ggt2, ggt3, ggt4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

    # Save plots
    plot_name <- paste0("./inst/12_mc_simulate_data/graphs/mean/mean_risk_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".png")
    ggsave(filename = plot_name, plot = gplot,
           width = 7.5, height = 5.97, units = "in", dpi = 300, bg = "white")
    rm(ggt1, ggt2, ggt3, ggt4)

  } else if (white_noise == "fBm") {
    dt_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk),
                                "mutitle_mse" = mean(mutitle_mse)),
                            by = c("t", "h", "Htrue")]
    dt_risk_melt <- data.table::melt(data = dt_risk, value.name = "risk",
                                     measure.vars = c("mean_risk", "mutitle_mse"),
                                     id.vars = c("t", "h", "Htrue"))

    fplot <- function(ti, Hi = 0.4){
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
    for(Hi in dt_risk_melt[, unique(Htrue)]){
      # build plots
      ggt1 <- fplot(ti = 0.2, Hi = Hi)
      ggt2 <- fplot(ti = 0.4, Hi = Hi)
      ggt3 <- fplot(ti = 0.7, Hi = Hi)
      ggt4 <- fplot(ti = 0.8, Hi = Hi)
      gplot <- ggpubr::ggarrange(ggt1, ggt2, ggt3, ggt4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")

      # Save plots
      plot_name <- paste0("./inst/12_mc_simulate_data/graphs/mean/mean_risk_",
                          process,"_", white_noise, "_H=", Hi, "_N=", N, "_lambda=", lambda, "_", design,".png")
      ggsave(filename = plot_name, plot = gplot,
             width = 7.5, height = 5.97, units = "in", dpi = 300, bg = "white")
      rm(ggt1, ggt2, ggt3, ggt4, gplot)
    }
  }
  # return(gplot)
  return(paste0("Done : plot for ", process," with ", white_noise, " for ", "N=", N, " and lambda=", lambda, " for ", design))
}

# Plot local regularity parameters ----
## design 1 ----
### FAR ----
ggplot_mean_risk(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
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


# Mean function estimates graph ----
ggplot_mean <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_mean <- readRDS(file_name)

  ## Add process true mean
  data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(data_file_name)
  dt_mean <- data.table::merge.data.table(x = dt_mean, y = unique(dt[, .("t" = tobs, "mutrue" = process_mean)]), by = "t")
  rm(dt) ; gc()

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(plot.title = element_text(size = 9),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 12, margin = margin(t = 10, r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'),
          legend.position = "top")

  if (white_noise == "mfBm") {
    ## define segment and set scale label
    dt_pr <- unique(dt_mean[, .("t" = as.factor(t), "x" = as.factor(t - 0.05),
                                  "xend" = as.factor(t + 0.05), "mutrue" =  mutrue)])
    scale_label <- c(dt_pr[, t], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(scale_label))
    scale_label[- which(scale_label %in% as.character(dt_pr[, t]))] <- ""

    ## set t as factor
    dt_mean[, t := as.factor(t)]

    title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
    # y_lim <- c(0.2, 0.9)
    x_lab <- "t"
    y_lab <-  latex2exp::TeX("$\\widehat{\\mu}(t)$")
    geom_true_param <- geom_segment(
      data = dt_pr, mapping = aes(x = x, xend = xend, y = mutrue, yend = mutrue),
      linetype = 2)
    ggplt <- ggplot(data = dt_mean, mapping = aes(x = t, y = muhat)) +
      geom_boxplot() +
      # ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      geom_theme

  } else if (white_noise == "fBm") {

    ## define segment and set scale label
    dt_pr <- unique(dt_mean[, .("t" = as.factor(t), "x" = as.factor(t - 0.06),
                                "xend" = as.factor(t + 0.06), "mutrue" =  mutrue)])
    scale_label <- c(dt_pr[, t], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(scale_label))
    scale_label[- which(scale_label %in% as.character(dt_pr[, t]))] <- ""

    ## set t and Htrue as factors
    dt_mean[, t := as.factor(t)]
    dt_mean[, Htrue := as.factor(Htrue)]

    title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
    # y_lim <- c(0.2, 0.9)
    x_lab <- "t"
    y_lab <-  latex2exp::TeX("$\\widehat{\\mu}(t)$")
    geom_true_param <- geom_segment(
      data = dt_pr, mapping = aes(x = x, xend = xend, y = mutrue, yend = mutrue),
      linetype = 2)
    ggplt <- ggplot(data = dt_mean, mapping = aes(x = t, y = muhat, fill = Htrue, color = Htrue)) +
      geom_boxplot() +
      # ylim(y_lim) +
      ggtitle(latex2exp::TeX(title_exp)) +
      xlab(x_lab) +
      ylab(y_lab) +
      geom_true_param +
      scale_x_discrete(labels = scale_label) +
      geom_theme

    ## define segment and set scale label
    dt_pr <- unique(dt_mean[, .("t" = as.factor(t), "x" = as.factor(t - 0.05),
                                "xend" = as.factor(t + 0.05), Htrue)])
    scale_label <- c(dt_pr[, as.factor(Htrue)], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(unique(scale_label)))
    scale_label[- which(scale_label %in% as.character(dt_pr[, as.factor(Htrue)]))] <- ""

    ## set t and Htrue as factor
    dt_mean[, t := as.factor(t)]
    dt_mean[, Htrue := as.factor(Htrue)]

    if (param == "Ht") {
      title_exp <- paste0("Zero-mean ", process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.2, 0.9)
      x_lab <- latex2exp::TeX("True $H_t$")
      y_lab <- latex2exp::TeX("$\\widehat{H}_t$")
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = 2)
    } else if (param == "Ht_plus_mean") {
      title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.2, 0.9)
      x_lab <- latex2exp::TeX("True $H_t$")
      y_lab <-  latex2exp::TeX("$\\widehat{H}_t$")
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = 2)
    } else if (param == "Lt") {
      title_exp <- paste0("Zero-mean ", process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 15)
      x_lab <- latex2exp::TeX("True $H_t$")
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = 2)
      y_lab <- latex2exp::TeX("$\\widehat{L}_t^2$")
      scale_label <- scale_label[! scale_label == ""]
    } else if (param == "Lt_plus_mean"){
      title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 15)
      x_lab <- latex2exp::TeX("True $H_t$")
      y_lab <- latex2exp::TeX("$\\widehat{L}_t^2$")
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = 2)
      scale_label <- scale_label[! scale_label == ""]
    }

    ggplt <- ggplot(data = dt_mean, mapping = aes(x = Htrue, y = get(param), fill = t)) +
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

  return(ggplt)
}

ggplot_mean(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")

