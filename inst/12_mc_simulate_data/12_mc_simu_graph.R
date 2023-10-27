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
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

# Local regularity graph function ----

ggplot_locreg <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt_locreg <- readRDS(file_name)

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
    dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(t - 0.05),
                                  "xend" = as.factor(t + 0.05), "Htrue" =  Hlogistic(t))])
    scale_label <- c(dt_pr[, t], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(scale_label))
    scale_label[- which(scale_label %in% as.character(dt_pr[, t]))] <- ""

    ## set t as factor
    dt_locreg[, t := as.factor(t)]
    if (param == "Ht") {
      title_exp <- paste0("Zero-mean ", process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.2, 0.9)
      x_lab <- "t"
      y_lab <- latex2exp::TeX("$\\widehat{H}_t$")
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = 2)
    } else if (param == "Ht_plus_mean") {
      title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(0.2, 0.9)
      x_lab <- "t"
      y_lab <-  latex2exp::TeX("$\\widehat{H}_t$")
      geom_true_param <- geom_segment(
        data = dt_pr, mapping = aes(x = x, xend = xend, y = Htrue, yend = Htrue),
        linetype = 2)
    } else if (param == "Lt") {
      title_exp <- paste0("Zero-mean ", process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 15)
      x_lab <- "t"
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = 2)
      y_lab <- latex2exp::TeX("$\\widehat{L}_t^2$")
      scale_label <- scale_label[! scale_label == ""]
    } else if (param == "Lt_plus_mean") {
      title_exp <- paste0(process, "(1) - WN = ", white_noise, " - N = ", N , ", $\\lambda$=", lambda)
      y_lim <- c(-2, 15)
      x_lab <- "t"
      y_lab <- latex2exp::TeX("$\\widehat{L}_t^2$")
      geom_true_param <- geom_hline(yintercept = 4, color = "#283747", linetype = 2)
      scale_label <- scale_label[! scale_label == ""]
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
    dt_pr <- dt_pr <- unique(dt_locreg[, .("t" = as.factor(t), "x" = as.factor(Htrue - 0.05),
                                           "xend" = as.factor(Htrue + 0.05), Htrue)])
    scale_label <- c(dt_pr[, as.factor(Htrue)], dt_pr[, x], dt_pr[, xend])
    scale_label <- sort(as.character(unique(scale_label)))
    scale_label[- which(scale_label %in% as.character(dt_pr[, as.factor(Htrue)]))] <- ""

    ## set t and Htrue as factor
    dt_locreg[, t := as.factor(t)]
    dt_locreg[, Htrue := as.factor(Htrue)]

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

  return(ggplt)
}

# Plot local regularity parameters ----
## FAR ----
## d1
g_locreg_far_mfBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)

g_locreg_far_fBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)

## FMA ----
## d1
g_locreg_fma_mfBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)

g_locreg_fma_fBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)


