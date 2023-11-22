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
ggplot_mean_risk_by_t <- function(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
                                  ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  dt_risk <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
    file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                        process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
    dt_mean_risk <- readRDS(file_name)
    dt_mean_risk <- dt_mean_risk[! (is.nan(mutitle_mse) | is.nan(mean_risk))]
    dt_mean_risk <- dt_mean_risk[, .("mean_risk" = mean(mean_risk)), by = c("N", "lambda", "t", "h")]
    dt_mean_risk[, N_lambda := paste0("(", N, ",", lambda, ")")]
    return(dt_mean_risk)
  }))

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))

  if (white_noise == "mfBm") {
    title_exp <- paste0("t=", ti)
    ggplt <- ggplot(data = dt_risk[t == ti], mapping = aes(x = h, y = mean_risk, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
      geom_line(linewidth = 0.9) +
      ylim(0, 0.4) +
      ggtitle(latex2exp::TeX(title_exp)) +
      labs(y = "", x = "h") +
      scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                            values = c("(1000,1000)" = "solid", "(400,300)" = "dotted", "(1000,40)" = "dotdash", "(150,40)" = "dashed"),
                            labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                          values = c("(1000,1000)" = "#273746", "(400,300)" = "#34495E", "(1000,40)" = "#707B7C", "(150,40)" = "#909497"),
                          labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      geom_theme
    return(ggplt)
  }
}

g_mean_risk_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_mean_risk_by_t(ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_risk_by_t(ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_risk_by_t(ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_risk_by_t(ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
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


# Mean function table ----
table_mean <- function(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
                       process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  if (white_noise == "mfBm") {
    dt_mean <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
      file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                          process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
      dt_mean <- readRDS(file_name)
      dt_mean <- dt_mean[! (is.nan(muhat) | is.nan(muhat))]
      dt_mean <- dt_mean[, .("bias" = as.character(round(mean(muhat) - mean(mutrue), 4)),
                             "sd" = as.character(round(sd(muhat), 4))), by = c("N", "lambda", "t")]
      dt_mean[, N_lambda := paste0("(", N, ",", lambda, ")")]
      return(dt_mean)
    }))
    dt_mean_dcast <- data.table::dcast(data = dt_mean, formula = N + lambda ~ t, value.var = c("bias", "sd"))
    data.table::setcolorder(x = dt_mean_dcast,
                            neworder = c("N", "lambda", "bias_0.2", "sd_0.2", "bias_0.4", "sd_0.4",
                                         "bias_0.7", "sd_0.7", "bias_0.8", "sd_0.8"))
    dt_mean_dcast <- dt_mean_dcast[order(lambda)]

    # library("Hmisc")
    # library("htmlTable")
    # tab <- htmlTable::htmlTable(x = dt_mean_dcast, n.cgroup = c(1, 1, 2, 2, 2, 2) ,
    #                             cgroup = c("", "", "t = 0.2", "t = 0.4", "t = 0.7", "t = 0.8"),
    #                             header = c("$N$", "$\\lmabda$", "$Bias(\\widehat(\\mu))$", "$Sd(\\widehat(\\mu))$", "$Bias(\\widehat(\\mu))$",
    #                                        "$Sd(\\widehat())$", "$Bias(\\widehat(\\mu))$", "$Sd(\\widehat(\\mu))$", "$Bias(\\widehat(\\mu))$", "$Sd(\\widehat(\\mu))$"))
    Hmisc::latex(
      object = dt_mean_dcast,
      file =  paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/latex_mean_estimates_", process,"_", white_noise, "_", design,".tex"),
      n.cgroup = c(1, 1, 2, 2, 2, 2),
      cgroup = c("", "", "t = 0.2", "t = 0.4", "t = 0.7", "t = 0.8"),
      colheads = c("$N$", "$\\lambda$", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd"),
      rowname = NULL)
  } else if (white_noise == "fBm") {
    dt_mean <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
      file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_estimates_",
                          process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
      dt_mean <- readRDS(file_name)
      dt_mean <- dt_mean[! (is.nan(muhat) | is.nan(muhat))]
      dt_mean <- dt_mean[, .("bias" = as.character(round(mean(muhat) - mean(mutrue), 4)),
                             "sd" = as.character(round(sd(muhat), 4))), by = c("N", "lambda", "t", "Htrue")]
      dt_mean[, N_lambda := paste0("(", N, ",", lambda, ")")]
      return(dt_mean)
    }))
    dt_mean_dcast <- data.table::dcast(data = dt_mean, formula = Htrue + N + lambda ~ t, value.var = c("bias", "sd"))
    data.table::setcolorder(x = dt_mean_dcast,
                            neworder = c("Htrue", "N", "lambda", "bias_0.2", "sd_0.2", "bias_0.4", "sd_0.4",
                                         "bias_0.7", "sd_0.7", "bias_0.8", "sd_0.8"))
    dt_mean_dcast <- dt_mean_dcast[order(Htrue, lambda)]

    # library("Hmisc")
    # library("htmlTable")
    # htmlTable::htmlTable(x = dt_mean_dcast,
    #                      n.cgroup = c(1, 1, 1, 2, 2, 2, 2) ,
    #                      cgroup = c("", "", "", "t = 0.2", "t = 0.4", "t = 0.7", "t = 0.8"),
    #                      # n.rgroup = c(4, 4, 4),
    #                      # rgoup = c("$H_t = 0.4$", "$H_t = 0.5$", "$H_t = 0.7$"),
    #                      # rnames = c("$H_t = 0.4$", "$H_t = 0.5$", "$H_t = 0.7$"),
    #                      header = c("$H_t$", "$N$", "$\\almbda$", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd"))
    Hmisc::latex(
      object = dt_mean_dcast,
      file =  paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/latex_mean_estimates_", process,"_", white_noise, "_", design,".tex"),
      n.cgroup = c(1, 1, 1, 2, 2, 2, 2) ,
      cgroup = c("", "", "", "t = 0.2", "t = 0.4", "t = 0.7", "t = 0.8"),
      colheads = c("$H_t$", "$N$", "$\\lambda$", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd", "Bias", "Sd"),
      rowname = NULL)
  }


}

# Estimate autocovariance function ----

## Risk function of the autocovariance function
ggplot_autocov_risk_by_st <- function(N_vec = c(150, 1000, 400, 1000), lambda_vec = c(40, 40, 300, 1000),
                                      si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"){

  ## Load data, remove NaN values and reshape data
  dt_risk <- data.table::rbindlist(lapply(1:length(N_vec), function(i){
    ## Load data, remove NaN values and reshape data
    file_name <- paste0("./inst/12_mc_simulate_data/", process, "/autocov_estimates/dt_auto_risk_",
                        process,"_", white_noise, "_", "N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design,".RDS")
    dt_autocov_risk <- readRDS(file_name)
    dt_autocov_risk <- dt_autocov_risk[! (is.nan(gammatilde_mse) | is.nan(autocov_risk))]
    dt_autocov_risk <- dt_autocov_risk[, .("autocov_risk" = mean(autocov_risk)), by = c("N", "lambda", "s", "t", "h")]
    dt_autocov_risk[, N_lambda := paste0("(", N, ",", lambda, ")")]
    return(dt_autocov_risk)
  }))

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = -10),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))
  if (white_noise == "mfBm") {
    title_exp <- paste0("(s,t)=(", si, ",", ti, ")")
    ggplt <- ggplot(data = dt_risk[s == si & t == ti], mapping = aes(x = h, y = autocov_risk, group = N_lambda, linetype = N_lambda, color = N_lambda)) +
      geom_line(linewidth = 0.9) +
      ylim(0, 65) +
      ggtitle(latex2exp::TeX(title_exp)) +
      labs(y = "", x = "h") +
      scale_linetype_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                            values = c("(1000,1000)" = "solid", "(400,300)" = "dotted", "(1000,40)" = "dotdash", "(150,40)" = "dashed"),
                            labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      scale_colour_manual(name = latex2exp::TeX(" $(N,\\lambda)$"),
                          values = c("(1000,1000)" = "#273746", "(400,300)" = "#34495E", "(1000,40)" = "#707B7C", "(150,40)" = "#909497"),
                          labels = c("(1000,1000)" = "(1000,1000)  ", "(400,300)" = "(400,300)  ", "(1000,40)" = "(1000,40)  ", "(150,40)" = "(150,40)  ")) +
      geom_theme
  } else if (white_noise == "mfBm") {
    # comming ...
  }

  return(ggplt)
}

g_autocov_risk_far_mfBm_d1  <- ggpubr::ggarrange(
  ggplot_autocov_risk_by_st(si = 0.2, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_autocov_risk_by_st(si = 0.8, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_autocov_risk_by_st(si = 0.4, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_autocov_risk_by_st(si = 0.7, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
g_autocov_risk_far_mfBm_d1

ggsave(filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/autocov_risk_far_mfBm_d1.png",
       plot = g_autocov_risk_far_mfBm_d1, width = 9, height = 6, units = "in", bg = "white")




