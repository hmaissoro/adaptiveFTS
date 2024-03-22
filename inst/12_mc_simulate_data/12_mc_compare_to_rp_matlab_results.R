library(ggplot2)

# N = 150, lambda = 40

graph_mean_compare <- function(N = 150, lambda = 40, design = "d1", version = "", element = "mean"){
  # Load and prepare data
  dt_mean <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N, "_lambda=", lambda, "_", design, version, ".RDS"))
  # dt_gammatilde <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/autocov_estimates/dt_autocovtilde_FAR_mfBm_N=5000_", design, ".RDS"))
  dt_mean_rp <- data.table::fread(paste0("../../matlab/rubin_mean_autocov_method/N", N, "lambda", lambda, design, ".csv"))

  # Add mean function "true" estimates
  # dt_mean <- data.table::merge.data.table(
  #   x = dt_mean,
  #   y = unique(dt_gammatilde[, .(t, "mutilde" = mutilde_t)]),
  #   by = "t")

  dt_mean_rp_melt <- data.table::melt(
    data = dt_mean_rp[, .SD, .SDcols = ! "bdth"], id.vars = "id_mc",
    value.vars = c("ft2", "ft4", "ft7", "ft8"), value.name = "muhat_RP")
  dt_mean_rp_melt[, t := as.numeric(gsub("ft", "", variable)) / 10]
  dt_mean_rp_melt[, variable := NULL]
  dt_merge <- data.table::merge.data.table(x = dt_mean, y = dt_mean_rp_melt, by = c("id_mc", "t"))

  # Theme global parameters
  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
          axis.text.x =  element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.width= unit(0.8, 'cm'),
          legend.position = "top")

  # Calculate the error
  dt_merge[, error_ratio := log(abs(muhat - mutrue) / abs(muhat_RP - mutrue))]
  dt_merge[, t := as.factor(t)]

  if (element == "mean") {
    if (design == "d1") {
      y_lim <- c(-6, 6)
    } else {
      y_lim <- c(-6, 6)
    }
    g <- ggplot(data = dt_merge, mapping = aes(x = t, y = error_ratio)) +
      geom_boxplot() +
      ylab("") +
      ylim(y_lim) +
      ggtitle(latex2exp::TeX(paste0("$N=", N , ", $\\lambda$=", lambda))) +
      geom_hline(yintercept = 0, linetype = 2, size = 0.8) +
      geom_theme
  } else if (element == "bw") {
    dt_optbw <- rbind(dt_mean[, .(id_mc, "t" = paste("t=", t), optbw)],
                      dt_mean_rp[, .(id_mc, "t" = "RP", "optbw" = bdth)])
    dt_optbw <- dt_optbw[order(t, decreasing = FALSE)]
    if (design == "d1") {
      y_lim <- c(0, 0.12)
    } else {
      y_lim <- c(0, 0.12)
    }

    g <- ggplot(data = dt_optbw, mapping = aes(x = t, y = optbw)) +
      geom_boxplot() +
      ylim(y_lim) +
      ylab("") +
      xlab("") +
      ggtitle(latex2exp::TeX(paste0("$N=", N , ", $\\lambda$=", lambda))) +
      geom_theme
  } else if (element == "sd") {
    N_vec <- c(150, 1000, 400, 1000)
    lambda_vec <- c(40, 40, 300, 1000)

    dt_res <- data.table::rbindlist(lapply(1:4, function(i){
      # Load and prepare data
      dt_mean <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design, ".RDS"))
      dt_mean_rp <- data.table::fread(paste0("../../matlab/rubin_mean_autocov_method/N", N_vec[i], "lambda", lambda_vec[i], design, ".csv"))

      dt_mean_rp_melt <- data.table::melt(
        data = dt_mean_rp[, .SD, .SDcols = ! "bdth"], id.vars = "id_mc",
        value.vars = c("ft2", "ft4", "ft7", "ft8"), value.name = "muhat_RP")
      dt_mean_rp_melt[, t := as.numeric(gsub("ft", "", variable)) / 10]
      dt_mean_rp_melt[, variable := NULL]
      dt_merge <- data.table::merge.data.table(x = dt_mean, y = dt_mean_rp_melt, by = c("id_mc", "t"))
      dt_merge[, .("SD_ratio" =  round(sd(muhat) / sd(muhat_RP), 4)), by = c("N", "lambda", "t")]
    }))
    dt_res <- data.table::dcast(dt_res, formula = N + lambda ~ t, value.var = "SD_ratio")
    g <- dt_res[order(lambda)]
  } else if (element == "bias") {
    N_vec <- c(150, 1000, 400, 1000)
    lambda_vec <- c(40, 40, 300, 1000)

    dt_res <- data.table::rbindlist(lapply(1:4, function(i){
      # Load and prepare data
      dt_mean <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design, ".RDS"))
      dt_mean_rp <- data.table::fread(paste0("../../matlab/rubin_mean_autocov_method/N", N_vec[i], "lambda", lambda_vec[i], design, ".csv"))

      dt_mean_rp_melt <- data.table::melt(
        data = dt_mean_rp[, .SD, .SDcols = ! "bdth"], id.vars = "id_mc",
        value.vars = c("ft2", "ft4", "ft7", "ft8"), value.name = "muhat_RP")
      dt_mean_rp_melt[, t := as.numeric(gsub("ft", "", variable)) / 10]
      dt_mean_rp_melt[, variable := NULL]
      dt_merge <- data.table::merge.data.table(x = dt_mean, y = dt_mean_rp_melt, by = c("id_mc", "t"))
      dt_merge[, .("BIAS_ratio" =  round(abs(mean(muhat - mutrue)) - abs(mean(muhat_RP - mutrue)), 4)), by = c("N", "lambda", "t")]
    }))
    dt_res <- data.table::dcast(dt_res, formula = N + lambda ~ t, value.var = "BIAS_ratio")
    g <- dt_res[order(lambda)]
  } else  if(element == "mse") {
    N_vec <- c(150, 1000, 400, 1000)
    lambda_vec <- c(40, 40, 300, 1000)

    dt_res <- data.table::rbindlist(lapply(1:4, function(i){
      # Load and prepare data
      dt_mean <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design, ".RDS"))
      dt_mean_rp <- data.table::fread(paste0("../../matlab/rubin_mean_autocov_method/N", N_vec[i], "lambda", lambda_vec[i], design, ".csv"))

      dt_mean_rp_melt <- data.table::melt(
        data = dt_mean_rp[, .SD, .SDcols = ! "bdth"], id.vars = "id_mc",
        value.vars = c("ft2", "ft4", "ft7", "ft8"), value.name = "muhat_RP")
      dt_mean_rp_melt[, t := as.numeric(gsub("ft", "", variable)) / 10]
      dt_mean_rp_melt[, variable := NULL]
      dt_merge <- data.table::merge.data.table(x = dt_mean, y = dt_mean_rp_melt, by = c("id_mc", "t"))
      # dt_merge[, .("MSE_ratio" =  round(mean((muhat - mutrue)**2) / mean((muhat_RP - mutrue)**2), 4)), by = c("N", "lambda", "t")]
      dt_merge[, .("MSE_ratio" =  round(median((muhat - mutrue)**2 / (muhat_RP - mutrue)**2), 4)), by = c("N", "lambda", "t")]
    }))
    dt_res <- data.table::dcast(dt_res, formula = N + lambda ~ t, value.var = "MSE_ratio")
    g <- dt_res[order(lambda)]
  } else if (element == "mae") {
    N_vec <- c(150, 1000, 400, 1000)
    lambda_vec <- c(40, 40, 300, 1000)

    dt_res <- data.table::rbindlist(lapply(1:4, function(i){
      # Load and prepare data
      dt_mean <- readRDS(file = paste0("./inst/12_mc_simulate_data/FAR/mean_estimates/dt_mean_estimates_FAR_mfBm_N=", N_vec[i], "_lambda=", lambda_vec[i], "_", design, ".RDS"))
      dt_mean_rp <- data.table::fread(paste0("../../matlab/rubin_mean_autocov_method/N", N_vec[i], "lambda", lambda_vec[i], design, ".csv"))

      dt_mean_rp_melt <- data.table::melt(
        data = dt_mean_rp[, .SD, .SDcols = ! "bdth"], id.vars = "id_mc",
        value.vars = c("ft2", "ft4", "ft7", "ft8"), value.name = "muhat_RP")
      dt_mean_rp_melt[, t := as.numeric(gsub("ft", "", variable)) / 10]
      dt_mean_rp_melt[, variable := NULL]
      dt_merge <- data.table::merge.data.table(x = dt_mean, y = dt_mean_rp_melt, by = c("id_mc", "t"))
      dt_merge[, .("MSE_ratio" =  round(mean(abs(muhat - mutrue)) / mean(abs(muhat_RP - mutrue)), 4)), by = c("N", "lambda", "t")]
    }))
    dt_res <- data.table::dcast(dt_res, formula = N + lambda ~ t, value.var = "MSE_ratio")
    g <- dt_res[order(lambda)]
  }

  return(g)
}

g_mean_compare_rp_far_mfBm_d2  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d1", element = "mean"),
  graph_mean_compare(N = 400, lambda = 300, design = "d1", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d1", element = "mean"),
  nrow = 2, ncol = 2)
g_mean_compare_rp_far_mfBm_d2
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_compare_rp_far_mfBm_d2.png", plot = g_mean_compare_rp_far_mfBm_d2,
  width = 9, height = 6, units = "in", bg = "white")


g_mean_compare_rp_far_mfBm_d3  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d3", element = "mean"),
  graph_mean_compare(N = 400, lambda = 300, design = "d3", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d3", element = "mean"),
  nrow = 2, ncol = 2)
g_mean_compare_rp_far_mfBm_d3
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_compare_rp_far_mfBm_d3.png", plot = g_mean_compare_rp_far_mfBm_d3,
  width = 9, height = 6, units = "in", bg = "white")

g_mean_compare_rp_far_mfBm_d5  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d5", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d5", element = "mean"),
  graph_mean_compare(N = 400, lambda = 300, design = "d5", element = "mean"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d5", element = "mean"),
  nrow = 2, ncol = 2)
g_mean_compare_rp_far_mfBm_d5
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_compare_rp_far_mfBm_d5.png",
  plot = g_mean_compare_rp_far_mfBm_d5,
  width = 9, height = 6, units = "in", bg = "white")


# Bandwidth comparison
g_bw_compare_rp_far_mfBm_d2  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d1", element = "bw"),
  graph_mean_compare(N = 400, lambda = 300, design = "d1", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d1", element = "bw"),
  nrow = 2, ncol = 2)
g_bw_compare_rp_far_mfBm_d2
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_bw_compare_rp_far_mfBm_d2.png", plot = g_bw_compare_rp_far_mfBm_d2,
  width = 9, height = 6, units = "in", bg = "white")

g_bw_compare_rp_far_mfBm_d3  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d3", element = "bw"),
  graph_mean_compare(N = 400, lambda = 300, design = "d3", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d3", element = "bw"),
  nrow = 2, ncol = 2)
g_bw_compare_rp_far_mfBm_d3
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_bw_compare_rp_far_mfBm_d3.png", plot = g_bw_compare_rp_far_mfBm_d3,
  width = 9, height = 6, units = "in", bg = "white")


g_bw_compare_rp_far_mfBm_d5  <- ggpubr::ggarrange(
  graph_mean_compare(N = 150, lambda = 40, design = "d5", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 40, design = "d5", element = "bw"),
  graph_mean_compare(N = 400, lambda = 300, design = "d5", element = "bw"),
  graph_mean_compare(N = 1000, lambda = 1000, design = "d5", element = "bw"),
  nrow = 2, ncol = 2)
g_bw_compare_rp_far_mfBm_d5
ggsave(
  filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/mean_bw_compare_rp_far_mfBm_d5.png", plot = g_bw_compare_rp_far_mfBm_d5,
  width = 9, height = 6, units = "in", bg = "white")

# Standard deviation comparaison ----
graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "sd")
graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "sd")

# Bias comparaison ----
graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "bias")
graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "bias")

# MSE comparaison ----
graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "mse")
graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "mse")
graph_mean_compare(N = 150, lambda = 40, design = "d5", element = "mse")

# MAE comparaison ----
graph_mean_compare(N = 150, lambda = 40, design = "d1", element = "mae")
graph_mean_compare(N = 150, lambda = 40, design = "d3", element = "mae")
