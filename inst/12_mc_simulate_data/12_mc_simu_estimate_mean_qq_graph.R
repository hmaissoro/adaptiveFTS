library(ggplot2)
t0 <- c(0.2, 0.4, 0.7, 0.8)
# Mean function estimates for qqplot ----

estim_mean_qqplot_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0, n_cores = 30){
  # Load mean risk data
  mean_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                                process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design, ".RDS")
  dt_mean_risk <- readRDS(mean_risk_file_name)

  if (white_noise == "mfBm") {
    # Get data file and MC index in the file name
    data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
    id_mc_data <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                       replacement = "", x = data_file_list)
    id_mc_data <- as.numeric(id_mc_data)
    index_mc <- id_mc_data
    # # Get the already done MC repettion
    # Nmc_done <- length(list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda)))
    # if (Nmc_done < length(id_mc_data)) {
    #   data_file_done <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/", design, "_mean_risk/N", N, "lambda", lambda))
    #   id_mc_done <- gsub(pattern = paste0("dt_mean_risk_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", "id_mc=", "|", "_", design,".RDS"),
    #                      replacement = "", x = data_file_done)
    #   id_mc_done <- as.numeric(id_mc_done)
    #   index_mc <- id_mc_done
    # } else {
    #   index_mc <- id_mc_data
    # }

    # Get optimal bandwidth smoothing parameter
    dt_optbw <- dt_mean_risk[, .("optbw" = (h[which.min(mean_risk)]) ** (1.1), "time_taken" = unique(time_taken)), by = c("id_mc", "t")]
    # Estimate mean by mc
    dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_optbw, Ni, lambdai, t0){
      # Load data
      data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                               process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
      dt <- readRDS(data_file_name)

      # Extract data
      dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
      optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
      t0 <- dt_optbw[id_mc == mc_i][order(t), t]

      # Estimate the mean function
      start_time <- Sys.time()
      dt_mean <- estimate_mean(
        data = dt_random_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, optbw = optbw)
      end_time <- Sys.time()
      time_taken <- as.numeric(end_time - start_time)
      time_taken <- time_taken + dt_optbw[id_mc == mc_i, unique(time_taken)]

      # Return and clean
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai,
                                       "time_taken" = time_taken, dt_mean[, .(t, optbw, PN, muhat)])
      dt_res <- data.table::merge.data.table(
        x = dt_res,
        y = unique(dt[ttag == "tcommon", .("t" = tobs, "mutrue" = process_mean)]),
        by = "t")
      rm(optbw, dt_mean) ; gc()
      return(dt_res)
    }, mc.cores = n_cores, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    # Not need for the paper
  }

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_qqplot_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(mean_risk_file_name, file_name, dt_mean_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_mean_qqplot_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Scenario 1
estim_mean_qqplot_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_qqplot_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_qqplot_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_mean_qqplot_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)

# Scenario 3
estim_mean_qqplot_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_qqplot_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_qqplot_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_mean_qqplot_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)


# qqplot graph ----

ggplot_mean_qq_by_t <- function(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_qqplot_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design, ".RDS")
  dt_mean <- readRDS(file_name)
  dt_mean <- unique(dt_mean)

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))

  title_exp <- paste0("t=", ti, ", ", "N=", N, ", $\\lambda$=", lambda)
  dt_mean[, "mu" := sqrt(PN) * (muhat - mutrue), by = t]
  dt_qq <- dt_mean[t == ti]
  dt_qq <- dt_qq[!is.nan(mu)]
  ggplt <-  ggplot(data = dt_qq, aes(sample = mu, group = 1)) +
    # ylim(-18, 18) +
    stat_qq(data = dt_qq, aes(sample = mu, group = 1)) +
    stat_qq_line(data = dt_qq, aes(sample = mu, group = 1)) +
    labs(y = "", x = "N(0,1)") +
    ggtitle(latex2exp::TeX(title_exp)) +
    geom_theme
  return(ggplt)
}

# Scenario 2 :
## t = 0.2
g_mean_qq_far_mfBm_d10.2  <- ggpubr::ggarrange(
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2)

ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_t=", 0.2,"_d2.png"),
       plot = g_mean_qq_far_mfBm_d10.2, width = 9, height = 6, units = "in", bg = "white")
rm(g_mean_qq_far_mfBm_d10.2)

## t = 0.4
g_mean_qq_far_mfBm_d10.4  <- ggpubr::ggarrange(
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.4, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2)

ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_t=", 0.4,"_d2.png"),
       plot = g_mean_qq_far_mfBm_d10.4, width = 9, height = 6, units = "in", bg = "white")
rm(g_mean_qq_far_mfBm_d10.4)

## t = 0.7
g_mean_qq_far_mfBm_d10.7  <- ggpubr::ggarrange(
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.7, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2)

ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_t=", 0.7,"_d2.png"),
       plot = g_mean_qq_far_mfBm_d10.7, width = 9, height = 6, units = "in", bg = "white")
rm(g_mean_qq_far_mfBm_d10.7)

## t = 0.8
g_mean_qq_far_mfBm_d10.8  <- ggpubr::ggarrange(
  ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = 0.8, process = "FAR", white_noise = "mfBm", design = "d1"),
  nrow = 2, ncol = 2)

ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_t=", 0.8,"_d2.png"),
       plot = g_mean_qq_far_mfBm_d10.8, width = 9, height = 6, units = "in", bg = "white")
rm(g_mean_qq_far_mfBm_d10.8)

# Scenario 3 :
for(ti in c(0.2, 0.4, 0.7, 0.8)){
  g_mean_qq_far_mfBm_d3  <- ggpubr::ggarrange(
    ggplot_mean_qq_by_t(N = 150, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_qq_by_t(N = 1000, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_qq_by_t(N = 400, lambda = 300, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_qq_by_t(N = 1000, lambda = 1000, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    nrow = 2, ncol = 2)

  ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_qq_far_mfBm_t=", ti,"_d3.png"),
         plot = g_mean_qq_far_mfBm_d3, width = 9, height = 6, units = "in", bg = "white")
  rm(g_mean_qq_far_mfBm_d3)
}

# Sandardize mean distribution version ----
estim_standardise_mean_qqplot_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", t0, n_cores = 30){
  # Load mean risk data
  mean_risk_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_risk_",
                                process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design, ".RDS")
  dt_mean_risk <- readRDS(mean_risk_file_name)

  # Get data file and MC index in the file name
  data_file_list <- list.files(paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda))
  index_mc <- gsub(pattern = paste0("dt_mc_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_id_mc=", "|", "_", design, ".RDS"),
                     replacement = "", x = data_file_list)
  index_mc <- as.numeric(index_mc)

  # Get optimal bandwidth smoothing parameter
  dt_optbw <- dt_mean_risk[, .("optbw" = (h[which.min(mean_risk)]) ** (1.1), "time_taken" = unique(time_taken)), by = c("id_mc", "t")]
  # Estimate mean by mc
  dt_mean_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt_optbw, Ni, lambdai, t0){
    # Load data
    data_file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/", design, "_slice/N", N, "lambda", lambda,"/dt_mc_",
                             process,"_", white_noise, "_", "N=", Ni, "_lambda=", lambdai, "_", "id_mc=", mc_i, "_", design,".RDS")
    dt <- readRDS(data_file_name)

    # Extract data
    dt_random_mc <- dt[ttag == "trandom"][id_mc == mc_i]
    id_curves <- dt_random_mc[, unique(id_curve)]
    optbw <- dt_optbw[id_mc == mc_i][order(t), optbw]
    t0 <- dt_optbw[id_mc == mc_i][order(t), t]

    # Estimate sigma
    dt_sig <- estimate_sigma(data = dt_random_mc, idcol = "id_curve", tcol = "tobs", ycol = "X", t = t0)
    dt_sig_mutrue <- data.table::merge.data.table(
      x = dt_sig,
      y = unique(dt[ttag == "tcommon"][id_mc == mc_i, .("t" = tobs, "mutrue" = process_mean)]),
      by = "t")

    # Smooth curves with optimal bandwidth parameters
    dt_Xhat <- data.table::rbindlist(lapply(id_curves, function(curve_index, t, data, optbw, smooth_ker){
      # \pi_n(t,h)
      Tn <- data[id_curve == curve_index, tobs]
      Yn <- data[id_curve == curve_index, X]
      pi_n <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw_mean){
        as.numeric(abs(Tn - t[tidx]) <= optbw_mean[tidx])
      }, Tn = Tn, t = t, optbw_mean = optbw)
      pi_n <- t(pi_n)
      pi_n <- as.numeric(rowSums(pi_n, na.rm = TRUE) >= 1)

      # \sum{i=1}^{M_n} W_{n,i}(t,h)^2
      Wn_square <- sapply(X = 1:length(t), function(tidx, Tn, t, optbw){
        K_ni <- smooth_ker((Tn - t[tidx]) / optbw[tidx])
        kn_sum <- sum(K_ni)
        if (kn_sum == 0) {
          W_ni <- K_ni
        } else {
          W_ni <- {K_ni / kn_sum} ** 2
        }
        return(W_ni)
      }, Tn = Tn, t = t, optbw = optbw)
      Wn_square <- t(Wn_square)
      Wn_square <- rowSums(Wn_square)

      # \widehat X(t;h)
      Xhat <- mapply(function(t, h, Yn, Tn, ker){
        estimate_nw(y = Yn, t = Tn, tnew = t, h = h, smooth_ker = ker)$yhat
      }, t = t, h = optbw,
      MoreArgs = list(Yn = Yn, Tn = Tn, ker = smooth_ker))

      dt_res <- data.table::data.table("id_curve" = curve_index, "t" = t, "Xhat" = Xhat, "pi_n" = pi_n, "Wn_square" = Wn_square, "sigma" = dt_sig$sig)
      return(dt_res)
    }, data = dt_random_mc, t = t0, optbw = optbw, smooth_ker = epanechnikov))

    # Estimate mean function \widehat \mu(t)
    dt_Xhat[is.nan(Xhat) & pi_n == 0, Xhat := 0]

    # Add sigma estimates and true mean function
    dt_Xhat <- data.table::merge.data.table(x = dt_Xhat, y = dt_sig_mutrue, by = "t")
    dt_Xhat[, PN := sum(pi_n), by = t]
    dt_Xhat <- dt_Xhat[!is.nan(Xhat)]
    dt_Xhat[, "muhat" := sum(pi_n * Xhat) / PN, by = "t"]
    dt_Xhat[, "Sigma_element" := sum(pi_n * (Wn_square ** 2) * (sigma ** 2)) / PN, by = "t"]
    dt_Xhat[, "SSmu_element" := sum(pi_n * (Xhat - mutrue) / sqrt(PN)), by = "t"]

    # Add the optimal bandwidth
    dt_muhat <- unique(dt_Xhat[, list(t, muhat, mutrue, PN, Sigma_element, SSmu_element)])
    dt_muhat[, "id_mc" := mc_i]
    dt_muhat <- data.table::merge.data.table(x = dt_optbw[, .(id_mc, t, "hN" = optbw)], y = dt_muhat, by = c("id_mc", "t"))

    return(dt_muhat)
  }, mc.cores = n_cores, dt_optbw = dt_optbw, Ni = N, lambdai = lambda))

  # Estimate the variance part : \mathbb S_mu(t)
  dt_mean_mc[, SSmu_t := var(SSmu_element), by = "t"]

  # Estimate the variance part : \Sigma(t)
  dt_mean_mc[, Sigma_t := mean(Sigma_element), by = "t"]

  # Estimate the standarize mean estimation
  dt_mean_mc[, mean_standardise := sqrt(PN) * (muhat - mutrue) / sqrt(Sigma_t + SSmu_t),  by = "t"]

  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_standardise_qqplot_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  saveRDS(object = dt_mean_mc, file = file_name)
  rm(mean_risk_file_name, file_name, dt_mean_mc, dt_optbw) ; gc()

  return(paste0("Done : dt_mean_standardise_qqplot_estimates_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Scenario 1
estim_standardise_mean_qqplot_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d1", n_cores = 30)

# Scenario 3
estim_standardise_mean_qqplot_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)
estim_standardise_mean_qqplot_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", design = "d3", n_cores = 30)

# Q-Q plots
ggplot_mean_standardised_qq_by_t <- function(N = 400, lambda = 300, ti = 0.2, process = "FAR", white_noise = "mfBm", design = "d1"){
  ## Load data
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/mean_estimates/dt_mean_standardise_qqplot_estimates_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design, ".RDS")
  dt_mean <- readRDS(file_name)
  dt_mean <- unique(dt_mean)

  ## ggplot parameters
  geom_theme <- theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width= unit(0.8, 'cm'))

  title_exp <- paste0("t=", ti, ", ", "N=", N, ", $\\lambda$=", lambda)
  dt_qq <- dt_mean[t == ti]
  dt_qq <- dt_qq[!is.nan(mean_standardise)]
  ggplt <-  ggplot(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    # ylim(-18, 18) +
    stat_qq(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    # stat_qq_line(data = dt_qq, aes(sample = mean_standardise, group = 1)) +
    geom_abline() +
    labs(y = "", x = "N(0,1)") +
    ggtitle(latex2exp::TeX(title_exp)) +
    geom_theme
  return(ggplt)
}

# Scenario 1 :
for(ti in c(0.2, 0.4, 0.7, 0.8)){
  g_mean_qq_sd_far_mfBm_d2  <- ggpubr::ggarrange(
    ggplot_mean_standardised_qq_by_t(N = 150, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d1"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d1"),
    ggplot_mean_standardised_qq_by_t(N = 400, lambda = 300, ti = ti, process = "FAR", white_noise = "mfBm", design = "d1"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 1000, ti = ti, process = "FAR", white_noise = "mfBm", design = "d1"),
    nrow = 2, ncol = 2)

  ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_standardised_qq_far_mfBm_t=", ti,"_d2.png"),
         plot = g_mean_qq_sd_far_mfBm_d2, width = 9, height = 7, units = "in", bg = "white")
  rm(g_mean_qq_sd_far_mfBm_d2)
}

# Scenario 2 :
for(ti in c(0.2, 0.4, 0.7, 0.8)){
  g_mean_qq_sd_far_mfBm_d3  <- ggpubr::ggarrange(
    ggplot_mean_standardised_qq_by_t(N = 150, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 40, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_standardised_qq_by_t(N = 400, lambda = 300, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    ggplot_mean_standardised_qq_by_t(N = 1000, lambda = 1000, ti = ti, process = "FAR", white_noise = "mfBm", design = "d3"),
    nrow = 2, ncol = 2)

  ggsave(filename = paste0("./inst/12_mc_simulate_data/graphs/paper_graphs/mean_standardised_qq_far_mfBm_t=", ti,"_d3.png"),
         plot = g_mean_qq_sd_far_mfBm_d3, width = 9, height = 7, units = "in", bg = "white")
  rm(g_mean_qq_sd_far_mfBm_d3)
}
