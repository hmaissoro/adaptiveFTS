library(data.table)

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

## Logistic constant hurst function
Hvec <- c(0.4, 0.5, 0.7)

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

# zero-mean function
zero_mean_func <- function(t) 0 * t


# Local regularity estimation new function
estimate_locreg_test <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                                 t = 1/2, Delta = NULL, h = NULL,
                                 smooth_ker = epanechnikov, center = TRUE){
  # Control easy checkable arguments
  if (! (methods::is(t, "numeric") & all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.")
  if (! methods::is(smooth_ker, "function"))
    stop("'smooth_ker' must be a function.")
  if (! methods::is(center, "logical"))
    stop("'center' must be a TRUE or FALSE.")

  # Control and format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  N <- data[, length(unique(id_curve))]

  # Control and set arguments depending on the data
  if (! is.null(Delta)) {
    if (! (methods::is(Delta, "numeric") & data.table::between(Delta, 0, 1) & length(Delta) == 1))
      stop("'Delta' must be a numeric scalar value between 0 and 1.")
  } else {
    lambdahat <- mean(data[, .N, by = "id_curve"][, N])
    Delta <- 2 * exp(-log(lambdahat) ** 0.72)
  }

  # Control on the pre-smoothing bandwidth
  if (! is.null(h)) {
    # h is a vector or a scalar
    if (! all(methods::is(h, "numeric") & data.table::between(h, 0, 1))){
      stop("'h' must be a numeric vector or scalar value(s) between 0 and 1.")
    } else if (length(h) > 1 & length(h) != N) {
      stop("If 'h' is given as a vector, its length must be equal to the number of curves in 'data'.")
    }
  } else {
    # If h = NULL, choose the bandwidth by CV
    if (N > 50) {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = 30, smooth_ker = smooth_ker)
      h <- dt_optbw[, median(optbw)]
      rm(dt_optbw) ; gc()
    } else {
      dt_optbw <- get_nw_optimal_bw(
        data = data, idcol = "id_curve", tcol = "tobs", ycol = "X",
        bw_grid = NULL, nsubset = NULL, smooth_ker = smooth_ker)
      h <- dt_optbw[, optbw]
      rm(dt_optbw) ; gc()
    }
  }

  # If the bandwidth is given as scalar or computed
  if (length(h) == 1) h <- rep(h, N)

  # Step 1: estimate the curve at
  ## Take into account the cases where t-Delta/2 < 0 and t + Delta/2
  dt_t <- data.table::rbindlist(lapply(t, function(ti, D){
    if ((ti - D / 2 ) <= 0) {
      t1 <- ti
      t2 <- ti + D / 2
      t3 <- ti + D
    } else if ((ti + D / 2 ) >= 1) {
      t3 <- ti
      t2 <- ti - D / 2
      t1 <- ti - D
    } else {
      t1 <- ti - Delta / 2
      t2 <- ti
      t3 <- ti + Delta / 2
    }
    return(data.table::data.table("t1" = t1, "t2" = t2, "t3" = t3))
  }, D = Delta))
  t1 <- dt_t[, t1]
  t2 <- dt_t[, t2]
  t3 <- dt_t[, t3]
  rm(dt_t)

  dt_smooth <- data.table::rbindlist(lapply(data[, unique(id_curve)], function(i, data, h, t1, t2, t3){
    ## smooth an estimate
    dt1 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t1, h = h[i],
                       smooth_ker = smooth_ker)
    dt2 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t2, h = h[i],
                       smooth_ker = smooth_ker)
    dt3 <- estimate_nw(y = data[id_curve == i, X],
                       t = data[id_curve == i, tobs],
                       tnew = t3, h = h[i],
                       smooth_ker = smooth_ker)

    dt_out <- data.table::data.table(
      id_curve = i, t1 = t1, xt1 = dt1$yhat,
      t2 = t2, xt2 = dt2$yhat,
      t3 = t3, xt3 = dt3$yhat
    )
    # rm(xtilde, inKernelSupp)

    return(dt_out)
  }, h = h, data = data, t1 = t1, t2 = t2, t3 = t3))

  # Step 2 : estimate regularity parameters

  dt_reg <- data.table::rbindlist(lapply(1:length(t2), function(i, dt_smooth, t1, t2, t3, t, center){
    ## Extract X_1(g),...,X_N(g) where g = t1, t2 or t3
    xt1 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt1]
    xt2 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt2]
    xt3 <- dt_smooth[t1 == t1[i] & t2 == t2[i] & t3 == t3[i], xt3]

    ## Remove NaN values
    any_nan <- (is.nan(xt1) | is.nan(xt2) | is.nan(xt3))
    xt1 <- xt1[! any_nan]
    xt2 <- xt2[! any_nan]
    xt3 <- xt3[! any_nan]
    rm(any_nan)

    ## Remove extreme values
    rxt1 <- (xt1 >= quantile(xt1, 0.025, na.rm = TRUE)) & (xt1 <= quantile(xt1, 0.975, na.rm = TRUE))
    rxt2 <- (xt2 >= quantile(xt2, 0.025, na.rm = TRUE)) & (xt2 <= quantile(xt2, 0.975, na.rm = TRUE))
    rxt3 <- (xt3 >= quantile(xt3, 0.025, na.rm = TRUE)) & (xt3 <= quantile(xt3, 0.975, na.rm = TRUE))
    # xt1 <- xt1[rxt1 & rxt2 & rxt3]
    # xt2 <- xt2[rxt1 & rxt2 & rxt3]
    # xt3 <- xt3[rxt1 & rxt2 & rxt3]
    # Nused <- sum(rxt1 & rxt2 & rxt3)
    # rm(rxt1, rxt2, rxt3)

    ## Center data if the argument center = TRUE
    # xt1 <- xt1 - center * mean(xt1)
    # xt2 <- xt2 - center * mean(xt2)
    # xt3 <- xt3 - center * mean(xt3)

    ## Compute Unweighed local regularity parameters
    # theta_t1_t3 <- mean((xt1 - xt3) ** 2, na.rm = TRUE)
    # theta_t1_t2 <- mean((xt1 - xt2) ** 2, na.rm = TRUE)
    # theta_t2_t3 <- mean((xt2 - xt3) ** 2, na.rm = TRUE)
    theta_t1_t3 <- sum((xt1[rxt1 & rxt3] - xt3[rxt1 & rxt3]) ** 2, na.rm = TRUE)
    theta_t1_t2 <- sum((xt1 - xt2) ** 2, na.rm = TRUE)
    theta_t2_t3 <- sum((xt2[rxt2 & rxt3] - xt3[rxt2 & rxt3]) ** 2, na.rm = TRUE)
    Ht <- (log(theta_t1_t3) - log(theta_t2_t3))  / (2 * log(2))
    Lt <- mean((xt2[rxt2 & rxt3] - xt3[rxt2 & rxt3]) ** 2, na.rm = TRUE) / (abs(t2[i] - t3[i])**(2 * Ht))

    ## Return the result
    # dt_out <- data.table(t = t[i], Ht = Ht, Lt = Lt, Nused = Nused)
    dt_out <- data.table(t = t[i], Ht = Ht, Lt = Lt)
    rm(theta_t1_t3, theta_t2_t3, theta_t1_t2, Ht, Lt, xt1, xt2, xt3)

    return(dt_out)
  }, dt_smooth = dt_smooth, t1 = t1, t2 = t2, t3 = t3, t = t, center = center))

  dt_reg[, c("locreg_bw", "Delta") := list(median(h), Delta)]

  # data.table::setcolorder(x = dt_reg, neworder = c("t", "Delta", "Nused", "locreg_bw", "Ht", "Lt"))
  data.table::setcolorder(x = dt_reg, neworder = c("t", "Delta", "locreg_bw", "Ht", "Lt"))

  return(dt_reg)
}

# local regularity estimation ----
estim_locreg_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1"){
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  dt <- dt[ttag == "trandom"][, .SD, .SDcols = ! c("ttag", "process_mean")]
  index_mc <- dt[, unique(id_mc)]

  if (white_noise == "mfBm") {
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      dt_mc <- dt_mc[order(id_curve, tobs)]
      bw <- unique(dt_mc[, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
      bw_plus_mean <- unique(dt_mc[, .(id_curve, presmooth_bw_plus_mean)])[order(id_curve), presmooth_bw_plus_mean]

      # Estimate the local regularity
      ## Delta
      lambdahat <- mean(dt_mc[, .N, by = id_curve][, N])
      delta <- exp(- log(lambdahat) ** 0.35)

      ## Centered process
      dt_locreg <- estimate_locreg_test(
        data = dt_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X",
        t = t0, Delta = delta, h = bw,
        smooth_ker = epanechnikov, center = FALSE)

      ## process plus mean function
      dt_locreg_plus_mean <- estimate_locreg_test(
        data = dt_mc, idcol = "id_curve",
        tcol = "tobs", ycol = "X_plus_mean",
        t = t0, Delta = delta, h = bw_plus_mean,
        smooth_ker = epanechnikov, center = FALSE)
      data.table::setnames(x = dt_locreg_plus_mean,
                           old = c("locreg_bw", "Delta", "Ht", "Lt"),
                           new = c("locreg_bw_plus_mean", "Delta_plus_mean", "Ht_plus_mean", "Lt_plus_mean"))

      ## Merge local regularity parameters estimates and return the obtained result
      dt_locreg_res <- data.table::merge.data.table(x = dt_locreg, y = dt_locreg_plus_mean, by = "t")
      dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg_res)
      rm(dt_locreg_res, dt_locreg, dt_locreg_plus_mean) ; gc()
      return(dt_res)
    }, mc.cores = 75, dt = dt, Ni = N, lambdai = lambda))

  } else if (white_noise == "fBm") {
    # Estimate local regularity by mc
    dt_reg_mc <- data.table::rbindlist(parallel::mclapply(index_mc, function(mc_i, dt, Ni, lambdai, t0){
      # Extract and sort data
      dt_mc <- dt[id_mc == mc_i]
      dt_mc <- dt_mc[order(id_curve, tobs)]
      Hvec <- dt_mc[, sort(unique(Htrue))]

      dt_by_Hvec <- data.table::rbindlist(lapply(Hvec, function(Hi, dt_mc, t0){
        ## Set Delta
        lambdahat <- mean(dt_mc[Htrue == Hi, .N, by = id_curve][, N])
        delta <- exp(- log(lambdahat) ** 0.25)

        ## Extract bandwidth
        bw <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw)])[order(id_curve), presmooth_bw]
        bw_plus_mean <- unique(dt_mc[Htrue == Hi, .(id_curve, presmooth_bw_plus_mean)])[order(id_curve), presmooth_bw_plus_mean]

        ## Centered process
        dt_locreg <- estimate_locreg_test(
          data = dt_mc[Htrue == Hi], idcol = "id_curve",
          tcol = "tobs", ycol = "X",
          t = t0, Delta = delta, h = bw,
          smooth_ker = epanechnikov, center = FALSE)

        ## process plus mean function
        dt_locreg_plus_mean <- estimate_locreg_test(
          data = dt_mc[Htrue == Hi], idcol = "id_curve",
          tcol = "tobs", ycol = "X_plus_mean",
          t = t0, Delta = delta, h = bw_plus_mean,
          smooth_ker = epanechnikov, center = FALSE)
        data.table::setnames(x = dt_locreg_plus_mean,
                             old = c("locreg_bw", "Delta", "Ht", "Lt"),
                             new = c("locreg_bw_plus_mean", "Delta_plus_mean", "Ht_plus_mean", "Lt_plus_mean"))

        ## Merge local regularity parameters estimates and return the obtained result
        dt_locreg_res <- data.table::merge.data.table(x = dt_locreg, y = dt_locreg_plus_mean, by = "t")
        dt_res <- data.table::data.table("id_mc" = mc_i, "N" = Ni, "lambda" = lambdai, "lambdahat" = lambdahat, dt_locreg_res)
        dt_res[, Htrue := Hi]

        ## Clean
        rm(dt_locreg_res, dt_locreg, dt_locreg_plus_mean) ; gc()
        return(dt_res)
      }, dt_mc = dt_mc, t0 = t0))

      return(dt_by_Hvec)
    }, mc.cores = 75, dt = dt, Ni = N, lambdai = lambda, t0 = t0))
  }
  ## Save
  file_name <- paste0("./inst/12_mc_simulate_data/", process, "/locreg_estimates/dt_locreg_",
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_test.RDS")
  saveRDS(object = dt_reg_mc, file = file_name)
  rm(dt, file_title, file_name, dt_reg_mc) ; gc() ; gc()

  return(paste0("Done : dt_locreg_", process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS at ", Sys.time()))
}

# Estimate local regularity ----
## design 1 ----
### FAR ----
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1")

### FMA ----
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1")
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1")

## design 2 ----
### FAR ----
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2")

### FMA ----
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2")
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2")

## design 3 ----
### FAR ----
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3")
estim_locreg_fun(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3")

### FMA ----
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3")
estim_locreg_fun(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3")

# Plot local regularity estimates ----
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
                      process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_test.RDS")
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
## design 1 ----
### FAR ----
g_locreg_far_mfBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_mfBm_d1_test.png", plot = g_locreg_far_mfBm_d1,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_far_fBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_fBm_d1_test.png", plot = g_locreg_far_fBm_d1,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

### FMA ----
g_locreg_fma_mfBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_mfBm_d1_test.png", plot = g_locreg_fma_mfBm_d1,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_fma_fBm_d1 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d1", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_fBm_d1_test.png", plot = g_locreg_fma_fBm_d1,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

## design 2 ----
### FAR ----
g_locreg_far_mfBm_d2 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d2", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_mfBm_d2_test.png", plot = g_locreg_far_mfBm_d2,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_far_fBm_d2 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d2", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_fBm_d2_test.png", plot = g_locreg_far_fBm_d2,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

### FMA ----
g_locreg_fma_mfBm_d2 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d2", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_mfBm_d2_test.png", plot = g_locreg_fma_mfBm_d2,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_fma_fBm_d2 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d2", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_fBm_d2_test.png", plot = g_locreg_fma_fBm_d2,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

## design 3 ----
### FAR ----
g_locreg_far_mfBm_d3 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_mfBm_d3_test.png", plot = g_locreg_far_mfBm_d3,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_far_fBm_d3 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FAR", white_noise = "fBm", design = "d3", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_far_fBm_d3_test.png", plot = g_locreg_far_fBm_d3,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

### FMA ----
g_locreg_fma_mfBm_d3 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "mfBm", design = "d3", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_mfBm_d3_test.png", plot = g_locreg_fma_mfBm_d3,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

g_locreg_fma_fBm_d3 <- gridExtra::grid.arrange(
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Ht"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Ht_plus_mean"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Lt"),
  ggplot_locreg(N = 400, lambda = 300, process = "FMA", white_noise = "fBm", design = "d3", param = "Lt_plus_mean"),
  nrow = 2, ncol = 2
)
ggsave(filename = "./inst/12_mc_simulate_data/graphs/locreg_fma_fBm_d3_test.png", plot = g_locreg_fma_fBm_d3,
       width = 9.98, height = 8.5, units = "in", dpi = 300)

