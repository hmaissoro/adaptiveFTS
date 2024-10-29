library(ggplot2)
source("./R/03_estimate_regularity.R")
Rcpp::sourceCpp("./src/08_estimate_curve_cpp.cpp")

# Import the data
dt <- readRDS("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=1_fts_model_2.RDS")
dt_common <- dt[ttag == "tcommon"]
dt <- dt[ttag == "trandom"]
data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")

# Estimation parameters
t0 <- c(0.2, 0.4, 0.7, 0.8)
path <- "../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/"
file_list <- list.files(path)
sort(as.numeric(gsub("dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=|_fts_model_2.RDS", "" ,file_list)))

dt_res_blup <- data.table::rbindlist(parallel::mclapply(1:100, function(id){
  dt <- readRDS(paste0("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N150lambda40/dt_mc_FAR_mfBm_N=150_lambda=40_id_mc=", id, "_fts_model_2.RDS"))
  dt_common <- dt[ttag == "tcommon"]
  dt <- dt[ttag == "trandom"]
  dt[! id_curve == 150 & tobs > 0.5]
  data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
  data_to_use <- data_prepared[-which(id_curve == 150 & tobs > 0.5)]

  res_blup_one <- estimate_curve(data = data_to_use, id_curve = 150, t = t0, kernel_name = "epanechnikov")

  dt_res <- data.table::as.data.table(res_blup_one$res_blup)
  names(dt_res) <- c("t", "muhat", "blup")
  dt_res <- data.table::merge.data.table(x = dt_res, y = dt_common[id_curve == 150, .("t" = tobs, "Xtrue" = X)], by = "t")
  dt_res[, id_mc := id]
  return(dt_res)
}, mc.cores = 20))

saveRDS(dt_res_blup, "../../curve_reconstruction/adaptive_blup_emp_study/blup/dt_res_blup_N=150_lambda=40_fts_model_2.RDS")

dt_res_blup[, error := Xtrue - blup]
dt_res_blup[, t := as.factor(t)]
ggplt <- ggplot(data = dt_res_blup[error < 10 & error > -10], mapping = aes(x = t, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red") +
  ggtitle(latex2exp::TeX("N = 150, $\\lambda$ = 40")) +
  ylab(latex2exp::TeX("$X_{n_0}(t) - \\hat{X}_{n_0}(t)$")) +
  xlab("t") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.width= unit(0.8, 'cm'),
        legend.position = "top")
ggsave(plot = ggplt, filename = file.path("../../curve_reconstruction/adaptive_blup_emp_study/blup/blup_N=150_lambda=40_fts_model_2.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")







res_blup_one <- estimate_curve(data = data_prepared, id_curve = 150, t = t0, kernel_name = "epanechnikov")
res_blup_two <- estimate_curve(data = data_prepared, id_curve = 150, t = t0, kernel_name = "epanechnikov")

res_blup$Cov_n0 %>% View()
res_blup$Sigma_n0
res_blup$cov_n0_lag %>% View()
res_blup$autocov_n0_lag__n0
res_blup$autocov_n0_lag_tvec
res_blup$res_blup
dt_common[id_curve == 150]

res_blup$res_blup
res_blup_one$res_blup
res_blup_two$res_blup


# N = 400 lambda = 300 ----
dt_res_blup_N400lambda300 <- data.table::rbindlist(parallel::mclapply(1:20, function(id){
  dt <- readRDS(paste0("../../curve_reconstruction/adaptive_blup_emp_study/data/fts_model_2/N400lambda300/dt_mc_FAR_mfBm_N=400_lambda=300_id_mc=", id, "_fts_model_2.RDS"))
  dt_common <- dt[ttag == "tcommon"]
  dt <- dt[ttag == "trandom"]
  # dt[! id_curve == 400 & tobs > 0.5]
  data_prepared <- format_data(data = dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
  data_to_use <- data_prepared[-which(id_curve == 400 & tobs > 0.5)]

  res_blup_one <- estimate_curve(data = data_to_use, id_curve = 400, t = t0, kernel_name = "epanechnikov")

  dt_res <- data.table::as.data.table(res_blup_one$res_blup)
  names(dt_res) <- c("t", "muhat", "blup")
  dt_res <- data.table::merge.data.table(x = dt_res, y = dt_common[id_curve == 400, .("t" = tobs, "Xtrue" = X)], by = "t")
  dt_res[, id_mc := id]
  return(dt_res)
}, mc.cores = 20))

saveRDS(dt_res_blup_N400lambda300, "../../curve_reconstruction/adaptive_blup_emp_study/blup/dt_res_blup_N=400_lambda=300_fts_model_2.RDS")

dt_res_blup[, error := Xtrue - blup]
dt_res_blup[, t := as.factor(t)]
ggplt <- ggplot(data = dt_res_blup[error < 10 & error > -10], mapping = aes(x = t, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red") +
  ggtitle(latex2exp::TeX("N = 150, $\\lambda$ = 40")) +
  ylab(latex2exp::TeX("$X_{n_0}(t) - \\hat{X}_{n_0}(t)$")) +
  xlab("t") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(size = 11, margin = margin(t = 10, r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.width= unit(0.8, 'cm'),
        legend.position = "top")
ggsave(plot = ggplt, filename = file.path("../../curve_reconstruction/adaptive_blup_emp_study/blup/blup_N=150_lambda=40_fts_model_2.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")



