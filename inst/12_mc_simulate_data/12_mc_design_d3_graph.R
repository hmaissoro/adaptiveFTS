library(data.table)
library(ggplot2)

# Simulation global parameters----
sig <- 0.25
mc <- 100
t0 <- c(0.2, 0.4, 0.7, 0.8)

# (N, \lambda) \in {(400, 300), (1000, 1000), (150, 40), (1000, 40)}
# (N, \lambda) \in {(150, 40), (1000, 40)}
# (N, \lambda) = (150, 1000)

# Note that (N, \lambda) = (1000, 1000) for design = "d1" not computed

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}


## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}

ker_d3 <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.2)
dt_gen <- simulate_far(N = 50, lambda = 400,
                       tdesign = "common",
                       Mdistribution = NULL,
                       tdistribution = NULL,
                       tcommon = seq(0.01, 0.99, len = 100),
                       hurst_fun = Hlogistic,
                       L = 4,
                       far_kernel = ker_d3,
                       far_mean = get_real_data_mean,
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)

saveRDS(dt_gen, "./inst/12_mc_simulate_data/FAR/data/dt_common_d3.RDS")


## All curves
dt_gen <- readRDS("./inst/12_mc_simulate_data/FAR/data/dt_common_d3.RDS")
dt <- copy(dt_gen)
dt[, id_curve := as.factor(id_curve)]
dt <- dt[, .(id_curve, "t" = tobs, X)]
g_d3 <- ggplot(dt, aes(x = t, y = X, color = id_curve, group = id_curve)) +
  geom_line() +
  ylim(233, 246) +
  scale_color_grey() +
  geom_abline(intercept = 0.7, linetype = 2, size = 0.9) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
g_d3
# Save and clean
ggsave(plot = g_real_data_all_curves, filename = file.path(figures_path, "real_data_all_curves.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
rm(g_real_data_all_curves) ; gc()

g_d3 <- ggplot(dt, aes(x = t, y = X, color = id_curve, group = id_curve)) +
  geom_line(size = 0.9)  +
  ylim(233, 247) +
  geom_vline(xintercept = 0.7, linetype = 2, size = 0.9) +
  geom_text(x = 0.7, label = "t=0.7", y = 234, angle = 90, vjust = -0.2, size = 7) +
  xlab(label = "t") +
  ylab(latex2exp::TeX("$X(t)$  ")) +
  scale_color_grey(start = 0.2, end = 0.5) +
  theme_minimal() +
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_d3
ggsave(plot = g_d3, filename = file.path("./inst/12_mc_simulate_data/graphs/paper_graphs/real_curves_d3.png"),
       width = 12, height = 5, units = "in", dpi = 300, bg = "white")



