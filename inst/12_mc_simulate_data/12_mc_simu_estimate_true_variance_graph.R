library(ggplot2)

# design 1 : var ----
dt <- readRDS("./inst/12_mc_simulate_data/FAR/data/dt_mc_common_for_variance_FAR_mfBm_N=5000_d1.RDS")
dt_var <- dt[, .("var_by_mc" = var(X)), by = c("id_mc", "tobs")]
dt_var <- dt_var[, .("var" = mean(var_by_mc)), by = "tobs"]

## Process variance
g_var_d2 <- ggplot(data = dt_var, aes(x = tobs, y = var)) +
  geom_line(size = 0.5)  +
  ylim(0.5, 5.5) +
  xlab(label = "t") +
  ylab("") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_var_d2
ggsave(plot = g_var_d2, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/process_var_d2.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# design 1 : zero mean autocov ----
dt_gammatilde <- readRDS("./inst/12_mc_simulate_data/FAR/data/dt_true_gammatilde_FAR_mfBm_N=5000_zero_mean_d1.RDS")

g_gammatilde_d2 <- ggplot(dt_gammatilde, aes(x = s, y = t, z = gammatilde)) +
  geom_contour_filled(breaks = c(0, 0.733, 1.100, 1.467, 1.833, 2.200, max(dt_gammatilde$gammatilde))) +
  xlab(label = "s") +
  ylab(label = "t") +
  labs(fill = latex2exp::TeX("$\\widetilde{\\gamma}$ $(s,t)$  ")) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom",
        # legend.key.width = unit(2, "cm"),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
g_gammatilde_d2
ggsave(plot = g_gammatilde_d2, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/gammatilde_far_zero_mean_d4.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# design 1 : mean function
dt_mu <- data.table::data.table("t" = seq(0, 1, len  = 200), "mutrue" = 4 * sin(1.5 * pi * seq(0, 1, len  = 200)))

## Process variance
g_mutrue_d2 <- ggplot(data = dt_mu, aes(x = t, y = mutrue)) +
  geom_line(size = 0.5)  +
  ylim(-5, 5) +
  xlab(label = "t") +
  ylab("") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_mutrue_d2
ggsave(plot = g_mutrue_d2, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/fts_true_mean_d2.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")




# design 3 ----
dt <- readRDS("./inst/12_mc_simulate_data/FAR/data/dt_mc_common_for_variance_FAR_mfBm_N=5000_d3.RDS")
dt_var <- dt[, .("var_by_mc" = var(X)), by = c("id_mc", "tobs")]
dt_var <- dt_var[, .("var" = mean(var_by_mc)), by = "tobs"]

## Process variance
g_var_d3 <- ggplot(data = dt_var, aes(x = tobs, y = var)) +
  geom_line(size = 0.5)  +
  ylim(1.5, 6.5) +
  xlab(label = "t") +
  ylab("") +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_var_d3
ggsave(plot = g_var_d3, filename = "./inst/12_mc_simulate_data/graphs/paper_graphs/process_var_d3.png",
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
