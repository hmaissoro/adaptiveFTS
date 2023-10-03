library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Simulation parameters ----
N <- c(100L, 200L, 400L)
lambda <- c(25L, 50L, 100L, 200L, 300L)
sig <- 0.5
mc <- 500
t0 <- seq(0.2, 0.8, len = 6)

## Logistic constant hurst function
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.8,
                 change_point_position = 0.5, slope = 5)
}

# Plot true local regularity parameters ----
## Generate data
dt_Hlogistic <- data.table::data.table(
  "t" = seq(0.01, 0.99, len = 100),
  "Ht" = Hlogistic(t = seq(0.01, 0.99, len = 100))
)
dt_Lconstant <- data.table::data.table(
  "t" = seq(0.01, 0.99, len = 100),
  "Lt" = rep(4, 100)
)

## Plot and save
figures_path <- "../../../report/learning-smmoothness/Learning-smoothness/figures/"
theme_set(theme_minimal())

### Local exponent
ggplot(dt_Hlogistic, aes(x = t, y = Ht)) +
  geom_line(size = 1.5) +
  ylim(0.4, 0.8) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$H_t$")) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "local_exponent.png"), units = "px", dpi = 300)

### Hölder function
ggplot(dt_Lconstant, aes(x = t, y = Lt)) +
  geom_line(size = 1.5) +
  ylim(3, 5) +
  xlab(label = "t") +
  ylab(label = latex2exp::TeX("$L_t^2$")) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "holder_constant.png"), units = "px", dpi = 300)

# Delta expo vs Delta poly ----
dt_locreg <- readRDS("./inst/locreg_estimates/many_delta/dt_Hlogistic_bw_sig05_ml.RDS")

## Estimates of the local regularity parameters ----
dt_locreg[, c("N", "lambda") := .(as.factor(N), as.factor(lambda))]

### Plot function
gplot <- function(data = dt_locreg, i, delta_formula = "exponential", param = "H"){
  data[, c("N", "lambda") := .(as.factor(N), as.factor(lambda))]

  if (param == "H") {
    var <- "H_mc"
    ylabel <- latex2exp::TeX("$H_t$")
    yconst <- Hlogistic(t0[i])
  } else if (param == "L") {
    var <- "L_mc"
    ylabel <- latex2exp::TeX("$L_t^2$")
    yconst <- 4
  } else if (param == "Ht") {
    var <- "Ht"
    ylabel <- latex2exp::TeX("$H_t$")
    yconst <- Hlogistic(t0[i])
  } else if (param == "Lt") {
    var <- "Lt"
    ylabel <- latex2exp::TeX("$L_t^2$")
    yconst <- 4
  } else if (param == "gH_mc") {
    var <- "gH_mc"
    ylabel <- latex2exp::TeX("$\\gamma$ for $H_t$")
    yconst <- 0
  } else if (param == "gL_mc") {
    var <- "gL_mc"
    ylabel <- latex2exp::TeX("$\\gamma$ for $L_t^2$")
    yconst <- 0
  } else if (param == "DeltaH_mc") {
    var <- "DeltaH_mc"
    ylabel <- latex2exp::TeX("$\\Delta$ for $H_t$")
    yconst <- 0
  } else if (param == "DeltaL_mc") {
    var <- "DeltaL_mc"
    ylabel <- latex2exp::TeX("$\\Delta$ for $L_t^2$")
    yconst <- 0
  }
  if (delta_formula == "exponential"){
    title_exp <- paste0("$\\Delta=e^{-\\log(\\lambda)^\\gamma}$, t=", t0[i] , ", $H_t$=", round(Hlogistic(t0[i]), 3), ", $L_t^2 =4$")
    dt <- data[t == t0[i] & delta_formula == "exponential"]
  } else if (delta_formula == "polynomial"){
    title_exp <- paste0("$\\Delta=\\lambda^{-\\gamma}$, t=", t0[i] , ", $H_t$=", round(Hlogistic(t0[i]), 3), ", $L_t^2 =4$")
    dt <- data[t == t0[i] & delta_formula == "polynomial"]
  } else if (delta_formula == "delta_expo"){
    title_exp <- paste0("$\\Delta=e^{-\\log(\\lambda)^\\gamma}$, t=", t0[i] , ", $H_t$=", round(Hlogistic(t0[i]), 3), ", $L_t^2 =4$")
    dt <- data[t == t0[i]]
  } else if (delta_formula == "delta_poli"){
    title_exp <- paste0("$\\Delta=\\lambda^{-\\gamma}$, t=", t0[i] , ", $H_t$=", round(Hlogistic(t0[i]), 3), ", $L_t^2 =4$")
    dt <- data[t == t0[i]]
  }

  gplt <- ggplot(dt, aes(x = lambda, y = get(var), fill = N)) +
    geom_boxplot() +
    ggtitle(latex2exp::TeX(title_exp)) +
    xlab(latex2exp::TeX("$\\lambda$")) +
    ylab(ylabel) +
    geom_hline(aes(yintercept = yconst), colour = "black") +
    scale_fill_grey() +
    theme_minimal() +
    theme(legend.position = "top", title = element_text(size = 7))
  return(gplt)
}

### Plot Local exponent H ----
gHexp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "H"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "H"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "H"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "H"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "H"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "H"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/H_Delta_expo.png"),
  plot = gHexp, units = "px", dpi = 300, bg = "white")

gHpoli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "H"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "H"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "H"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "H"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "H"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "H"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/H_Delta_poli.png"),
  plot = gHpoli, units = "px", dpi = 300, bg = "white")

### Hölder constant ----
gLexp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "L"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "L"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "L"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "L"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "L"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "L"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/L_Delta_expo.png"),
  plot = gLexp, units = "px", dpi = 300, bg = "white")

gLpoli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "L"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "L"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "L"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "L"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "L"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "L"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/L_Delta_poli.png"),
  plot = gLpoli, units = "px", dpi = 300, bg = "white")

### Delta values ----
#### For H
gDeltaH_exp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "DeltaH_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/deltaH_Delta_expo.png"),
  plot = gDeltaH_exp, units = "px", dpi = 300, bg = "white")

gDeltaH_poli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "DeltaH_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "DeltaH_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/deltaH_Delta_poli.png"),
  plot = gDeltaH_poli, units = "px", dpi = 300, bg = "white")

### For L
gDeltaL_exp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "DeltaL_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/deltaL_Delta_expo.png"),
  plot = gDeltaL_exp, units = "px", dpi = 300, bg = "white")

gDeltaL_poli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "DeltaL_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "DeltaL_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/deltaL_Delta_poli.png"),
  plot = gDeltaL_poli, units = "px", dpi = 300, bg = "white")
### Delta values ----
#### For H
ggH_exp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "gH_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "gH_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "gH_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "gH_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "gH_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "gH_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/gH_Delta_expo.png"),
  plot = ggH_exp, units = "px", dpi = 300, bg = "white")

ggH_poli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "gH_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "gH_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "gH_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "gH_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "gH_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "gH_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/gH_Delta_poli.png"),
  plot = ggH_poli, units = "px", dpi = 300, bg = "white")

### For L
ggL_exp <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "exponential", param = "gL_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "exponential", param = "gL_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "exponential", param = "gL_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "exponential", param = "gL_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "exponential", param = "gL_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "exponential", param = "gL_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/gL_Delta_expo.png"),
  plot = ggL_exp, units = "px", dpi = 300, bg = "white")

ggL_poli <- gridExtra::grid.arrange(
  gplot(data = dt_locreg, i = 1, delta_formula = "polynomial", param = "gL_mc"),
  gplot(data = dt_locreg, i = 2, delta_formula = "polynomial", param = "gL_mc"),
  gplot(data = dt_locreg, i = 3, delta_formula = "polynomial", param = "gL_mc"),
  gplot(data = dt_locreg, i = 4, delta_formula = "polynomial", param = "gL_mc"),
  gplot(data = dt_locreg, i = 5, delta_formula = "polynomial", param = "gL_mc"),
  gplot(data = dt_locreg, i = 6, delta_formula = "polynomial", param = "gL_mc"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/gL_Delta_poli.png"),
  plot = ggL_poli, units = "px", dpi = 300, bg = "white")

# Graphs for one Delta ----
dt_locreg_delta_expo <- readRDS("./inst/locreg_estimates/best_delta/dt_Hlogistic_bw_sig05_Delta_expo.RDS")
dt_locreg_delta_expo[, c("N", "lambda") := .(as.factor(N), as.factor(lambda))]

## Delta exponential
gHexp_best <- gridExtra::grid.arrange(
  gplot(data = dt_locreg_delta_expo, i = 1, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_expo, i = 2, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_expo, i = 3, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_expo, i = 4, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_expo, i = 5, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_expo, i = 6, delta_formula = "delta_expo", param = "Ht"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/best_delta/H_bestDelta_expo.png"),
  plot = gHexp_best, units = "px", dpi = 300, bg = "white")

gLexp_best <- gridExtra::grid.arrange(
  gplot(data = dt_locreg_delta_expo, i = 1, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_expo, i = 2, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_expo, i = 3, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_expo, i = 4, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_expo, i = 5, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_expo, i = 6, delta_formula = "delta_expo", param = "Lt"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/best_delta/L_bestDelta_expo.png"),
  plot = gLexp_best, units = "px", dpi = 300, bg = "white")

## Delta constant
dt_locreg_delta_const <- readRDS("./inst/locreg_estimates/best_delta/dt_Hlogistic_bw_sig05_Delta015.RDS")
dt_locreg_delta_const[, c("N", "lambda") := .(as.factor(N), as.factor(lambda))]

gHconstant_best <- gridExtra::grid.arrange(
  gplot(data = dt_locreg_delta_const, i = 1, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_const, i = 2, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_const, i = 3, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_const, i = 4, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_const, i = 5, delta_formula = "delta_expo", param = "Ht"),
  gplot(data = dt_locreg_delta_const, i = 6, delta_formula = "delta_expo", param = "Ht"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/H_bestDelta_constant.png"),
  plot = gHconstant_best, units = "px", dpi = 300, bg = "white")

gLconstant_best <- gridExtra::grid.arrange(
  gplot(data = dt_locreg_delta_const, i = 1, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_const, i = 2, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_const, i = 3, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_const, i = 4, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_const, i = 5, delta_formula = "delta_expo", param = "Lt"),
  gplot(data = dt_locreg_delta_const, i = 6, delta_formula = "delta_expo", param = "Lt"),
  ncol = 3, nrow = 2
)
ggsave(
  filename = file.path("./inst/locreg_estimates/many_delta/L_bestDelta_constant.png"),
  plot = gHconstant_best, units = "px", dpi = 300, bg = "white")



