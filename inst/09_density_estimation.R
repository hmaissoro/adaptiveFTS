## Design-density estimation demo: the leave-one-out Parzen-Rosenblatt estimator
## (estimate_density) and the subset bandwidth selection (get_density_optimal_bw)
## used for the independent-design weights of the adaptive BLUP.

library(adaptiveFTS)
library(data.table)
library(ggplot2)

theme_set(theme_minimal(base_size = 13))
blup_theme <- theme(legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5, face = "bold"))
col_true <- "#154360"
col_pred <- "#C0392B"

data("data_far")
dt <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")

Tn0 <- dt[id_curve == 150, sort(unique(tobs))]

## Bandwidth by least-squares cross-validation on a single curve.
fit <- estimate_density(x = Tn0, bw_grid = exp(seq(log(0.02), log(0.5), length.out = 30)))
fit$h_star

ggplot(data.table(bandwidth = fit$bw_grid, cv = fit$cv_curve), aes(x = bandwidth, y = cv)) +
  geom_line(colour = col_true) +
  geom_point(colour = col_true, size = 1.2) +
  geom_vline(xintercept = fit$h_star, linetype = "dashed", colour = col_pred) +
  scale_x_log10() + blup_theme +
  labs(x = "bandwidth", y = "LSCV score", title = "Design-density LSCV")

ggplot(data.table(t = Tn0, density = fit$estimate), aes(x = t, y = density)) +
  geom_line(colour = col_true) +
  geom_point(colour = col_true, size = 1) +
  blup_theme +
  labs(x = "t", y = "density", title = "Leave-one-out density estimate")

## Design weights.
rho <- 1 / (length(Tn0) * pmax(fit$estimate, 1e-6))
rho <- rho / sum(rho)

## A single bandwidth selected on a subset of curves, reused via the fixed-h path.
h_sub <- get_density_optimal_bw(data = dt, nsubset = 30L)
ghat  <- estimate_density(x = Tn0, h = h_sub)$estimate
