# Adaptive functional BLUP through the package API.

library(adaptiveFTS)
library(data.table)
library(ggplot2)

theme_set(theme_minimal(base_size = 13))
blup_theme <- theme(legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5, face = "bold"))
col_true <- "#154360"  # true curve
col_pred <- "#C0392B"  # prediction
pred_scales <- list(
  scale_colour_manual(name = NULL,
                      values = c(Xtrue = col_true, prediction = col_pred),
                      labels = c(Xtrue = "True curve", prediction = "Adaptive prediction")),
  scale_linetype_manual(name = NULL,
                        values = c(Xtrue = "solid", prediction = "dashed"),
                        labels = c(Xtrue = "True curve", prediction = "Adaptive prediction")),
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2))))

data("data_far")
data_prepared <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")

# Hold out the last curve as the target; condition on its predecessor.
n0 <- data_prepared[, max(id_curve)]
data_train <- data_prepared[id_curve != n0]
data_test  <- data_prepared[id_curve == n0]
t0 <- data_test[, sort(unique(tobs))]

# Bandwidth grid for the adaptive risk.
N <- data_train[, length(unique(id_curve))]
lambdahat <- data_train[, .N, by = id_curve][, mean(N)]
K <- 15
b0 <- 0.5 / ((N * lambdahat) ** (1 / 2))
bK <- 0.05
bw_grid_blup <- b0 * exp((log(bK) - log(b0)) / K) ** seq_len(K)

# Density bandwidth, selected once and reused.
set.seed(1)
h_density <- get_density_optimal_bw(data = data_train, nsubset = 30L)

# One-step-ahead prediction
# tikhonov defaults to NULL, so blup_fit() selects it by cross-validation.
fit <- blup_fit(data = data_train, bw_grid = bw_grid_blup, density_bw = h_density,
                n_cv_tikhonov = 20L)
summary(fit)
pred <- predict(fit, t = t0)

dt_cmp <- merge(data_test[, .(t = tobs, Xtrue = X)], pred[, .(t, prediction)], by = "t")
dt_long <- rbind(
  dt_cmp[, .(t, curve = "Xtrue",      value = Xtrue)],
  dt_cmp[, .(t, curve = "prediction", value = prediction)])

ggplot(dt_long, aes(x = t, y = value, colour = curve, linetype = curve)) +
  geom_line(linewidth = 0.7) +
  pred_scales + blup_theme +
  labs(x = "t", y = NULL, title = "One-step-ahead adaptive BLUP")

## The Tikhonov cross-validation is carried in the fit.
cv <- fit$tikhonov_cv
ggplot(data.table(tikhonov = cv$tikhonov_grid, cv = cv$cv_curve), aes(x = tikhonov, y = cv)) +
  geom_line(colour = col_true) +
  geom_point(colour = col_true, size = 1.2) +
  geom_vline(xintercept = cv$tikhonov_star, linetype = "dashed", colour = col_pred) +
  scale_x_log10() + blup_theme +
  labs(x = expression(alpha), y = expression(CV(alpha)),
       title = "Tikhonov parameter cross-validation")

# One-call wrapper
## Reuse the selected Tikhonov to avoid re-running the cross-validation.
res <- blup(data = data_train, t = t0, tikhonov = fit$tikhonov,
            bw_grid = bw_grid_blup, density_bw = h_density)
summary(res)
stopifnot(isTRUE(all.equal(pred$prediction, res$prediction$prediction)))

# Multi-step-ahead prediction
pred_multi <- predict(fit, t = t0, horizon = 3L)

ggplot(pred_multi, aes(x = t, y = prediction, colour = factor(horizon))) +
  geom_line(linewidth = 0.7) +
  scale_colour_viridis_d(name = "Horizon", end = 0.85) +
  blup_theme +
  labs(x = "t", y = "prediction", title = "Multi-step-ahead adaptive BLUP")
