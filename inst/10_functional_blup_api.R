## Adaptive functional BLUP demo using the package API: get_density_optimal_bw /
## estimate_density, blup_fit / predict, the blup() wrapper, and
## select_tikhonov_parameter. Mirrors the steps of inst/08_functional_blup.R.

library(adaptiveFTS)
library(data.table)
library(ggplot2)

data("data_far")
data_prepared <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")

## Hold out the last curve as the target; condition on its predecessor.
n0 <- data_prepared[, max(id_curve)]
data_train <- data_prepared[id_curve != n0]
data_test  <- data_prepared[id_curve == n0]
t0 <- data_test[, sort(unique(tobs))]

## Bandwidth grid for the adaptive risk (blup_fit uses a comparable default when
## bw_grid = NULL).
N <- data_train[, length(unique(id_curve))]
lambdahat <- data_train[, .N, by = id_curve][, mean(N)]
K <- 15
b0 <- 0.5 / ((N * lambdahat) ** (1 / 2))
bK <- 0.05
bw_grid_blup <- b0 * exp((log(bK) - log(b0)) / K) ** seq_len(K)

## Density bandwidth, selected once and reused (set.seed for a reproducible subset).
set.seed(1)
h_density <- get_density_optimal_bw(data = data_train, nsubset = 30L)

## ---- One-step-ahead prediction ------------------------------------------
fit  <- blup_fit(data = data_train, bw_grid = bw_grid_blup, density_bw = h_density)
pred <- predict(fit, t = t0)

dt_cmp <- merge(data_test[, .(t = tobs, Xtrue = X)], pred[, .(t, prediction)], by = "t")
dt_long <- rbind(
  dt_cmp[, .(t, Quantity = "prediction", value = prediction)],
  dt_cmp[, .(t, Quantity = "Xtrue",      value = Xtrue)])

ggplot(dt_long, aes(x = t, y = value, colour = Quantity)) +
  geom_line() + theme_minimal() + theme(legend.position = "top") +
  labs(x = "t", y = NULL, title = "Adaptive BLUP: one-step-ahead prediction")

## blup() is the one-call equivalent of blup_fit() + predict().
pred_wrap <- blup(data = data_train, t = t0, bw_grid = bw_grid_blup, density_bw = h_density)
stopifnot(isTRUE(all.equal(pred$prediction, pred_wrap$prediction)))

## ---- Multi-step-ahead prediction ----------------------------------------
pred_multi <- predict(fit, t = t0, horizon = 3L)

ggplot(pred_multi, aes(x = t, y = prediction, colour = factor(horizon))) +
  geom_line() + theme_minimal() + theme(legend.position = "top") +
  labs(x = "t", y = "prediction", colour = "horizon",
       title = "Adaptive BLUP: multi-step-ahead prediction")

## ---- Tikhonov parameter by cross-validation -----------------------------
## Slow: a rolling one-step-ahead CV over the last n_val curves.
cv <- select_tikhonov_parameter(data = data_train, method = "cv",
                                n_val = 30L, bw_grid = bw_grid_blup)
cv$tikhonov_star

ggplot(data.table(tikhonov = cv$tikhonov_grid, cv = cv$cv_curve),
       aes(x = tikhonov, y = cv)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = cv$tikhonov_star, linetype = 2, colour = "red") +
  scale_x_log10() + theme_minimal() +
  labs(x = expression(alpha), y = expression(CV(alpha)),
       title = "Holdout CV for the Tikhonov parameter")

fit_cv  <- blup_fit(data = data_train, tikhonov = cv$tikhonov_star,
                    bw_grid = bw_grid_blup, density_bw = h_density)
pred_cv <- predict(fit_cv, t = t0)

dt_cmp_cv <- merge(data_test[, .(t = tobs, Xtrue = X)], pred_cv[, .(t, prediction)], by = "t")
ggplot(dt_cmp_cv, aes(x = t)) +
  geom_line(aes(y = prediction, colour = "prediction")) +
  geom_line(aes(y = Xtrue, colour = "Xtrue")) +
  theme_minimal() + theme(legend.position = "top") +
  labs(x = "t", y = NULL, colour = NULL,
       title = "Adaptive BLUP: prediction with CV-selected Tikhonov")
