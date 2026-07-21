## Adaptive functional BLUP demo using the package API: get_density_optimal_bw /
## estimate_density, blup_fit / predict (with automatic Tikhonov selection), the
## blup() wrapper, and the summary methods. Mirrors inst/08_functional_blup.R.

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
## tikhonov defaults to NULL, so blup_fit() selects it by cross-validation.
fit <- blup_fit(data = data_train, bw_grid = bw_grid_blup, density_bw = h_density,
                n_cv_tikhonov = 20L)
summary(fit)
pred <- predict(fit, t = t0)

dt_cmp <- merge(data_test[, .(t = tobs, Xtrue = X)], pred[, .(t, prediction)], by = "t")
dt_long <- rbind(
  dt_cmp[, .(t, Quantity = "prediction", value = prediction)],
  dt_cmp[, .(t, Quantity = "Xtrue",      value = Xtrue)])

ggplot(dt_long, aes(x = t, y = value, colour = Quantity)) +
  geom_line() + theme_minimal() + theme(legend.position = "top") +
  labs(x = "t", y = NULL, title = "Adaptive BLUP: one-step-ahead prediction")

## The Tikhonov cross-validation is carried in the fit.
cv <- fit$tikhonov_cv
ggplot(data.table(tikhonov = cv$tikhonov_grid, cv = cv$cv_curve), aes(x = tikhonov, y = cv)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = cv$tikhonov_star, linetype = 2, colour = "red") +
  scale_x_log10() + theme_minimal() +
  labs(x = expression(alpha), y = expression(CV(alpha)),
       title = "Holdout CV for the Tikhonov parameter")

## ---- One-call wrapper ---------------------------------------------------
## Reuse the selected Tikhonov to avoid re-running the cross-validation.
res <- blup(data = data_train, t = t0, tikhonov = fit$tikhonov,
            bw_grid = bw_grid_blup, density_bw = h_density)
summary(res)
stopifnot(isTRUE(all.equal(pred$prediction, res$prediction$prediction)))

## ---- Multi-step-ahead prediction ----------------------------------------
pred_multi <- predict(fit, t = t0, horizon = 3L)

ggplot(pred_multi, aes(x = t, y = prediction, colour = factor(horizon))) +
  geom_line() + theme_minimal() + theme(legend.position = "top") +
  labs(x = "t", y = "prediction", colour = "horizon",
       title = "Adaptive BLUP: multi-step-ahead prediction")
