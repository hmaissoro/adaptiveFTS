## Adaptive functional BLUP — package API demo.
##
## Mirrors the logic and steps of inst/08_functional_blup.R, but using the
## exported package API: the design-density estimator
## (get_density_optimal_bw() / estimate_density()), the fit/predict engine
## (blup_fit() / predict()), the one-call wrapper blup(), and the Tikhonov
## parameter selection select_tikhonov_parameter().

library(adaptiveFTS)
library(data.table)
library(ggplot2)

## ---- Data ---------------------------------------------------------------
data("data_far")
data_prepared <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")

## Hold out the last curve as the target to predict; condition on its predecessor.
n0 <- data_prepared[, max(id_curve)]
data_train <- data_prepared[id_curve != n0]
data_test  <- data_prepared[id_curve == n0]

## Prediction points: the observation grid of the held-out curve.
t0 <- data_test[, sort(unique(tobs))]

## Bandwidth grid for the adaptive mean/(auto)covariance risk (as in inst/08;
## blup_fit() computes a comparable default when bw_grid = NULL).
N <- data_train[, length(unique(id_curve))]
lambdahat <- data_train[, .N, by = id_curve][, mean(N)]
K <- 15
b0 <- 0.5 / ((N * lambdahat) ** (1 / 2))
bK <- 0.05
bw_grid_blup <- b0 * exp((log(bK) - log(b0)) / K) ** seq_len(K)

## ---- Design density (independent design) --------------------------------
## The design weights rely on a single density bandwidth selected once on a
## subset of curves and reused for every curve (set.seed for a reproducible
## subset).
set.seed(1)
h_density <- get_density_optimal_bw(data = data_train, nsubset = 30L)
ghat <- estimate_density(x = t0, h = h_density)$estimate
plot(t0, ghat, type = "l", main = "Estimated design density", xlab = "t", ylab = "g(t)")
## ---- One-step-ahead prediction ------------------------------------------
fit <- blup_fit(data = data_train, bw_grid = bw_grid_blup,
                density_bw = h_density, kernel_name = "epanechnikov")
pred <- predict(fit, t = t0)

## Compare the prediction with the (held-out) truth.
dt_cmp <- merge(data_test[, .(t = tobs, Xtrue = X)], pred[, .(t, prediction)], by = "t")
dt_long <- rbind(
  dt_cmp[, .(t, Quantity = "prediction", value = prediction)],
  dt_cmp[, .(t, Quantity = "Xtrue",      value = Xtrue)])

ggplot(dt_long, aes(x = t, y = value, colour = Quantity)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "t", y = NULL, title = "Adaptive BLUP: one-step-ahead prediction")

## The one-call wrapper is equivalent to blup_fit() + predict().
pred_wrap <- blup(data = data_train, t = t0, bw_grid = bw_grid_blup, density_bw = h_density)
stopifnot(isTRUE(all.equal(pred$prediction, pred_wrap$prediction)))

## ---- Multi-step-ahead prediction ----------------------------------------
## horizon = 3 returns the horizon 1, 2 and 3 predictions stacked, with a
## `horizon` column; the operators are held fixed and only the conditioning
## curve is fed back at each step.
pred_multi <- predict(fit, t = t0, horizon = 3L)

ggplot(pred_multi, aes(x = t, y = prediction, colour = factor(horizon))) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(x = "t", y = "prediction", colour = "horizon",
       title = "Adaptive BLUP: multi-step-ahead prediction")

## ---- Selecting the Tikhonov parameter by cross-validation ---------------
## (Slow: a rolling one-step-ahead CV over the last n_val curves.)
if (FALSE) {
  cv <- select_tikhonov_parameter(
    data = data_train, method = "cv",
    tikhonov_grid = exp(seq(-5, 0, length.out = 25)),
    n_val = 30L, bw_grid = bw_grid_blup)
  cv$tikhonov_star

  ggplot(data.table(tikhonov = cv$tikhonov_grid, cv = cv$cv_curve),
         aes(x = tikhonov, y = cv)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = cv$tikhonov_star, linetype = 2, colour = "red") +
    scale_x_log10() +
    theme_minimal() +
    labs(x = expression(alpha), y = expression(CV(alpha)),
         title = "Holdout CV for the Tikhonov parameter")

  ## Refit with the selected parameter and predict.
  fit_cv <- blup_fit(data = data_train, tikhonov = cv$tikhonov_star, bw_grid = bw_grid_blup)
  pred_cv <- predict(fit_cv, t = t0)
  
  dt_cmp_cv <- merge(data_test[, .(t = tobs, Xtrue = X)], pred_cv[, .(t, prediction)], by = "t")
  ggplot(dt_cmp_cv, aes(x = t)) +
    geom_line(aes(y = prediction, colour = "prediction")) +
    geom_line(aes(y = Xtrue, colour = "Xtrue")) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "t", y = NULL, colour = NULL,
         title = "Adaptive BLUP: one-step-ahead prediction with CV-selected Tikhonov")

  
}
