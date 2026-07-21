# Behaviour and edge-case tests for the adaptive BLUP engine and the design-
# density estimator (properties and error handling, not exact numerics, which
# are locked in test-blup-refs.R).

test_that("estimate_density validates input and honours the fixed-bandwidth path", {
  x <- fixture_data_far(3L)[id_curve == 1L][order(tobs), tobs]
  expect_error(estimate_density(x[1]), "at least two")
  expect_error(estimate_density(x, h = -1), "positive")
  expect_error(estimate_density(x, kernel_name = "nope"))

  fit_cv <- suppressWarnings(estimate_density(x, kernel_name = "epanechnikov"))
  expect_false(is.null(fit_cv$cv_curve))
  expect_length(fit_cv$estimate, length(x))

  fit_fixed <- estimate_density(x, h = fit_cv$h_star)
  expect_null(fit_fixed$cv_curve)
  expect_equal(fit_fixed$h_star, fit_cv$h_star)
  # Fixed h reproduces the LSCV estimate when given the selected bandwidth.
  expect_equal(fit_fixed$estimate, fit_cv$estimate)
  expect_true(all(is.finite(fit_fixed$estimate) & fit_fixed$estimate >= 0))
})

test_that("get_density_optimal_bw is deterministic and validates nsubset", {
  dt <- fixture_data_far(12L)
  N <- length(unique(dt$id_curve))
  expect_error(get_density_optimal_bw(dt, nsubset = N + 5L), "positive integer")
  expect_error(get_density_optimal_bw(dt, nsubset = 2.5), "positive integer")

  h_all <- get_density_optimal_bw(dt)
  expect_length(h_all, 1L)
  expect_true(is.finite(h_all) && h_all > 0)

  set.seed(1); h1 <- get_density_optimal_bw(dt, nsubset = 6L)
  set.seed(1); h2 <- get_density_optimal_bw(dt, nsubset = 6L)
  expect_equal(h1, h2)
})

test_that("blup_fit returns a well-formed object", {
  dt <- fixture_data_far(12L)
  bwg <- seq(0.05, 0.2, length.out = 6)
  fit <- blup_fit(dt, bw_grid = bwg)
  expect_s3_class(fit, "blup_fit")
  expect_true(all(c("V", "c0hat", "muhat_Tn0", "rho", "Tn0",
                    "opt_mean", "opt_cov", "opt_autocov") %in% names(fit)))
  M <- fit$Mn0
  expect_equal(dim(fit$V), c(M, M))
  expect_equal(dim(fit$c0hat), c(M, M))
  # V is symmetric positive-definite: C0 is projected onto the PSD cone, so
  # V = D half %*% C0 %*% D half + noise + tikhonov I is PD.
  expect_equal(fit$V, t(fit$V), tolerance = 1e-10)
  expect_true(all(eigen(fit$V, symmetric = TRUE, only.values = TRUE)$values > 0))
  # design weights sum to one
  expect_equal(sum(fit$rho), 1, tolerance = 1e-10)
})

test_that("predict.blup_fit handles t, newdata and h, and validates them", {
  dt <- fixture_data_far(12L)
  bwg <- seq(0.05, 0.2, length.out = 6)
  fit <- blup_fit(dt, bw_grid = bwg)
  tt <- c(0.25, 0.5, 0.75)

  p1 <- predict(fit, t = tt)
  expect_true(data.table::is.data.table(p1))
  expect_equal(names(p1), c("horizon", "t", "muhat", "prediction"))
  expect_equal(p1$t, tt)
  expect_true(all(p1$horizon == 1L))
  expect_true(all(is.finite(p1$prediction)))

  # default t is the conditioning-curve design
  pd <- predict(fit)
  expect_equal(pd$t, fit$Tn0)

  # multi-step returns every intermediate horizon, and its horizon-1 block
  # equals the standalone one-step prediction.
  p2 <- predict(fit, t = tt, horizon = 2L)
  expect_equal(nrow(p2), 2L * length(tt))
  expect_equal(sort(unique(p2$horizon)), c(1L, 2L))
  expect_equal(p2[horizon == 1L, prediction], p1$prediction, tolerance = 1e-10)
  expect_true(all(is.finite(p2$prediction)))

  # newdata overrides the conditioning values
  pn <- predict(fit, t = tt, newdata = fit$Yn0)
  expect_equal(pn$prediction, p1$prediction, tolerance = 1e-10)

  expect_error(predict(fit, t = 1.5), "between 0 and 1")
  expect_error(predict(fit, horizon = 0L), "positive integer")
  expect_error(predict(fit, newdata = c(1, 2, 3)), "length equal")
})

test_that("blup() equals predict(blup_fit()) and cv_blup_alpha validates n_cv_tikhonov", {
  dt <- fixture_data_far(12L)
  bwg <- seq(0.05, 0.2, length.out = 6)
  tt <- c(0.3, 0.6)
  fit <- blup_fit(dt, bw_grid = bwg)
  expect_equal(blup(dt, t = tt, bw_grid = bwg), predict(fit, t = tt), tolerance = 1e-10)

  expect_error(select_tikhonov_parameter(dt, n_cv_tikhonov = 100L), "smaller than the number")
  cv <- suppressWarnings(select_tikhonov_parameter(
    dt, tikhonov_grid = exp(seq(-2, 0, length.out = 6)), n_cv_tikhonov = 3L, bw_grid = bwg))
  expect_true(all(c("tikhonov_star", "cv_curve", "cv_matrix", "val_ids") %in% names(cv)))
  expect_length(cv$cv_curve, 6L)
  expect_equal(dim(cv$cv_matrix), c(3L, 6L))
  expect_true(is.finite(cv$tikhonov_star))
})
