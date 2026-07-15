# Regression tests for the adaptive functional BLUP fit/predict engine and the
# design-density estimator. These lock the numerical output of the new exported
# functions (R/10_density_estimator.R, R/11_blup.R) to the committed references
# in tests/testthat/_refs/, captured on the deterministic build. They also serve
# as the sanity check that the Rcpp port (Phase 3) must reproduce.

# TOL locks values that are either pure C++ or diagonal-scaled (reproducible
# across platforms). TOL_BLAS is used for results that flow through a general
# LAPACK solve()/eigen() (the prediction and the cross-validation), which are
# not bit-reproducible across BLAS/LAPACK implementations: on Linux they differ
# from the Windows-captured references by ~1e-9, far below any real regression.
TOL <- 1e-10
TOL_BLAS <- 1e-6

test_that("design-density estimator matches references", {
  p <- ref_inputs()
  expect_equal(suppressWarnings(estimate_density(p$one$tobs, kernel_name = "epanechnikov")),
               read_ref("estimate_density_lscv"), tolerance = TOL)
  # Fixed-bandwidth path skips CV and reuses the supplied h.
  expect_equal(estimate_density(p$one$tobs, h = 0.1, kernel_name = "epanechnikov")$estimate,
               read_ref("estimate_density_fixed"), tolerance = TOL)
  expect_equal(get_density_optimal_bw(p$dt, kernel_name = "epanechnikov"),
               read_ref("get_density_optimal_bw"), tolerance = TOL)
})

test_that("blup_fit components match references", {
  p <- ref_inputs()
  fit <- blup_fit(p$dt, bw_grid = p$bwg, kernel_name = "epanechnikov")
  ref <- read_ref("blup_fit_components")
  expect_equal(fit$c0hat, ref$c0hat, tolerance = TOL)
  expect_equal(fit$V, ref$V, tolerance = TOL)
  expect_equal(fit$muhat_Tn0, ref$muhat_Tn0, tolerance = TOL)
  expect_equal(fit$rho, ref$rho, tolerance = TOL)
  expect_equal(fit$sigma2, ref$sigma2, tolerance = TOL)
  expect_equal(fit$Tn0, ref$Tn0, tolerance = TOL)
  expect_identical(fit$is_common_design, ref$is_common_design)
})

test_that("blup() prediction matches references (h = 1 and h = 2)", {
  p <- ref_inputs()
  expect_equal(blup(p$dt, t = p$tt, bw_grid = p$bwg, kernel_name = "epanechnikov"),
               read_ref("blup_predict"), tolerance = TOL_BLAS)
  expect_equal(blup(p$dt, t = p$tt, bw_grid = p$bwg, horizon = 2L, kernel_name = "epanechnikov"),
               read_ref("blup_predict_h2"), tolerance = TOL_BLAS)
})

test_that("blup() wrapper equals predict(blup_fit())", {
  p <- ref_inputs()
  fit <- blup_fit(p$dt, bw_grid = p$bwg, kernel_name = "epanechnikov")
  expect_equal(blup(p$dt, t = p$tt, bw_grid = p$bwg, kernel_name = "epanechnikov"),
               predict(fit, t = p$tt), tolerance = TOL)
})

test_that("select_tikhonov_parameter matches references", {
  p <- ref_inputs()
  tikhonov_grid_ref <- exp(seq(-2, 0, length.out = 8))
  cv <- suppressWarnings(select_tikhonov_parameter(
    p$dt, tikhonov_grid = tikhonov_grid_ref, n_val = 3L,
    bw_grid = p$bwg, kernel_name = "epanechnikov"))
  expect_equal(cv, read_ref("select_tikhonov_parameter"), tolerance = TOL_BLAS)
})
