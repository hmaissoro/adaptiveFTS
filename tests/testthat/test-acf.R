# Adaptive functional autocorrelation function (FACF).

test_that("estimate_acf returns a classed fts_acf with a valid FACF", {
  d <- fixture_data_far(12L)
  acf <- estimate_acf(d, lag.max = 3L, n_grid = 6L)
  expect_s3_class(acf, "fts_acf")
  expect_s3_class(acf, "adaptiveFTS_est")
  expect_true(data.table::is.data.table(acf))
  expect_equal(acf$lag, 1:3)
  expect_true(all(c("lag", "norm", "facf") %in% names(acf)))
  expect_true(all(is.finite(acf$facf)))
  expect_true(all(acf$facf >= 0))                 # positive L2 norm / variance
  expect_output(summary(acf), "Functional autocorrelation")
  expect_invisible(summary(acf))
})

test_that("estimate_acf clamps lag.max and validates inputs", {
  d <- fixture_data_far(6L)
  expect_warning(estimate_acf(d, lag.max = 10L, n_grid = 5L))
  expect_error(estimate_acf(d, lag.max = 0L))
  expect_error(estimate_acf(d, t = c(-0.1, 0.5)))
})

test_that("estimate_acf autoplot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  d <- fixture_data_far(12L)
  acf <- estimate_acf(d, lag.max = 3L, n_grid = 6L)
  expect_s3_class(ggplot2::autoplot(acf), "ggplot")
})
