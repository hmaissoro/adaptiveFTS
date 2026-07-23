# Adaptive functional autocorrelation function (FACF).

test_that("estimate_facf returns a classed fts_acf with a valid FACF", {
  d <- fixture_data_far(12L)
  facf <- estimate_facf(d, lag.max = 3L, n_grid = 6L)
  expect_s3_class(facf, "fts_acf")
  expect_s3_class(facf, "adaptiveFTS_est")
  expect_true(data.table::is.data.table(facf))
  expect_equal(facf$lag, 1:3)
  expect_true(all(c("lag", "norm", "facf") %in% names(facf)))
  expect_true(all(is.finite(facf$facf)))
  expect_true(all(facf$facf >= 0))                 # positive L2 norm / variance
  expect_output(summary(facf), "Functional autocorrelation")
  expect_invisible(summary(facf))
})

test_that("estimate_facf clamps lag.max and validates inputs", {
  d <- fixture_data_far(6L)
  expect_warning(estimate_facf(d, lag.max = 10L, n_grid = 5L))
  expect_error(estimate_facf(d, lag.max = 0L))
  expect_error(estimate_facf(d, t = c(-0.1, 0.5)))
})

test_that("estimate_facf autoplot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  d <- fixture_data_far(12L)
  facf <- estimate_facf(d, lag.max = 3L, n_grid = 6L)
  expect_s3_class(ggplot2::autoplot(facf), "ggplot")
})
