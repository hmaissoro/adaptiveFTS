# Coverage for the Rubin-Panaretos (RP) estimators and their helpers. These are
# the slow comparison-method wrappers; tests use tiny inputs and skip the
# cross-validation variants on CRAN to keep runtime within limits.

tiny_far <- function(seed = 5L, ncur = 8L, npts = 12L) {
  set.seed(seed)
  ids <- rep(seq_len(ncur), each = npts)
  data.table::data.table(id_curve = ids, tobs = stats::runif(length(ids)), X = stats::rnorm(length(ids)))
}

test_that(".Spq_fun and .Qpq_fun return finite scalars", {
  dt <- fixture_data_far(10L)
  spq <- adaptiveFTS::.Spq_fun(dt, s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                               h = 0.1, smooth_ker = epanechnikov)
  qpq <- adaptiveFTS::.Qpq_fun(dt, s = 1/4, t = 1/2, lag = 1, p = 1, q = 1,
                               h = 0.1, optbw_mean = 0.1, dt_mean_rp = NULL,
                               smooth_ker = epanechnikov)
  expect_length(spq, 1L); expect_true(is.finite(spq))
  expect_length(qpq, 1L); expect_true(is.finite(qpq))
})

test_that("estimate_mean_rp returns finite estimates and is deterministic", {
  dt <- fixture_data_far(10L)
  m1 <- estimate_mean_rp(dt, t = c(0.25, 0.5, 0.75), h = 0.1, smooth_ker = epanechnikov)
  m2 <- estimate_mean_rp(dt, t = c(0.25, 0.5, 0.75), h = 0.1, smooth_ker = epanechnikov)
  expect_true("muhat_RP" %in% names(m1))
  expect_true(all(is.finite(m1$muhat_RP)))
  expect_equal(m1, m2)
})

test_that("estimate_autocov_rp returns finite, deterministic estimates", {
  dt <- fixture_data_far(10L)
  f <- function() estimate_autocov_rp(dt, s = c(0.25, 0.5, 0.75), t = c(0.3, 0.5, 0.7),
                                      lag = 1L, h = 0.1, optbw_mean = 0.1,
                                      dt_mean_rp = NULL, smooth_ker = epanechnikov)
  a <- f(); b <- f()
  expect_true("autocovhat_rp" %in% names(a))
  expect_true(all(is.finite(a$autocovhat_rp)))
  expect_equal(a, b)
})

test_that("estimate_mean_bw_rp runs and returns a cv_error grid", {
  skip_on_cran()
  dt <- tiny_far()
  res <- estimate_mean_bw_rp(dt, Kfold = 2L, bw_grid = c(0.1, 0.2), smooth_ker = epanechnikov)
  expect_true(all(c("h", "cv_error") %in% names(res)))
  expect_equal(nrow(res), 2L)
})

test_that("estimate_autocov_bw_rp runs and returns a cv_error grid", {
  skip_on_cran()
  # The autocovariance cross-validation wrapper is very slow (per-pair
  # Rubin-Panaretos estimation); only run it on explicit opt-in.
  if (!identical(Sys.getenv("ADAPTIVEFTS_RUN_SLOW"), "true"))
    skip("slow RP autocovariance CV; set ADAPTIVEFTS_RUN_SLOW=true to run")
  dt <- tiny_far(ncur = 6L, npts = 7L)
  res <- estimate_autocov_bw_rp(dt, Kfold = 2L, bw_grid = c(0.1, 0.2),
                                optbw_mean = 0.15, dt_mean_rp = NULL, smooth_ker = epanechnikov)
  expect_true(all(c("h", "cv_error") %in% names(res)))
  expect_equal(nrow(res), 2L)
})
