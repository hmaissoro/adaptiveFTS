# The thread setting must not change any result (each parallel output is
# computed by one thread, in a fixed summation order).

test_that("set_threads() / get_threads() set and report the thread count", {
  old <- set_threads(1)
  on.exit(set_threads(old))
  expect_equal(get_threads(), 1L)
  expect_true(is.numeric(old) && old >= 1)
  expect_error(set_threads(0), "positive integer")
  expect_error(set_threads(c(1, 2)), "positive integer")
  expect_error(set_threads(1.5), "positive integer")
})

test_that("results are identical whatever the number of threads", {
  dt <- fixture_data_far(12L)
  bwg <- seq(0.05, 0.2, length.out = 6)
  tt <- c(0.2, 0.4, 0.6, 0.8)
  run <- function() list(
    mean = estimate_mean(dt, t = tt, bw_grid = bwg),
    autocov = estimate_autocov(dt, s = rep(tt, 4), t = rep(tt, each = 4), lag = 1, bw_grid = bwg),
    cov = estimate_autocov(dt, s = rep(tt, 4), t = rep(tt, each = 4), lag = 0, bw_grid = bwg),
    blup = blup(dt, t = tt, tikhonov = 1e-3, bw_grid = bwg)$prediction)
  old <- set_threads(1)
  on.exit(set_threads(old))
  r1 <- run()
  set_threads(2)
  r2 <- run()
  set_threads(old)
  r3 <- run()
  expect_identical(r1, r2)
  expect_identical(r1, r3)
})
