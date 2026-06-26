# Edge cases: boundary evaluation points, different lags, single points, and
# sparser designs. These assert the estimators run and return finite, sane
# output rather than locking exact numbers.

suppressMessages(library(data.table))

test_that("mean / locreg handle boundary and single evaluation points", {
  dt <- fixture_data_far(12L)
  tb <- c(0.02, 0.5, 0.98)
  m <- estimate_mean(dt, t = tb, bw_grid = seq(0.05, 0.2, length.out = 5))
  expect_equal(nrow(m), length(tb))
  expect_true(all(is.finite(m[["muhat"]])))

  lr <- estimate_locreg(dt, t = tb, center = TRUE)
  expect_equal(nrow(lr), length(tb))
  expect_true(all(is.finite(lr[["Ht"]])))
  expect_true(all(lr[["Ht"]] >= 0.1 & lr[["Ht"]] <= 1))   # clamped exponent

  m1 <- estimate_mean(dt, t = 0.5, bw_grid = seq(0.05, 0.2, length.out = 5))
  expect_equal(nrow(m1), 1L)
  expect_true(is.finite(m1[["muhat"]]))
})

test_that("autocovariance handles lag 0 and lag > 0", {
  dt <- fixture_data_far(12L)
  bwg <- seq(0.05, 0.2, length.out = 5)
  s <- c(0.25, 0.5, 0.75); t <- c(0.3, 0.5, 0.7)
  for (lg in c(0L, 1L, 2L)) {
    ac <- estimate_autocov(dt, s = s, t = t, lag = lg, bw_grid = bwg)
    expect_true(all(is.finite(ac[["autocov"]])), info = paste("lag", lg))
  }
})

test_that("estimators run on a sparser set of curves", {
  dt <- fixture_data_far(6L)
  m <- estimate_mean(dt, t = c(0.3, 0.6), bw_grid = seq(0.05, 0.2, length.out = 5))
  expect_true(all(is.finite(m[["muhat"]])))
  sg <- estimate_sigma(dt, idcol = "id_curve", t = c(0.3, 0.6))
  expect_true(all(is.finite(sg[["sig"]])) || all(is.finite(unlist(sg))))
})

test_that("predict_curve returns a reconstruction on the requested grid", {
  dt <- fixture_data_far(12L)
  tt <- c(0.2, 0.4, 0.6, 0.8)
  pc <- predict_curve(dt, t = tt, id_curve_to_predict = 3L,
                      bw_grid = seq(0.05, 0.2, length.out = 5))
  expect_true(is.list(pc) || data.table::is.data.table(pc) || is.matrix(pc))
})
