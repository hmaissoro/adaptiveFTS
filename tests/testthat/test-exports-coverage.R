# Coverage tests for exported functions not pinned by the regression references:
# kernels (R wrappers), format_data, Hurst functions, simulators, NW bandwidth
# selectors, and the real-data helpers. Properties + determinism, not exact
# numerics (those are locked in test-regression-refs.R).

# data.table functions are accessed via data.table::

test_that("exported kernel functions have the right support and values", {
  kerns <- list(epanechnikov = epanechnikov, biweight = biweight, triweight = triweight,
                tricube = tricube, triangular = triangular, uniform = uniform)
  u <- seq(-1.5, 1.5, by = 0.25)
  for (nm in names(kerns)) {
    f <- kerns[[nm]]
    val <- f(u)
    expect_length(val, length(u))
    expect_true(all(is.finite(val)), info = nm)
    # zero outside the [-1, 1] support
    expect_true(all(val[abs(u) > 1] == 0), info = nm)
    # non-negative on the support
    expect_true(all(val[abs(u) <= 1] >= 0), info = nm)
  }
  expect_equal(epanechnikov(0), 0.75)
  expect_equal(uniform(0), 0.5)
})

test_that("format_data normalizes to id_curve/tobs/X", {
  raw <- data.frame(curve = rep(1:3, each = 4),
                    tobs = rep(seq(0, 1, length.out = 4), 3),
                    X = rnorm(12))
  out <- format_data(raw, idcol = "curve", tcol = "tobs", ycol = "X")
  expect_true(data.table::is.data.table(out))
  expect_true(all(c("id_curve", "tobs", "X") %in% names(out)))
  expect_equal(length(unique(out$id_curve)), 3L)
  expect_equal(nrow(out), 12L)
})

test_that("Hurst functions return values in (0,1) on the unit interval", {
  tg <- seq(0.05, 0.95, length.out = 25)
  for (h in list(hurst_arctan(tg), hurst_linear(tg), hurst_logistic(tg))) {
    expect_length(h, length(tg))
    expect_true(all(is.finite(h)))
    expect_true(all(h > 0 & h < 1))
  }
})

test_that("simulators are reproducible and well-formed", {
  sim <- function() {
    set.seed(123)
    simulate_far(N = 6L, lambda = 30L, tdesign = "random",
                 Mdistribution = stats::rpois, tdistribution = stats::runif,
                 tcommon = NULL, hurst_fun = hurst_logistic, L = 4,
                 far_kernel = function(s, t) 9/4 * exp(-(t + 2 * s)^2),
                 far_mean = function(t) 4 * sin(1.5 * pi * t),
                 int_grid = 100L, burnin = 50L, remove_burnin = TRUE)
  }
  a <- sim(); b <- sim()
  expect_true(data.table::is.data.table(a))
  expect_true(all(c("id_curve", "tobs", "X") %in% names(a)))
  expect_equal(length(unique(a$id_curve)), 6L)
  expect_equal(a, b)                     # reproducible under set.seed
  expect_true(all(is.finite(a$X)))

  set.seed(7)
  fma <- simulate_fma(N = 5L, lambda = 25L, tdesign = "random",
                      Mdistribution = stats::rpois, tdistribution = stats::runif,
                      tcommon = NULL, hurst_fun = hurst_logistic, L = 4,
                      fma_kernel = function(s, t) exp(-(t + 2 * s)^2),
                      fma_mean = function(t) sin(2 * pi * t),
                      int_grid = 100L, burnin = 50L, remove_burnin = TRUE)
  expect_true(all(c("id_curve", "tobs", "X") %in% names(fma)))
  expect_equal(length(unique(fma$id_curve)), 5L)
})

test_that("simulate_fBm / simulate_mfBm are reproducible and finite", {
  tg <- seq(0.1, 0.9, length.out = 20)
  set.seed(11); a <- simulate_fBm(t = tg, hurst = 0.5, L = 1)
  set.seed(11); b <- simulate_fBm(t = tg, hurst = 0.5, L = 1)
  expect_equal(a, b)
  expect_true(all(is.finite(unlist(a))))

  set.seed(13); m1 <- simulate_mfBm(t = tg, hurst_fun = hurst_logistic, L = 1)
  set.seed(13); m2 <- simulate_mfBm(t = tg, hurst_fun = hurst_logistic, L = 1)
  expect_equal(m1, m2)
  expect_true(all(is.finite(unlist(m1))))
})

test_that("NW bandwidth selectors return a valid scalar bandwidth", {
  p <- ref_inputs()
  bw1 <- estimate_nw_bw(p$one$X, p$one$tobs, bw_grid = p$bwg, kernel_name = "epanechnikov")
  expect_length(bw1, 1L); expect_true(is.finite(bw1) && bw1 > 0)

  set.seed(42)
  bw2 <- get_nw_optimal_bw(p$dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
                           bw_grid = p$bwg, nsubset = 8L, kernel_name = "epanechnikov")
  expect_length(bw2, 1L); expect_true(is.finite(bw2) && bw2 > 0)
  # reproducible under the same seed
  set.seed(42)
  bw2b <- get_nw_optimal_bw(p$dt, idcol = "id_curve", tcol = "tobs", ycol = "X",
                            bw_grid = p$bwg, nsubset = 8L, kernel_name = "epanechnikov")
  expect_equal(bw2, bw2b)
})

test_that("real-data helpers return finite numeric output", {
  m <- get_real_data_mean(t = seq(0.1, 0.9, length.out = 10))
  expect_length(as.numeric(m), 10L)
  expect_true(all(is.finite(as.numeric(m))))
  k <- get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
  expect_true(all(is.finite(as.numeric(unlist(k)))))
})
