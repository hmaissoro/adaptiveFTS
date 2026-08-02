# `intercept_var` (renamed from `shift_var`): a per-curve random intercept added
# to the innovation. The intercept is constant in `t`, so it cancels in the
# increments and leaves the local regularity (H_t, L_t) untouched; it only moves
# the realisation on the ordinate axis. These tests pin the four things that
# matter: the default path is bit-identical, the increments are invariant, the
# intercept is calibrated to sqrt(L * intercept_var), and `tied = TRUE` refuses it.

# data.table functions are accessed via data.table::

TG <- seq(0.05, 1, length.out = 12)
LC <- 4
IV <- 0.05

sim_far <- function(...) {
  simulate_far(N = 4L, lambda = 25L, design = "random",
               M_distribution = stats::rpois, t_distribution = stats::runif,
               t_common = NULL, hurst_fun = hurst_logistic, L = LC,
               far_kernel = function(s, t) 9 / 4 * exp(-(t + 2 * s)^2),
               far_mean = function(t) 4 * sin(1.5 * pi * t),
               n_int_grid = 60L, n_burnin = 40L, remove_burnin = TRUE, ...)
}

sim_fma <- function(...) {
  simulate_fma(N = 4L, lambda = 25L, design = "random",
               M_distribution = stats::rpois, t_distribution = stats::runif,
               t_common = NULL, hurst_fun = hurst_logistic, L = LC,
               fma_kernel = function(s, t) 9 / 4 * exp(-(t + 2 * s)^2),
               fma_mean = function(t) 4 * sin(1.5 * pi * t),
               n_int_grid = 60L, n_burnin = 40L, remove_burnin = TRUE, ...)
}

test_that("intercept_var = 0 reproduces the default path bit-for-bit", {
  # The regression guard: under the zero default no extra draw is taken, so the
  # RNG stream and every value are unchanged from before intercept_var existed.
  set.seed(101); a <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC, tied = FALSE)
  set.seed(101); b <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                                    intercept_var = 0, tied = FALSE)
  expect_identical(a[["mfBm"]], b[["mfBm"]])

  set.seed(202); fa <- sim_far()
  set.seed(202); fb <- sim_far(intercept_var = 0)
  expect_identical(fa[["X"]], fb[["X"]])
  expect_identical(fa[["tobs"]], fb[["tobs"]])

  set.seed(303); ma <- sim_fma()
  set.seed(303); mb <- sim_fma(intercept_var = 0)
  expect_identical(ma[["X"]], mb[["X"]])
})

test_that("the intercept cancels in the increments", {
  # This is the property that leaves H_t and L_t unchanged: with a shared seed
  # the two paths differ by a constant, so their first differences coincide.
  set.seed(7); p0 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                                   intercept_var = 0, tied = FALSE)[["mfBm"]]
  set.seed(7); p1 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                                   intercept_var = IV, tied = FALSE)[["mfBm"]]
  expect_equal(diff(p1), diff(p0), tolerance = 1e-12)
  expect_equal(stats::sd(p1 - p0), 0, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(p1, p0)))
})

test_that("the intercept is calibrated to sqrt(L * intercept_var)", {
  # Exact check: the path draw consumes the same stream in both calls, so the
  # displacement is sqrt(L * intercept_var) times the very next standard normal.
  set.seed(21)
  p0 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                      intercept_var = 0, tied = FALSE)[["mfBm"]]
  z <- stats::rnorm(1)
  set.seed(21)
  p1 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                      intercept_var = IV, tied = FALSE)[["mfBm"]]
  expect_equal(mean(p1 - p0), sqrt(LC * IV) * z, tolerance = 1e-10)

  # Distributional cross-check on the same quantity.
  delta <- vapply(seq_len(300), function(r) {
    set.seed(r)
    x0 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                        intercept_var = 0, tied = FALSE)[["mfBm"]]
    set.seed(r)
    x1 <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                        intercept_var = IV, tied = FALSE)[["mfBm"]]
    mean(x1 - x0)
  }, numeric(1))
  expect_equal(stats::sd(delta), sqrt(LC * IV), tolerance = 0.2)
})

test_that("the intercept widens the spread at the smallest design point", {
  # The marginal variance there is L * (Var xi(t1) + intercept_var), so the
  # variance *difference* is L * intercept_var; the level itself is not.
  spread <- function(iv) {
    vapply(seq_len(400), function(r) {
      set.seed(2000 + r)
      simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                    intercept_var = iv, tied = FALSE)[["mfBm"]][1]
    }, numeric(1))
  }
  s0 <- spread(0)
  s1 <- spread(IV)
  expect_gt(stats::sd(s1), stats::sd(s0))
  expect_equal(stats::var(s1) - stats::var(s0), LC * IV, tolerance = 0.35)
})

test_that("intercept_var is refused when tied = TRUE", {
  # A tied-down path with a non-zero intercept is neither tied down at the
  # origin nor an intercept-shifted mfBm, so the intercept is dropped.
  set.seed(11)
  a <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                     intercept_var = 0, tied = TRUE)
  set.seed(11)
  expect_warning(
    b <- simulate_mfBm(t = TG, hurst_fun = hurst_logistic, L = LC,
                       intercept_var = IV, tied = TRUE),
    "ignored when 'tied = TRUE'")
  expect_identical(a[["mfBm"]], b[["mfBm"]])
})

test_that("the old shift_var name is rejected", {
  # simulate_mfBm forwards ... to hurst_fun, so without the guard a stale
  # shift_var would surface as an error raised inside the user's Hurst function.
  expect_error(simulate_mfBm(t = TG, hurst_fun = hurst_logistic, shift_var = IV),
               "shift_var")
  expect_error(sim_far(shift_var = IV), "shift_var")
  expect_error(sim_fma(shift_var = IV), "shift_var")
})

test_that("intercept_var reaches simulate_mfBm in both design branches", {
  t_com <- seq(0.2, 0.8, length.out = 20)
  set.seed(9); r0 <- sim_far()[["X"]]
  set.seed(9); r1 <- sim_far(intercept_var = IV)[["X"]]
  expect_false(isTRUE(all.equal(r0, r1)))

  set.seed(9); m0 <- sim_fma()[["X"]]
  set.seed(9); m1 <- sim_fma(intercept_var = IV)[["X"]]
  expect_false(isTRUE(all.equal(m0, m1)))

  common_far <- function(...) {
    simulate_far(N = 4L, design = "common", M_distribution = NULL, t_distribution = NULL,
                 t_common = t_com, hurst_fun = hurst_logistic, L = LC,
                 far_kernel = function(s, t) 9 / 4 * exp(-(t + 2 * s)^2),
                 far_mean = function(t) 4 * sin(1.5 * pi * t),
                 n_int_grid = 60L, n_burnin = 40L, remove_burnin = TRUE, ...)
  }
  set.seed(9); c0 <- common_far()[["X"]]
  set.seed(9); c1 <- common_far(intercept_var = IV)[["X"]]
  expect_false(isTRUE(all.equal(c0, c1)))

  common_fma <- function(...) {
    simulate_fma(N = 4L, design = "common", M_distribution = NULL, t_distribution = NULL,
                 t_common = t_com, hurst_fun = hurst_logistic, L = LC,
                 fma_kernel = function(s, t) 9 / 4 * exp(-(t + 2 * s)^2),
                 fma_mean = function(t) 4 * sin(1.5 * pi * t),
                 n_int_grid = 60L, n_burnin = 40L, remove_burnin = TRUE, ...)
  }
  set.seed(9); k0 <- common_fma()[["X"]]
  set.seed(9); k1 <- common_fma(intercept_var = IV)[["X"]]
  expect_false(isTRUE(all.equal(k0, k1)))
})

test_that("intercept_var is validated", {
  expect_error(simulate_mfBm(t = TG, intercept_var = -1), "non-negative")
  # No message pattern: `&&` on a length-2 operand warns on R 4.2 and errors from
  # R 4.3 on, so only the rejection itself is portable.
  expect_error(suppressWarnings(simulate_mfBm(t = TG, intercept_var = c(0.1, 0.2))))
  expect_error(sim_far(intercept_var = -1), "non-negative")
  expect_error(sim_fma(intercept_var = "0.05"), "non-negative")
})
