# The local regularity knobs (`presmooth_bw`, `Delta`, `presmooth_bw_grid`,
# `presmooth_nsubset`) are consumed inside C++, several layers below the R entry
# points. These tests check that they actually travel all the way down: a
# silently dropped argument would leave the output identical to the default.
#
# What is asserted, and why:
#   * The `Ht` / `Lt2` / `locreg_bw` columns and the risk surfaces are the direct
#     product of the local regularity step, so they must move whenever a knob
#     changes. These are the real propagation tests.
#   * The point estimates (`muhat`, `autocov`, `cov_segment_hat`) are evaluated
#     at the risk-MINIMISING bandwidth. On a fixture this small the argmin often
#     sits on the boundary of `bw_grid` and does not move even though the whole
#     risk surface did, so asserting that they change would test an artifact of
#     the fixture rather than the plumbing. They are asserted only where the
#     argmin is genuinely interior (the BLUP block, on a fine grid).
#
# The complementary guarantee, that the defaults still reproduce the committed
# references bit-for-bit, is covered by test-regression-refs.R and
# test-blup-refs.R, which are deliberately left untouched.

TOL <- 1e-10

# TRUE when two results are numerically indistinguishable.
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = TOL))

test_that(".check_locreg_args rejects malformed input", {
  expect_error(.check_locreg_args(presmooth_bw = -0.1), "presmooth_bw")
  expect_error(.check_locreg_args(presmooth_bw = "a"), "presmooth_bw")
  expect_error(.check_locreg_args(Delta = 0), "Delta")
  expect_error(.check_locreg_args(Delta = 1.5), "Delta")
  expect_error(.check_locreg_args(Delta = c(0.1, 0.2)), "Delta")
  expect_error(.check_locreg_args(presmooth_bw_grid = 0.1), "presmooth_bw_grid")
  expect_error(.check_locreg_args(presmooth_bw_grid = c(-1, 0.2)), "presmooth_bw_grid")
  expect_error(.check_locreg_args(presmooth_nsubset = 0), "presmooth_nsubset")
  expect_error(.check_locreg_args(presmooth_nsubset = 2.5), "presmooth_nsubset")
  expect_silent(.check_locreg_args())
  expect_silent(.check_locreg_args(presmooth_bw = 0.05, Delta = 0.1,
                                   presmooth_bw_grid = c(0.02, 0.2),
                                   presmooth_nsubset = 5L))
})

test_that("estimate_locreg honours the new cross-validation knobs", {
  dt <- fixture_data_far(12)
  tt <- fixture_tgrid()

  # A degenerate grid pins the presmoothing bandwidth to that single value.
  expect_equal(unique(estimate_locreg(dt, t = tt,
                                      presmooth_bw_grid = c(0.08, 0.08))$locreg_bw),
               0.08, tolerance = TOL)

  # Restricting the cross-validation to fewer curves changes the selected
  # bandwidth, and with it the regularity estimates.
  lr_all <- estimate_locreg(dt, t = tt)
  lr_few <- estimate_locreg(dt, t = tt, presmooth_nsubset = 3L)
  expect_false(same(lr_all$locreg_bw, lr_few$locreg_bw))
  expect_false(same(lr_all$Ht, lr_few$Ht))

  # Supplying presmooth_bw short-circuits the cross-validation entirely, so the
  # two search knobs must become inert.
  expect_equal(estimate_locreg(dt, t = tt, presmooth_bw = 0.05,
                               presmooth_nsubset = 3L,
                               presmooth_bw_grid = c(0.01, 0.02)),
               estimate_locreg(dt, t = tt, presmooth_bw = 0.05),
               tolerance = TOL)

  # Delta is used and reported back verbatim, and is now validated.
  expect_equal(unique(estimate_locreg(dt, t = tt, Delta = 0.15)$Delta), 0.15,
               tolerance = TOL)
  expect_false(same(estimate_locreg(dt, t = tt, Delta = 0.15)$Ht, lr_all$Ht))
  expect_error(estimate_locreg(dt, t = tt, Delta = 2), "Delta")
})

test_that("presmooth_bw reaches the mean risk and its downstream moments", {
  dt <- fixture_data_far(12)
  tt <- c(0.25, 0.5, 0.75)
  bwg <- seq(0.05, 0.2, length.out = 6)

  risk_def <- estimate_mean_risk(dt, t = tt, bw_grid = bwg)
  risk_bw  <- estimate_mean_risk(dt, t = tt, bw_grid = bwg, presmooth_bw = 0.05)

  # The value is used for the regularity step ...
  expect_equal(unique(risk_bw$locreg_bw), 0.05, tolerance = TOL)
  expect_false(same(unique(risk_def$locreg_bw), 0.05))
  # ... it changes the estimated regularity ...
  expect_false(same(risk_def$Ht, risk_bw$Ht))
  # ... and, because the same bandwidth also feeds the empirical autocovariance
  # behind the dependence term, the whole risk surface moves.
  expect_false(same(risk_def$mean_risk, risk_bw$mean_risk))

  # estimate_mean takes the same risk path when bw is NULL, and reports the
  # regularity it used at the selected bandwidth.
  m_def <- estimate_mean(dt, t = tt, bw_grid = bwg)
  m_bw  <- estimate_mean(dt, t = tt, bw_grid = bwg, presmooth_bw = 0.05)
  expect_false(same(m_def$Ht, m_bw$Ht))
  expect_false(same(m_def$Lt2, m_bw$Lt2))

  # Supplying bw skips the risk path, so the knobs cannot bite there.
  expect_equal(estimate_mean(dt, t = tt, bw = rep(0.1, 3), presmooth_bw = 0.05),
               estimate_mean(dt, t = tt, bw = rep(0.1, 3)), tolerance = TOL)
})

test_that("the knobs reach the autocovariance and covariance-segment risks", {
  dt <- fixture_data_far(12)
  tt <- c(0.25, 0.5, 0.75)
  ss <- c(0.2, 0.4, 0.8)
  bwg <- seq(0.05, 0.2, length.out = 6)

  ac_def <- estimate_autocov_risk(dt, s = ss, t = tt, lag = 1, bw_grid = bwg)
  expect_false(same(
    ac_def$Ht,
    estimate_autocov_risk(dt, s = ss, t = tt, lag = 1, bw_grid = bwg,
                          Delta = 0.15)$Ht))
  expect_false(same(
    ac_def$autocov_risk,
    estimate_autocov_risk(dt, s = ss, t = tt, lag = 1, bw_grid = bwg,
                          presmooth_bw = 0.05)$autocov_risk))

  # estimate_autocov reports the regularity of the selected bandwidth pair.
  ac <- estimate_autocov(dt, s = ss, t = tt, lag = 1, bw_grid = bwg)
  expect_false(same(ac$Ht, estimate_autocov(dt, s = ss, t = tt, lag = 1,
                                            bw_grid = bwg, Delta = 0.15)$Ht))
  expect_false(same(ac$Ht, estimate_autocov(dt, s = ss, t = tt, lag = 1,
                                            bw_grid = bwg,
                                            presmooth_bw = 0.05)$Ht))

  cs_def <- estimate_cov_segment_risk(dt, t = tt, bw_grid = bwg)
  cs_bw  <- estimate_cov_segment_risk(dt, t = tt, bw_grid = bwg, presmooth_bw = 0.05)
  expect_equal(unique(cs_bw$locreg_bw), 0.05, tolerance = TOL)
  expect_false(same(cs_def$Ht, cs_bw$Ht))
  expect_false(same(cs_def$cov_segment_risk, cs_bw$cov_segment_risk))

  expect_false(same(estimate_cov_segment(dt, t = tt, bw_grid = bwg)$Ht,
                    estimate_cov_segment(dt, t = tt, bw_grid = bwg,
                                         presmooth_bw = 0.05)$Ht))
})

test_that("estimate_facf forwards the knobs to every lag", {
  dt <- fixture_data_far(12)
  f_def <- estimate_facf(dt, lag.max = 2L, n_grid = 5L)
  expect_false(same(f_def$facf,
                    estimate_facf(dt, lag.max = 2L, n_grid = 5L,
                                  presmooth_bw = 0.05)$facf))
  expect_false(same(f_def$facf,
                    estimate_facf(dt, lag.max = 2L, n_grid = 5L,
                                  Delta = 0.15)$facf))
})

test_that("blup_fit forwards the knobs to the three risk minimisations", {
  dt <- fixture_data_far(12)
  # A fine grid, so the risk-minimising bandwidths are interior and the change
  # propagates all the way to the prediction.
  bwg <- seq(0.02, 0.3, length.out = 40)

  fit_def <- blup_fit(dt, tikhonov = 1e-6, bw_grid = bwg)
  fit_bw  <- blup_fit(dt, tikhonov = 1e-6, bw_grid = bwg, presmooth_bw = 0.05)

  # One risk minimisation per cached bandwidth set; all three must move.
  expect_false(same(fit_def$opt_mean, fit_bw$opt_mean))
  expect_false(same(fit_def$opt_cov, fit_bw$opt_cov))
  expect_false(same(fit_def$opt_autocov, fit_bw$opt_autocov))

  # ... and so must the quantities built on them.
  expect_false(same(fit_def$muhat_Tn0, fit_bw$muhat_Tn0))
  expect_false(same(fit_def$c0hat, fit_bw$c0hat))

  p_def <- predict(fit_def, t = c(0.3, 0.6))
  p_bw  <- predict(fit_bw,  t = c(0.3, 0.6))
  expect_false(same(p_def$prediction, p_bw$prediction))

  # blup() is a thin wrapper and must agree with the two-step form.
  b <- blup(dt, t = c(0.3, 0.6), tikhonov = 1e-6, bw_grid = bwg,
            presmooth_bw = 0.05)
  expect_equal(b$prediction, p_bw, tolerance = TOL)

  expect_error(blup_fit(dt, tikhonov = 1e-6, bw_grid = bwg, Delta = 3), "Delta")
})

test_that("results with the knobs set are deterministic", {
  dt <- fixture_data_far(12)
  tt <- c(0.25, 0.5, 0.75)
  a <- estimate_mean_risk(dt, t = tt, presmooth_bw = 0.05, Delta = 0.12)
  b <- estimate_mean_risk(dt, t = tt, presmooth_bw = 0.05, Delta = 0.12)
  expect_equal(a, b, tolerance = TOL)
})
