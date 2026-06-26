# Regression tests: lock the numerical output of every exported *_cpp function
# and its R wrapper to the committed references in tests/testthat/_refs/.
# These references were captured on the deterministic build; any future change
# that alters results beyond 1e-10 will fail here.

suppressMessages(library(data.table))
TOL <- 1e-10

test_that("kernel functions match references", {
  p <- ref_inputs()
  for (k in c("epanechnikov_kernel", "biweight_kernel", "triweight_kernel",
              "tricube_kernel", "triangular_kernel", "uniform_kernel")) {
    expect_equal(g(k)(p$uvec), read_ref(k), tolerance = TOL, info = k)
  }
})

test_that("matrix / grid helpers match references", {
  p <- ref_inputs()
  expect_equal(g("build_grid")(p$ss, p$ttac), read_ref("build_grid"), tolerance = TOL)
  expect_equal(g("get_upper_tri_couple")(c(.2,.5,.5,.8), c(.5,.2,.8,.5)),
               read_ref("get_upper_tri_couple"), tolerance = TOL)
  set.seed(1); M <- matrix(rnorm(16), 4, 4); M <- M %*% t(M)
  expect_equal(g("ensure_positive_definite")(M, 1e-6),
               read_ref("ensure_positive_definite"), tolerance = TOL)
  expect_equal(g("combine_matrices")(M, M, M, M), read_ref("combine_matrices"), tolerance = TOL)
  set.seed(2); zz <- rnorm(20)
  expect_equal(g("estimate_numerator_DD_single_t")(zz, 3L),
               read_ref("estimate_numerator_DD_single_t"), tolerance = TOL)
})

test_that("smoothing (NW) matches references", {
  p <- ref_inputs()
  expect_equal(g("estimate_nw_cpp")(p$one$X, p$one$tobs, p$tt, rep(0.1, length(p$tt)), "epanechnikov"),
               read_ref("estimate_nw_cpp"), tolerance = TOL)
  expect_equal(g("estimate_nw_bw_cpp")(p$one$X, p$one$tobs, p$bwg, "epanechnikov"),
               read_ref("estimate_nw_bw_cpp"), tolerance = TOL)
  expect_equal(estimate_nw(p$one$X, p$one$tobs, p$tt, h = 0.1),
               read_ref("estimate_nw_wrapper"), tolerance = TOL)
})

test_that("local regularity matches references", {
  p <- ref_inputs()
  expect_equal(g("estimate_locreg_cpp")(p$dtf, p$tt, FALSE, "epanechnikov", NULL, NULL),
               read_ref("estimate_locreg_cpp"), tolerance = TOL)
  expect_equal(estimate_locreg(p$dt, t = p$tt, center = FALSE),
               read_ref("estimate_locreg_wrapper"), tolerance = TOL)
})

test_that("mean estimation matches references", {
  p <- ref_inputs()
  expect_equal(g("estimate_mean_risk_cpp")(p$dtf, p$tt, p$bwg, "epanechnikov"),
               read_ref("estimate_mean_risk_cpp"), tolerance = TOL)
  expect_equal(g("estimate_mean_cpp")(p$dtf, p$tt, NULL, p$bwg, "epanechnikov"),
               read_ref("estimate_mean_cpp"), tolerance = TOL)
  expect_equal(estimate_mean(p$dt, t = p$tt, bw_grid = p$bwg),
               read_ref("estimate_mean_wrapper"), tolerance = TOL)
})

test_that("autocovariance matches references (both branches, both lags)", {
  p <- ref_inputs()
  expect_equal(g("estimate_autocov_risk_cpp")(p$dtf, p$ss, p$ttac, 0L, p$bwg, FALSE, TRUE, "epanechnikov"),
               read_ref("estimate_autocov_risk_cpp_lag0"), tolerance = TOL)
  expect_equal(g("estimate_autocov_risk_cpp")(p$dtf, p$ss, p$ttac, 0L, p$bwg, TRUE, TRUE, "epanechnikov"),
               read_ref("estimate_autocov_risk_cpp_samebw"), tolerance = TOL)
  grid <- g("build_grid")(p$vg, p$vg)
  expect_equal(g("estimate_autocov_risk_cpp")(p$dtf, grid[,1], grid[,2], 0L, NULL, FALSE, TRUE, "epanechnikov"),
               read_ref("autocov_risk_grid361_lag0"), tolerance = TOL)
  expect_equal(g("estimate_autocov_risk_cpp")(p$dtf, grid[,1], grid[,2], 1L, NULL, FALSE, TRUE, "epanechnikov"),
               read_ref("autocov_risk_grid361_lag1"), tolerance = TOL)
  expect_equal(g("estimate_autocov_cpp")(p$dtf, p$ss, p$ttac, 0L, NULL, NULL, p$bwg, FALSE, TRUE, TRUE, "epanechnikov"),
               read_ref("estimate_autocov_cpp_lag0"), tolerance = TOL)
  expect_equal(g("estimate_autocov_cpp")(p$dtf, p$ss, p$ttac, 1L, NULL, NULL, p$bwg, FALSE, TRUE, FALSE, "epanechnikov"),
               read_ref("estimate_autocov_cpp_lag1"), tolerance = TOL)
  expect_equal(estimate_autocov(p$dt, s = p$ss, t = p$ttac, lag = 1L, bw_grid = p$bwg),
               read_ref("estimate_autocov_wrapper"), tolerance = TOL)
})

test_that("covariance segment matches references", {
  p <- ref_inputs()
  expect_equal(g("estimate_cov_segment_risk_cpp")(p$dtf, p$tt, p$bwg, TRUE, "epanechnikov"),
               read_ref("estimate_cov_segment_risk_cpp"), tolerance = TOL)
  expect_equal(g("estimate_cov_segment_cpp")(p$dtf, p$tt, NULL, p$bwg, TRUE, "epanechnikov"),
               read_ref("estimate_cov_segment_cpp"), tolerance = TOL)
  expect_equal(estimate_cov_segment(p$dt, t = p$tt, bw_grid = p$bwg),
               read_ref("estimate_cov_segment_wrapper"), tolerance = TOL)
})

test_that("constants (sigma, moments, empirical autocov) match references", {
  p <- ref_inputs(); hN <- rep(0.1, p$N)
  expect_equal(g("estimate_sigma_cpp")(p$dtf, p$tt), read_ref("estimate_sigma_cpp"), tolerance = TOL)
  expect_equal(estimate_sigma(p$dt, idcol = "id_curve", t = p$tt),
               read_ref("estimate_sigma_wrapper"), tolerance = TOL)
  expect_equal(g("estimate_empirical_mom_cpp")(p$dtf, p$tt, hN, 2, 0, "epanechnikov"),
               read_ref("estimate_empirical_mom_cpp"), tolerance = TOL)
  expect_equal(estimate_empirical_mom(p$dt, idcol = "id_curve", t = p$tt, h = 0.1, mom_order = 2, center = FALSE),
               read_ref("estimate_empirical_mom_wrapper"), tolerance = TOL)
  expect_equal(g("estimate_empirical_autocov_cpp")(p$dtf, p$tt, hN, c(0,1,2), "epanechnikov"),
               read_ref("estimate_empirical_autocov_cpp"), tolerance = TOL)
  expect_equal(estimate_empirical_autocov(p$dt, idcol = "id_curve", t = p$tt, lag = c(0,1,2), h = 0.1),
               read_ref("estimate_empirical_autocov_wrapper"), tolerance = TOL)
  expect_equal(g("estimate_empirical_XsXt_autocov_cpp")(p$dtf, p$ss, p$ttac, 0L, c(0,1), hN, "epanechnikov", TRUE),
               read_ref("estimate_empirical_XsXt_autocov_cpp"), tolerance = TOL)
  expect_equal(estimate_empirical_XsXt_autocov(p$dt, idcol = "id_curve", s = p$ss, t = p$ttac, cross_lag = 0L, lag = c(0,1), h = 0.1),
               read_ref("estimate_empirical_XsXt_autocov_wrapper"), tolerance = TOL)
  expect_equal(g("estimate_numerator_dependence_term_DD_cpp")(p$dtf, p$tt, p$bwg, hN, 3L, "epanechnikov", TRUE),
               read_ref("estimate_numerator_dependence_term_DD_cpp"), tolerance = TOL)
})

test_that("full BLUP pipeline matches references", {
  p <- ref_inputs()
  cur <- g("estimate_curve_cpp")(p$dtf, p$tt, 2L, p$bwg, FALSE, TRUE, TRUE, "epanechnikov")
  ref <- read_ref("estimate_curve_cpp")
  expect_equal(sort(names(cur)), sort(names(ref)))
  for (nm in names(ref)) expect_equal(cur[[nm]], ref[[nm]], tolerance = TOL, info = nm)
  expect_equal(predict_curve(p$dt, t = p$tt, id_curve_to_predict = 2L, bw_grid = p$bwg),
               read_ref("predict_curve_wrapper"), tolerance = TOL)
})
