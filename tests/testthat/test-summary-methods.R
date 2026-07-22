# S3 classing + text summaries for the adaptive-estimator outputs.
# Uses the deterministic packaged-data slice (helper-fixtures.R).

tt <- c(0.25, 0.5, 0.75)
ss <- c(0.2, 0.4, 0.6)

est_list <- function() {
  d <- fixture_data_far(12L)
  list(
    mean_est         = estimate_mean(d, t = tt),
    mean_risk        = estimate_mean_risk(d, t = tt),
    locreg_est       = estimate_locreg(d, t = tt),
    autocov_est      = estimate_autocov(d, s = ss, t = tt, lag = 1),
    autocov_risk     = estimate_autocov_risk(d, s = ss, t = tt, lag = 1),
    cov_segment_est  = estimate_cov_segment(d, t = tt),
    cov_segment_risk = estimate_cov_segment_risk(d, t = tt)
  )
}

test_that("adaptive estimators carry their S3 class and stay data.tables", {
  ests <- est_list()
  for (nm in names(ests)) {
    x <- ests[[nm]]
    expect_s3_class(x, nm)                       # per-estimator subclass
    expect_s3_class(x, "adaptiveFTS_est")        # shared parent
    expect_true(data.table::is.data.table(x), info = nm)
    expect_true(is.list(attr(x, "adaptive_meta")), info = nm)
  }
})

test_that("summary methods print and return the object invisibly", {
  ests <- est_list()
  for (nm in names(ests)) {
    x <- ests[[nm]]
    expect_output(summary(x))                    # prints a report
    expect_invisible(summary(x))                 # returns invisibly
    expect_identical(summary(x), x)              # returns its argument
  }
})

test_that("summary content reflects the estimator's metadata", {
  ests <- est_list()
  expect_output(summary(ests$mean_est), "Adaptive mean function estimate")
  expect_output(summary(ests$autocov_est), "lag = 1")
  expect_output(summary(ests$mean_risk), "Risk-minimising bandwidth")
  expect_output(summary(ests$locreg_est), "local regularity")
})

test_that("the parent-class fallback summary works on a bare adaptiveFTS_est", {
  dt <- data.table::data.table(t = tt, value = tt)
  x <- adaptiveFTS:::.as_adaptive_est(dt, "not_a_real_subclass")
  expect_output(summary(x), "adaptiveFTS estimator object")
  expect_invisible(summary(x))
})
