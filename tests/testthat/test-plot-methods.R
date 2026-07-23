# ggplot2 autoplot()/plot() methods for the adaptive-estimator outputs.
# ggplot2 is a suggested dependency; skip the whole file when it is absent.

skip_if_not_installed("ggplot2")

tt <- c(0.25, 0.5, 0.75)
ss <- c(0.2, 0.4, 0.6)

test_that("autoplot returns ggplot objects for every classed estimator", {
  d <- fixture_data_far(12L)
  ests <- list(
    estimate_mean(d, t = tt),
    estimate_mean_risk(d, t = tt),
    estimate_locreg(d, t = tt),
    estimate_autocov(d, s = ss, t = tt, lag = 1),   # off-diagonal -> surface
    estimate_autocov(d, s = tt, t = tt, lag = 0),   # diagonal -> curve
    estimate_autocov_risk(d, s = ss, t = tt, lag = 1),
    estimate_cov_segment(d, t = tt),
    estimate_cov_segment_risk(d, t = tt)
  )
  for (x in ests) {
    p <- ggplot2::autoplot(x)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot() draws and invisibly returns the ggplot", {
  d <- fixture_data_far(12L)
  x <- estimate_mean(d, t = tt)
  # Draw to a throw-away device so no Rplots.pdf is created.
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(x))
  expect_s3_class(plot(x), "ggplot")
})

test_that("locreg autoplot honours the `which` argument", {
  d <- fixture_data_far(12L)
  x <- estimate_locreg(d, t = tt)
  expect_s3_class(ggplot2::autoplot(x, which = "Ht"), "ggplot")
  expect_s3_class(ggplot2::autoplot(x, which = c("Ht", "Lt2")), "ggplot")
})
