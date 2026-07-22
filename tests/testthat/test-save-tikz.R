# Reusable TikZ export. tikzDevice needs a working LaTeX engine to compute
# string-width metrics; skip when it is absent or non-functional here.

skip_if_not_installed("tikzDevice")
skip_if_not_installed("ggplot2")

make_plot <- function() {
  d <- fixture_data_far(12L)
  ggplot2::autoplot(estimate_mean(d, t = c(0.25, 0.5, 0.75)))
}

# Can tikzDevice actually render this figure here? Renders via tikzDevice
# directly (not through save_plot_tikz), so only genuine environment/metric
# failures cause a skip -- a real bug in save_plot_tikz still surfaces.
tikz_render_ok <- function(p) {
  ff <- tempfile(fileext = ".tex")
  tryCatch({ tikzDevice::tikz(ff); print(p); TRUE },
           error = function(e) FALSE,
           finally = if (grDevices::dev.cur() > 1L) grDevices::dev.off())
}

test_that("save_plot_tikz writes a non-empty .tex file", {
  p <- make_plot()
  skip_if(!tikz_render_ok(p), "tikzDevice cannot render here (LaTeX metrics)")
  f <- tempfile(fileext = ".tex")
  out <- save_plot_tikz(p, file = f, width = 5, height = 4, compile = FALSE)
  expect_identical(out, f)
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
})

test_that("save_plot_tikz appends .tex and creates parent directories", {
  p <- make_plot()
  skip_if(!tikz_render_ok(p), "tikzDevice cannot render here (LaTeX metrics)")
  base <- file.path(tempdir(), paste0("tkz_", as.integer(runif(1, 1, 1e6))), "fig")
  out <- save_plot_tikz(p, file = base, compile = FALSE)
  expect_true(grepl("\\.tex$", out))
  expect_true(file.exists(out))
})
