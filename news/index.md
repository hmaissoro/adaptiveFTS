# Changelog

## adaptiveFTS 0.2.0

### New features

- Design-weighted, Tikhonov-regularised adaptive functional BLUP with an
  [`lm()`](https://rdrr.io/r/stats/lm.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)-style
  interface:
  - [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)
    estimates every prediction-point-independent component and caches
    the adaptive bandwidths;
    [`predict()`](https://rdrr.io/r/stats/predict.html) (method
    [`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md))
    evaluates the predictor and, for `horizon > 1`, returns every
    intermediate multi-step-ahead prediction;
    [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
    is a one-call wrapper.
  - The Tikhonov parameter is selected automatically by default
    (`tikhonov = NULL`), like `optbw` for `bw_grid`:
    [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)/[`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
    cross-validate over `tikhonov_grid` and return the selection as
    `tikhonov_cv`. Pass a numeric `tikhonov` to skip selection.
    `select_tikhonov_parameter(method = "cv")` exposes the selection
    directly (holdout under the common design, rolling origin under the
    independent design).
  - [`summary()`](https://rdrr.io/r/base/summary.html) methods for
    `blup_fit` and `blup` objects.
  - The numerical core runs in C++
    ([`blup_fit_cpp()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit_cpp.md),
    [`blup_predict_cpp()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_predict_cpp.md)).
- Design-density estimation for the independent-design weights:
  - [`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md)
    — leave-one-out Parzen–Rosenblatt estimator with a least-squares
    cross-validation or fixed-bandwidth path.
  - [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md)
    — selects the density bandwidth on a subset of curves (median of the
    per-curve LSCV optima).
- The adaptive estimators now return classed `data.table`s (see
  [`?adaptiveFTS_est`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS_est.md))
  so they gain [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods:
  - [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
    [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
    [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
    [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)
    and their `*_risk()` counterparts print a compact, design-aware
    [`summary()`](https://rdrr.io/r/base/summary.html) and draw a
    diagnostic `ggplot2` plot (mean/segment curves, a regularity panel,
    a covariance surface, and risk-vs-bandwidth curves marking the
    minimiser). The returned objects remain genuine `data.table`s, so
    existing code is unaffected. `ggplot2` stays a suggested dependency.
- Descriptive statistics:
  - [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md)
    — adaptive functional autocorrelation function (FACF),
    `rho_l = ||Gamma_l|| / integral Gamma_0(t, t) dt`, built on
    [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
    and handling both the common and independent designs. Returns a
    classed `fts_acf` object with
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
    [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
    (an ACF-style lag plot).
- [`save_plot_tikz()`](https://hmaissoro.github.io/adaptiveFTS/reference/save_plot_tikz.md)
  — export a `ggplot` (or any printable figure) to a standalone
  TikZ/LaTeX `.tex` file (optionally compiled to PDF), so the helper no
  longer has to be copied between projects. `tikzDevice` is a new
  suggested dependency.

### Deprecations

- [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)
  is deprecated in favour of
  [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  and
  [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md).
  It still works (with a warning) this release and will be removed in a
  future version. Note the new functions compute a different quantity
  (the design-weighted adaptive BLUP), so results are not
  interchangeable.

## adaptiveFTS 0.1.1

First CRAN-targeted release.

### Features

- Local regularity parameter estimation
  ([`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)).
- Adaptive mean function estimation
  ([`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
  [`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)).
- Adaptive autocovariance / covariance function estimation
  ([`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
  [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md),
  [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)).
- Adaptive Best Linear Unbiased Predictor
  ([`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)).
- Nadaraya–Watson smoothing and bandwidth selection, kernels, simulation
  of functional time series (FAR/FMA, fBm/mfBm) and Rubin–Panaretos
  estimators.
- Numerical core implemented in C++ via Rcpp/RcppArmadillo.

### Reproducibility and performance

- Local regularity estimation is now fully reproducible: the
  tie-breaking jitter that previously used Armadillo’s RNG (uncontrolled
  by [`set.seed()`](https://rdrr.io/r/base/Random.html)) is replaced by
  a deterministic offset.
- Substantial speed-ups of the C++ estimators (autocovariance risk and
  estimation, mean and covariance-segment risk) with bit-identical
  results, validated against committed regression references.

### Infrastructure

- Removed the `caret`, `fastmatrix` and `parallel` dependencies in
  favour of base R and `data.table`.
- Added a `testthat` (edition 3) test suite covering every exported
  function plus numerical regression tests.
- Added continuous integration (R-CMD-check on Linux/macOS/Windows, test
  coverage, lint, and a pkgdown site).
