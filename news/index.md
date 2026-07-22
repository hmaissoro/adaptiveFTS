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
