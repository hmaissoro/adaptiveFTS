# adaptiveFTS 0.2.0

## New features

* Design-weighted, Tikhonov-regularised adaptive functional BLUP with an
  `lm()`/`predict()`-style interface:
  * `blup_fit()` estimates every prediction-point-independent component and
    caches the adaptive bandwidths; `predict()` (method `predict.blup_fit()`)
    evaluates the one-step-ahead predictor and loops for h-step-ahead
    prediction; `blup()` is a one-call wrapper.
  * `select_tikhonov_parameter(method = "cv")` selects the Tikhonov parameter by
    one-step-ahead cross-validation (holdout under the common design, rolling
    origin under the independent design).
  * The numerical core runs in C++ (`blup_fit_cpp()`, `blup_predict_cpp()`).
* Design-density estimation for the independent-design weights:
  * `estimate_density()` — leave-one-out Parzen–Rosenblatt estimator with a
    least-squares cross-validation or fixed-bandwidth path.
  * `get_density_optimal_bw()` — selects the density bandwidth on a subset of
    curves (median of the per-curve LSCV optima).

## Deprecations

* `predict_curve()` is deprecated in favour of `blup_fit()`/`predict()` and
  `blup()`. It still works (with a warning) this release and will be removed in
  a future version. Note the new functions compute a different quantity (the
  design-weighted adaptive BLUP), so results are not interchangeable.

# adaptiveFTS 0.1.1

First CRAN-targeted release.

## Features

* Local regularity parameter estimation (`estimate_locreg()`).
* Adaptive mean function estimation (`estimate_mean()`, `estimate_mean_risk()`).
* Adaptive autocovariance / covariance function estimation
  (`estimate_autocov()`, `estimate_autocov_risk()`, `estimate_cov_segment()`).
* Adaptive Best Linear Unbiased Predictor (`predict_curve()`).
* Nadaraya–Watson smoothing and bandwidth selection, kernels, simulation of
  functional time series (FAR/FMA, fBm/mfBm) and Rubin–Panaretos estimators.
* Numerical core implemented in C++ via Rcpp/RcppArmadillo.

## Reproducibility and performance

* Local regularity estimation is now fully reproducible: the tie-breaking jitter
  that previously used Armadillo's RNG (uncontrolled by `set.seed()`) is replaced
  by a deterministic offset.
* Substantial speed-ups of the C++ estimators (autocovariance risk and
  estimation, mean and covariance-segment risk) with bit-identical results,
  validated against committed regression references.

## Infrastructure

* Removed the `caret`, `fastmatrix` and `parallel` dependencies in favour of
  base R and `data.table`.
* Added a `testthat` (edition 3) test suite covering every exported function plus
  numerical regression tests.
* Added continuous integration (R-CMD-check on Linux/macOS/Windows, test
  coverage, lint, and a pkgdown site).
