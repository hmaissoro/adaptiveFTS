# adaptiveFTS 0.2.0

## Breaking changes

* The Hölder-constant columns are renamed from `Lt`/`Ls` to `Lt2`/`Ls2` in the
  output of `estimate_locreg()`, `estimate_mean()`/`estimate_mean_risk()`,
  `estimate_autocov()`/`estimate_autocov_risk()` and
  `estimate_cov_segment()`/`estimate_cov_segment_risk()`, because the estimate is
  the squared constant (L_t^2 / L_s^2), not L_t / L_s. Update any code that
  referred to the `Lt`/`Ls` columns.

## New features

* Design-weighted, Tikhonov-regularised adaptive functional BLUP with an
  `lm()`/`predict()`-style interface:
  * `blup_fit()` estimates every prediction-point-independent component and
    caches the adaptive bandwidths; `predict()` (method `predict.blup_fit()`)
    evaluates the predictor and, for `horizon > 1`, returns every intermediate
    multi-step-ahead prediction; `blup()` is a one-call wrapper.
  * The Tikhonov parameter is selected automatically by default (`tikhonov =
    NULL`), like `optbw` for `bw_grid`: `blup_fit()`/`blup()` cross-validate over
    `tikhonov_grid` and return the selection as `tikhonov_cv`. Pass a numeric
    `tikhonov` to skip selection. `select_tikhonov_parameter(method = "cv")`
    exposes the selection directly (holdout under the common design, rolling
    origin under the independent design).
  * `summary()` methods for `blup_fit` and `blup` objects.
  * The numerical core runs in C++ (`blup_fit_cpp()`, `blup_predict_cpp()`).
* Design-density estimation for the independent-design weights:
  * `estimate_density()` — leave-one-out Parzen–Rosenblatt estimator with a
    least-squares cross-validation or fixed-bandwidth path.
  * `get_density_optimal_bw()` — selects the density bandwidth on a subset of
    curves (median of the per-curve LSCV optima).
* The adaptive estimators now return classed `data.table`s (see
  `?adaptiveFTS_est`) so they gain `summary()`, `plot()` and
  `ggplot2::autoplot()` methods:
  * `estimate_mean()`, `estimate_locreg()`, `estimate_autocov()`,
    `estimate_cov_segment()` and their `*_risk()` counterparts print a compact,
    design-aware `summary()` and draw a diagnostic `ggplot2` plot (mean/segment
    curves, a regularity panel, a covariance surface, and risk-vs-bandwidth
    curves marking the minimiser). The returned objects remain genuine
    `data.table`s, so existing code is unaffected. `ggplot2` stays a suggested
    dependency.
* Descriptive statistics:
  * `estimate_facf()` — adaptive functional autocorrelation function (FACF),
    `rho_l = ||Gamma_l|| / integral Gamma_0(t, t) dt`, built on
    `estimate_autocov()` and handling both the common and independent designs.
    Returns a classed `fts_acf` object with `summary()`, `plot()` and
    `ggplot2::autoplot()` (an ACF-style lag plot).
* `save_plot_tikz()` — export a `ggplot` (or any printable figure) to a
  standalone TikZ/LaTeX `.tex` file (optionally compiled to PDF), so the helper
  no longer has to be copied between projects. `tikzDevice` is a new suggested
  dependency.

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
