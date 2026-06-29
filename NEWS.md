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
