# Package index

## Data formatting

- [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  :

  Convert Data to a `data.table` Format

## Simulation

- [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
  : Functional Autoregressive process of order 1 (FAR(1)) simulation
- [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md)
  : Functional Moving Average process of order 1 (FMA(1)) simulation
- [`simulate_fBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fBm.md)
  : Draw a fractional Brownian motion sample path.
- [`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)
  : Draw a multifractional Brownian motion sample path.
- [`hurst_arctan()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md)
  : Arctan Hurst function
- [`hurst_linear()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md)
  : Linear Hurst function
- [`hurst_logistic()`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md)
  : Logistic Hurst function

## Kernels and smoothing

- [`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md)
  : Nadaraya-Watson Kernel Estimator
- [`estimate_nw_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw_bw.md)
  : Nadaraya-Watson Bandwidth Selection using Cross-Validation
- [`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md)
  : Estimate Optimal Bandwidth for Nadaraya-Watson Estimator on a Subset
  of Curves
- [`biweight()`](https://hmaissoro.github.io/adaptiveFTS/reference/biweight.md)
  : Biweight kernel function
- [`epanechnikov()`](https://hmaissoro.github.io/adaptiveFTS/reference/epanechnikov.md)
  : Epanechnikov kernel function
- [`triangular()`](https://hmaissoro.github.io/adaptiveFTS/reference/triangular.md)
  : Triangular kernel function
- [`tricube()`](https://hmaissoro.github.io/adaptiveFTS/reference/tricube.md)
  : Tricube kernel function
- [`triweight()`](https://hmaissoro.github.io/adaptiveFTS/reference/triweight.md)
  : Triweight kernel function
- [`uniform()`](https://hmaissoro.github.io/adaptiveFTS/reference/uniform.md)
  : Uniform kernel function

## Local regularity

- [`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md)
  : Local Regularity Parameters Estimation

## Mean function

- [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md)
  : Estimate Mean Function
- [`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)
  : Estimate the Risk Function of the Mean Function
- [`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md)
  : Estimate mean function using Rubín and Panaretos (2020) method.
- [`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md)
  : Bandwidth estimation using cross-validation for the Rubín and
  Panaretos (2020) mean function estimator.

## Autocovariance function

- [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
  : Estimate the Covariance or Autocovariance Function
- [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
  : Estimate the Risk of the Covariance or Autocovariance Function
- [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md)
  : Estimate lag-\\\ell\\ (\\\ell \leq 0\\) autocovariance function
  using Rubín and Panaretos (2020) method
- [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md)
  : Bandwidth estimation using cross-validation for the Rubín and
  Panaretos (2020) autocovariance function estimator.

## Covariance segment

- [`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)
  : Estimate Covariance Segment Function for Functional Data
- [`estimate_cov_segment_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md)
  : Estimate the Risk of the Covariance Segment Function

## Constants and empirical moments

- [`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md)
  : Estimate the the standard deviation of the observation error
- [`estimate_empirical_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md)
  : Estimate Empirical Autocovariance Function
- [`estimate_empirical_mom()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_mom.md)
  : Estimate empirical \\p\\-th order moment of \\X(t)\\.
- [`estimate_empirical_XsXt_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_XsXt_autocov.md)
  : Estimate Empirical \\X_0(s)X\_{\ell}(t)\\ Autocovariance Function
  for \\\ell\\ = 0, 1, ...

## Prediction (BLUP)

- [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)
  : Curve prediction using the Best Linear Unbiased Predictor (BLUP).

## Real-data helpers

- [`get_real_data_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_real_data_mean.md)
  : Mean function learned from the voltage curves of the electricity
- [`get_real_data_far_kenel()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_real_data_far_kenel.md)
  : FAR kernel learned from the voltage curves of the electricity

## Data

- [`data_far`](https://hmaissoro.github.io/adaptiveFTS/reference/data_far.md)
  : Sample for Functional Autoregressive Process of Order 1

## Internal helpers

- [`.Spq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Spq_fun.md)
  : \\S\_{pq}^{(\ell)}\\, (\\\ell \leq 0\\) function. See Rubín and
  Panaretos (2020) Equation (B.7)
- [`.Qpq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Qpq_fun.md)
  : \\Q\_{pq}^{(\ell)}\\, (\\\ell \leq 0\\) function. See Rubín and
  Panaretos (2020) Equation (B.7)
- [`.constant_d()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-constant_d.md)
  : Constant D(x,y) function
- [`.covariance_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-covariance_mfBm.md)
  : Covariance matrix of the multi-fractional Brownian Motion
- [`.random_design()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-random_design.md)
  : Generate a random design of the place where the observation is made.
