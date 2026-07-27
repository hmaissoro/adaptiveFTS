# Package index

## Data formatting

- [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  : Convert Raw Curve Observations to the Package Data Format

## Simulation

- [`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
  : Simulate a Functional Autoregressive Process of Order One
- [`simulate_fma()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_fma.md)
  : Simulate a Functional Moving Average Process of Order One
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
  : Estimate the Local Regularity Parameters

## Mean function

- [`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md)
  : Estimate the Mean Function
- [`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)
  : Estimate the Risk of the Mean Function Estimator
- [`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md)
  : Estimate the Mean Function by the Rubìn-Panaretos Method
- [`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md)
  : Select the Bandwidth of the Rubìn-Panaretos Mean Estimator

## Autocovariance function

- [`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md)
  : Estimate the Covariance or Autocovariance Function
- [`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
  : Estimate the Risk of the Autocovariance Function Estimator
- [`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md)
  : Estimate the Autocovariance Function by the Rubìn-Panaretos Method
- [`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md)
  : Select the Bandwidth of the Rubìn-Panaretos Autocovariance Estimator

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

## Descriptive statistics

- [`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md)
  : Estimate the functional autocorrelation function (FACF)

## Summaries and plots

- [`summary(`*`<locreg_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<mean_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<cov_segment_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<autocov_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<mean_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<cov_segment_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<autocov_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<adaptiveFTS_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  [`summary(`*`<fts_acf>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-summary.md)
  : Summarise an adaptive functional time series estimator
- [`autoplot.mean_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<mean_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.locreg_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<locreg_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.cov_segment_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<cov_segment_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.autocov_est()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<autocov_est>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<mean_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.cov_segment_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<cov_segment_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<autocov_risk>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`autoplot.fts_acf()`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  [`plot(`*`<fts_acf>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md)
  : Plot an adaptive functional time series estimator

## Export

- [`save_plot_tikz()`](https://hmaissoro.github.io/adaptiveFTS/reference/save_plot_tikz.md)
  : Save a plot as a (standalone) TikZ/LaTeX figure

## Prediction (BLUP)

- [`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md)
  : Fit the adaptive functional BLUP
- [`predict(`*`<blup_fit>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md)
  : Predict with an adaptive functional BLUP fit
- [`blup()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup.md)
  : Fit and predict the adaptive functional BLUP in one call
- [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md)
  : Select the Tikhonov regularisation parameter
- [`summary(`*`<blup_fit>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/summary.blup_fit.md)
  : Summarise an adaptive functional BLUP fit
- [`summary(`*`<blup>`*`)`](https://hmaissoro.github.io/adaptiveFTS/reference/summary.blup.md)
  : Summarise an adaptive functional BLUP prediction
- [`predict_curve()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict_curve.md)
  : Curve prediction using the Best Linear Unbiased Predictor (BLUP).

## Design density

- [`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md)
  : Leave-one-out Parzen-Rosenblatt density estimator
- [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md)
  : Select the design-density bandwidth on a subset of curves

## Real-data helpers

- [`get_real_data_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_real_data_mean.md)
  : Mean function learned from the voltage curves of the electricity
- [`get_real_data_far_kernel()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_real_data_far_kernel.md)
  : FAR kernel learned from the voltage curves of the electricity

## Data

- [`data_far`](https://hmaissoro.github.io/adaptiveFTS/reference/data_far.md)
  : Sample for Functional Autoregressive Process of Order 1

## Internal helpers

- [`.Spq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Spq_fun.md)
  : Weight Sum \\S\_{pq}^{(\ell)}\\ of the Rubìn-Panaretos
  Autocovariance Estimator
- [`.Qpq_fun()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Qpq_fun.md)
  : Weighted Cross-Product \\Q\_{pq}^{(\ell)}\\ of the Rubìn-Panaretos
  Estimator
- [`.constant_d()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-constant_d.md)
  : Constant D(x,y) function
- [`.covariance_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-covariance_mfBm.md)
  : Covariance matrix of the multi-fractional Brownian Motion
- [`.random_design()`](https://hmaissoro.github.io/adaptiveFTS/reference/dot-random_design.md)
  : Generate a random design of the place where the observation is made.
