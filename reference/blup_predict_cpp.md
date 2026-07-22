# Predict with the adaptive functional BLUP (C++ core)

Evaluates the adaptive BLUP at the prediction points, looping for
h-step-ahead prediction. Called by
[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md);
not intended to be used directly.

## Usage

``` r
blup_predict_cpp(
  data,
  opt_mean,
  opt_cov,
  opt_autocov,
  Tn0,
  muhat_Tn0,
  V,
  root_D,
  Yn0,
  density_bw,
  is_common,
  homoscedastic,
  tikhonov,
  t,
  horizon,
  kernel_name
)
```

## Arguments

- data:

  A DataFrame with columns `id_curve`, `tobs`, `X`.

- opt_mean, opt_cov, opt_autocov:

  Cached adaptive-bandwidth matrices.

- Tn0:

  Conditioning-curve design points.

- muhat_Tn0:

  Mean at the conditioning-curve design points.

- V:

  Regularised variance matrix of the fit.

- root_D:

  Square-root design-weight matrix of the fit.

- Yn0:

  Conditioning-curve values (possibly overriding the fit's).

- density_bw:

  Fixed design-density bandwidth.

- is_common:

  Whether the design is common across curves.

- homoscedastic:

  Whether to use a constant noise variance.

- tikhonov:

  Tikhonov regularisation parameter.

- t:

  Prediction points.

- horizon:

  Prediction horizon (steps ahead).

- kernel_name:

  Kernel name.

## Value

A matrix with columns `t`, `muhat`, `prediction`.
