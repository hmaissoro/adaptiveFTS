# Fit the adaptive functional BLUP (C++ core)

Estimates every prediction-point-independent component of the adaptive
Best Linear Unbiased Predictor and caches the adaptive bandwidths.
Called by the R function
[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md);
not intended to be used directly.

## Usage

``` r
blup_fit_cpp(
  data,
  id_lag,
  bw_grid,
  rho,
  homoscedastic,
  tikhonov,
  n_subgrid_bw,
  kernel_name
)
```

## Arguments

- data:

  A DataFrame with columns `id_curve`, `tobs`, `X`.

- id_lag:

  Integer id of the conditioning curve.

- bw_grid:

  Bandwidth grid for the adaptive risk.

- rho:

  Design weights of the conditioning curve.

- homoscedastic:

  Whether to use a constant noise variance.

- tikhonov:

  Tikhonov regularisation parameter.

- n_subgrid_bw:

  Number of points per axis of the bandwidth sub-grid.

- kernel_name:

  Kernel name.

## Value

A list with the cached bandwidths, the covariance operator, the mean,
the noise level, the regularised variance matrix and the residual.
