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
  id_conditioning_curve,
  bw_grid,
  rho,
  homoscedastic,
  tikhonov,
  bw_subgrid_size,
  kernel_name,
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL
)
```

## Arguments

- data:

  A DataFrame with columns `id_curve`, `tobs`, `X`.

- id_conditioning_curve:

  Integer id of the conditioning curve.

- bw_grid:

  Bandwidth grid for the adaptive risk.

- rho:

  Design weights of the conditioning curve.

- homoscedastic:

  Whether to use a constant noise variance.

- tikhonov:

  Tikhonov regularisation parameter.

- bw_subgrid_size:

  Number of points per axis of the bandwidth sub-grid.

- kernel_name:

  Kernel name.

- presmooth_bw:

  Numeric (positive vector or scalar). Bandwidth used to presmooth each
  curve in the local regularity step, see `estimate_locreg_cpp`. Default
  `NULL` selects it by cross-validation.

- Delta:

  Numeric (positive). Length of the neighborhood of each point used in
  the local regularity step. Default `NULL` estimates it from the data.

- presmooth_bw_grid:

  Numeric vector. Candidate bandwidths of the cross-validation that
  selects `presmooth_bw`. Default `NULL` uses the default grid.

- presmooth_nsubset:

  Integer (positive). Number of curves used by that cross-validation.
  Default `NULL` uses min(70, floor(N / 2)) curves.

## Value

A list with the cached bandwidths, the covariance operator, the mean,
the noise level, the regularised variance matrix and the residual.
