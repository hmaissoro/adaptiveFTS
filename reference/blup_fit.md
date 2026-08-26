# Fit the adaptive functional BLUP

Estimates every component of the adaptive Best Linear Unbiased Predictor
that does not depend on the prediction points.

## Usage

``` r
blup_fit(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  kernel_name = "epanechnikov",
  homoscedastic = TRUE,
  bw_subgrid_size = 10L,
  n_cv_curves = 30L,
  id_conditioning_curve = NULL,
  tikhonov = NULL,
  tikhonov_grid = NULL,
  bw_grid = NULL,
  density_bw = NULL,
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL
)
```

## Arguments

- data:

  Raw curve observations, as a `data.table` (or `data.frame`) in long
  format, or as a `list` with one element per curve. See
  [`format_data`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  for the accepted layouts and for the `id_curve` / `tobs` / `X` columns
  they are converted to.

- idcol:

  `character(1)` or `NULL`. Name of the column holding the curve index
  when `data` is a single table. Must be `NULL` when `data` is a list of
  curves.

- tcol:

  `character(1)`. Name of the column (or vector) holding the observation
  points of the curves.

- ycol:

  `character(1)`. Name of the column (or vector) holding the values
  observed at those points.

- kernel_name:

  Kernel name. Default `"epanechnikov"`.

- homoscedastic:

  If `TRUE` (default) a constant noise variance (median of the pointwise
  estimates) is used; otherwise the t-varying estimates.

- bw_subgrid_size:

  Number of points per axis of the coarse sub-grid on which the adaptive
  bandwidths are selected. Default `10`.

- n_cv_curves:

  Number of trailing curves held out for the Tikhonov cross-validation
  when `tikhonov` is `NULL`. Default `30`.

- id_conditioning_curve:

  Integer id of the curve conditioned on. Its successor is the curve
  that [`predict()`](https://rdrr.io/r/stats/predict.html) reconstructs.
  Default `NULL` uses the last curve in `data`.

- tikhonov:

  Tikhonov regularisation parameter \\\alpha\\. Default `NULL` selects
  it by cross-validation (see
  [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md))
  over `tikhonov_grid`; pass a numeric value to use it directly.

- tikhonov_grid:

  Candidate values for the cross-validation when `tikhonov` is `NULL`.
  Default `NULL` uses the default grid of
  [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md).

- bw_grid:

  Bandwidth grid for the adaptive mean/(auto)covariance risk. Default
  `NULL` sets a geometric grid from the data.

- density_bw:

  Optional fixed design-density bandwidth reused for every
  `estimate_density` call (independent design only). Default `NULL`
  selects it once via
  [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md).

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth of the
  Nadaraya-Watson estimator used to presmooth each curve before the
  regularity is estimated. A scalar applies the same bandwidth to every
  curve; a vector must hold one bandwidth per curve, in the order the
  curves appear in `data`. Default `NULL` selects a single bandwidth by
  cross-validation over every curve, as
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md)
  does.

- Delta:

  `numeric (positive)`. Length of the neighbourhood around each point of
  `t` used to estimate the local regularity. Default `NULL` sets it from
  the data; see Details.

- presmooth_bw_grid:

  `vector (numeric)`. Candidate bandwidths of the cross-validation that
  selects `presmooth_bw` when the latter is `NULL`. Default `NULL` uses
  the default grid of
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).
  Ignored when `presmooth_bw` is supplied.

- presmooth_nsubset:

  `integer (positive)`. Number of curves used by that cross-validation.
  Default `NULL` uses `min(70, floor(N / 2))` curves, where \\N\\ is the
  number of curves. Lower it to speed up the selection on large samples.
  Ignored when `presmooth_bw` is supplied.

## Value

An object of class `blup_fit`: a list whose main elements are:

- `Tn0`, `Yn0`: the conditioning-curve design points and values.

- `rho`, `root_Dn0`: the design weights and their square-root matrix.

- `muhat_Tn0`, `c0hat`, `sigma2`: the mean, covariance operator and
  noise level of the conditioning curve.

- `V`, `resid`: the regularised variance matrix and the conditioning
  residual.

- `opt_mean`, `opt_cov`, `opt_autocov`: the cached adaptive bandwidths.

- `data`, `kernel_name`, `is_common_design`, `density_bw`, `bw_grid`,
  `tikhonov`, `homoscedastic`: the information needed by
  [`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md).

- `tikhonov_cv`: the
  [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md)
  output when `tikhonov` was selected, otherwise `NULL`.

## Details

The fit conditions on a single curve (its immediate successor is what
[`predict()`](https://rdrr.io/r/stats/predict.html) reconstructs) and
caches the adaptive bandwidths, so that
[`predict()`](https://rdrr.io/r/stats/predict.html) only re-runs the
cheap plug-in estimates at the requested prediction points. The
covariance assembly is a single lag-1 block.

## See also

[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md),
[`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md).
