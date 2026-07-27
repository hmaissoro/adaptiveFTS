# Fit and predict the adaptive functional BLUP in one call

Convenience wrapper that fits the adaptive BLUP on `data` and
immediately predicts the curve following the conditioning curve at `t`.
Equivalent to `predict(blup_fit(data, ...), t = t, horizon = horizon)`,
plus the Tikhonov selection carried alongside the prediction.

## Usage

``` r
blup(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = seq(0.01, 0.99, length.out = 99),
  horizon = 1L,
  kernel_name = "epanechnikov",
  homoscedastic = TRUE,
  bw_subgrid_size = 10L,
  n_cv_curves = 30L,
  id_conditioning_curve = NULL,
  tikhonov = NULL,
  tikhonov_grid = NULL,
  bw_grid = NULL,
  density_bw = NULL
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

- t:

  Numeric vector of prediction points in \\\[0, 1\]\\.

- horizon:

  Integer prediction horizon (steps ahead). Default `1`.

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

## Value

An object of class `blup`: a list with

- `prediction`: a `data.table` with columns `horizon`, `t`, `muhat` and
  `prediction` (one block of rows per horizon; see
  [`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md)).

- `tikhonov`: the Tikhonov parameter used.

- `tikhonov_cv`: the
  [`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md)
  output when `tikhonov` was selected, otherwise `NULL`.

## See also

[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md),
[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md),
[`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md).
