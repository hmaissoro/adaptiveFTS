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
  n_subgrid_bw = 10L,
  n_cv_tikhonov = 30L,
  id_lag = NULL,
  tikhonov = NULL,
  tikhonov_grid = NULL,
  bw_grid = NULL,
  density_bw = NULL
)
```

## Arguments

- data:

  A `data.table` (or `data.frame`), a `list` of `data.table` (or
  `data.frame`), or a `list` of `list`.

  - If `data.table`: It should contain the raw curve observations in at
    least three columns.

    - `idcol` : The name of the column containing the curve index in the
      sample. Each curve index is repeated according to the number of
      observation points.

    - `tcol` : The name of the column with observation points associated
      with each curve index.

    - `ycol` : The name of the column with observed values at each
      observation point for each curve index.

  - If `list` of `data.table`: In this case, each element in the `list`
    represents the observation data of a curve in the form of a
    `data.table` or `data.frame`. Each `data.table` contains at least
    two columns.

    - `tcol` : The name of the column with observation points for the
      curve.

    - `ycol` : The name of the column with observed values for the
      curve.

  - If `list` of `list`: In this case, `data` is a list where each
    element is the observation data of a curve, given as a `list` of two
    vectors.

    - `tcol` : The vector containing observation points for the curve.

    - `ycol` : The vector containing observed values for the curve.

- idcol:

  `character`. If `data` is given as a `data.table` or `data.frame`,
  this is the name of the column that holds the curve index. Each curve
  index is repeated according to the number of observation points. If
  `data` is a `list` of `data.table` (or `data.frame`) or a `list` of
  `list`, set `idcol = NULL`.

- tcol:

  `character`. The name of the column (or vector) containing the
  observation points for the curves.

- ycol:

  `character`. The name of the column with observed values for the
  curves.

- t:

  Numeric vector of prediction points in \\\[0, 1\]\\.

- horizon:

  Integer prediction horizon (steps ahead). Default `1`.

- kernel_name:

  Kernel name. Default `"epanechnikov"`.

- homoscedastic:

  If `TRUE` (default) a constant noise variance (median of the pointwise
  estimates) is used; otherwise the t-varying estimates.

- n_subgrid_bw:

  Number of points per axis of the coarse sub-grid on which the adaptive
  bandwidths are selected. Default `10`.

- n_cv_tikhonov:

  Number of trailing curves used for the Tikhonov cross-validation when
  `tikhonov` is `NULL`. Default `30`.

- id_lag:

  Integer id of the conditioning curve. Its successor is the curve to be
  predicted. Default `NULL` uses the last curve in `data`.

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
