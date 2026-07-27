# Select the Bandwidth of the Rubìn-Panaretos Autocovariance Estimator

Selects the bandwidth of
[estimate_autocov_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md)
by \\K\\-fold cross-validation over the curves, as described in Rubìn
and Panaretos (2020). Each fold is scored by the squared error between
the empirical cross-products of the held-out curves and the lag-0
autocovariance estimated on the others.

## Usage

``` r
estimate_autocov_bw_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  n_folds = 10,
  bw_grid = seq(0.001, 0.15, len = 45),
  bw_mean = NULL,
  mean_rp = NULL,
  kernel_name = "epanechnikov"
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

- n_folds:

  `integer (positive)`. Number of cross-validation folds.

- bw_grid:

  `vector (numeric)`. Candidate bandwidths.

- bw_mean:

  `numeric (positive scalar)`. Bandwidth of the mean function estimator,
  used only when `mean_rp` is `NULL`.

- mean_rp:

  `data.table`. Mean function estimated at every observation point of
  every curve, with columns `id_curve`, `tobs` and `muhat_RP`. Default
  `NULL` estimates it from `bw_mean`.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per candidate bandwidth and columns:

- `bw`: the candidate bandwidth.

- `cv_error`: the cross-validation error at `bw`. The bandwidth
  minimising it is the one to pass to
  [estimate_autocov_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md).

## Details

Every candidate bandwidth requires a Rubìn-Panaretos estimate at every
pair of observation points of the held-out curves, so the runtime grows
with the fourth power of the number of points per curve. Keep `bw_grid`
short and the number of curves small.

## References

Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
series: estimation and prediction. *Electronic Journal of Statistics*,
14(1), 1137–1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690)

## See also

[`estimate_autocov_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_rp.md),
[`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md).
