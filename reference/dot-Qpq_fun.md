# Weighted Cross-Product \\Q\_{pq}^{(\ell)}\\ of the Rubìn-Panaretos Estimator

Computes the \\Q\_{pq}^{(\ell)}\\ term of Equation (B.7) of Rubìn and
Panaretos (2020), the counterpart of
[.Spq_fun](https://hmaissoro.github.io/adaptiveFTS/reference/dot-Spq_fun.md)
weighting the centred cross-products of the observed values.

## Usage

``` r
.Qpq_fun(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = 1/4,
  t = 1/2,
  lag = 1,
  p = 1,
  q = 1,
  bw,
  mean_rp = NULL,
  bw_mean = NULL,
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

- s:

  `numeric (scalar)`. First argument of the autocovariance function.

- t:

  `numeric (scalar)`. Second argument of the autocovariance function.

- lag:

  `integer (non-negative)`. Lag of the autocovariance.

- p, q:

  `numeric (integer)`. Exponents of the centred and scaled observation
  points in the sum.

- bw:

  `numeric (positive scalar)`. Bandwidth of the estimator.

- mean_rp:

  `data.table`. Mean function estimated at every observation point of
  every curve, with columns `id_curve`, `tobs` and `muhat_RP`. Default
  `NULL` estimates it from `bw_mean`.

- bw_mean:

  `numeric (positive scalar)`. Bandwidth of the mean function estimator,
  used only when `mean_rp` is `NULL`.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `numeric` scalar.

## References

Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
series: estimation and prediction. *Electronic Journal of Statistics*,
14(1), 1137–1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690)
