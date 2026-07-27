# Estimate the Autocovariance Function by the Rubìn-Panaretos Method

Estimates the lag-\\\ell\\ autocovariance function with the local-linear
estimator of Rubìn and Panaretos (2020), which smooths every pair of
observation points of curves \\\ell\\ apart with a single bandwidth. It
is provided for comparison with the adaptive estimator of
[estimate_autocov](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).

## Usage

``` r
estimate_autocov_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  bw,
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

- s:

  `vector (numeric)`. First argument of the autocovariance function: the
  points `s` of the pairs (`s`, `t`). Must have the same length as `t`.

- t:

  `vector (numeric)`. Second argument of the autocovariance function:
  the points `t` of the pairs (`s`, `t`). Must have the same length as
  `s`.

- lag:

  `integer (non-negative)`. Lag \\\ell\\ of the autocovariance.

- bw:

  `numeric (positive scalar)`. Bandwidth of the estimator, common to
  every pair. See
  [estimate_autocov_bw_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md)
  to select it by cross-validation.

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

A `data.table` with one row per pair (`s`, `t`) and columns:

- `s`, `t`: the arguments of the autocovariance function.

- `lag`: the lag \\\ell\\.

- `bw_mean`: the bandwidth used for the mean function.

- `bw`: the bandwidth used for the autocovariance.

- `autocovhat_rp`: the estimated autocovariance.

## References

Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
series: estimation and prediction. *Electronic Journal of Statistics*,
14(1), 1137–1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690)

## See also

[`estimate_autocov_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_bw_rp.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).

## Examples

``` r
# \donttest{
data("data_far")

dt_autocov_rp <- estimate_autocov_rp(
  data = data_far[data_far$id_curve <= 10, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1,
  bw = 0.1, bw_mean = 0.1, mean_rp = NULL, kernel_name = "epanechnikov")
dt_autocov_rp
#>        s     t   lag bw_mean autocovhat_rp
#>    <num> <num> <num>   <num>         <num>
#> 1:   0.2  0.25     1     0.1    0.29505726
#> 2:   0.4  0.50     1     0.1   -0.08008549
# }
```
