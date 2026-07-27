# Estimate the Mean Function by the Rubìn-Panaretos Method

Estimates the mean function with the local-linear estimator of Rubìn and
Panaretos (2020), which pools the observation points of all curves and
smooths them with a single bandwidth. It is provided for comparison with
the adaptive estimator of
[estimate_mean](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md).

## Usage

``` r
estimate_mean_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw,
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

- t:

  `vector (numeric)`. Points of \\\[0, 1\]\\ at which the mean function
  is estimated.

- bw:

  `numeric (positive scalar)`. Bandwidth of the estimator, common to
  every point of `t`. See
  [estimate_mean_bw_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md)
  to select it by cross-validation.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per point of `t` and columns:

- `t`: the point at which the mean function is estimated.

- `bw`: the bandwidth used.

- `muhat_RP`: the estimated mean function.

## References

Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
series: estimation and prediction. *Electronic Journal of Statistics*,
14(1), 1137–1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690)

## See also

[`estimate_mean_bw_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_bw_rp.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md).

## Examples

``` r
data("data_far")

dt_mean_rp <- estimate_mean_rp(
  data = data_far[data_far$id_curve <= 20, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw = 5/70, kernel_name = "epanechnikov")
dt_mean_rp
#>        t         bw muhat_RP
#>    <num>      <num>    <num>
#> 1:  0.25 0.07142857 243.3595
#> 2:  0.50 0.07142857 241.1444
#> 3:  0.75 0.07142857 241.4589
```
