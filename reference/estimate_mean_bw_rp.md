# Select the Bandwidth of the Rubìn-Panaretos Mean Estimator

Selects the bandwidth of
[estimate_mean_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md)
by \\K\\-fold cross-validation over the curves, as described in Rubìn
and Panaretos (2020). Curves are split into folds; each fold is
predicted from the mean function estimated on the others, and the
squared prediction errors are averaged.

## Usage

``` r
estimate_mean_bw_rp(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  n_folds = 10,
  bw_grid = seq(0.001, 0.15, len = 45),
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

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per candidate bandwidth and columns:

- `bw`: the candidate bandwidth.

- `cv_error`: the cross-validation error at `bw`. The bandwidth
  minimising it is the one to pass to
  [estimate_mean_rp](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md).

## References

Rubìn, T. and Panaretos, V. M. (2020). Sparsely observed functional time
series: estimation and prediction. *Electronic Journal of Statistics*,
14(1), 1137–1210.
[doi:10.1214/20-EJS1690](https://doi.org/10.1214/20-EJS1690)

## See also

[`estimate_mean_rp()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_rp.md).

## Examples

``` r
# \donttest{
data("data_far")
dt_small <- data_far[data_far$id_curve <= 10, ]

dt_bw <- estimate_mean_bw_rp(
  data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
  n_folds = 5, bw_grid = seq(0.02, 0.15, length.out = 5),
  kernel_name = "epanechnikov")

dt_mean_rp <- estimate_mean_rp(
  data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw = dt_bw[, bw[which.min(cv_error)]],
  kernel_name = "epanechnikov")
dt_mean_rp
#>        t    bw muhat_RP
#>    <num> <num>    <num>
#> 1:  0.25 0.085 242.8988
#> 2:  0.50 0.085 240.9457
#> 3:  0.75 0.085 240.7357
# }
```
