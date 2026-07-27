# Estimate Empirical Autocovariance Function

This function estimates the empirical autocovariance function used in
the empirical study section of the papers Maissoro, Patilea and Vimond
(2025) and Maissoro, Patilea and Vimond (2026).

## Usage

``` r
estimate_empirical_autocov(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  lag = c(0, 1, 2),
  presmooth_bw = NULL,
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

  `vector (numeric)`. Observation points at which we want to estimate
  the empirical autocovariance function.

- lag:

  `vector (integer)`. Lag of the autocovariance.

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth used to presmooth
  each curve before the estimation. A scalar applies the same bandwidth
  to every curve; a vector must hold one bandwidth per curve, in the
  order the curves appear in `data`. Default `NULL` selects it by
  cross-validation, see
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` with three columns: `t`, `lag`, and `autocov`
corresponding to the estimated autocovariance.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Examples

``` r
# Load data
data("data_far")

# Estimate empirical autocovariance with a specified bandwidth
dt_empirical_autocov <- estimate_empirical_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), lag = c(1, 2), presmooth_bw = 0.1,
  kernel_name = "epanechnikov")
dt_empirical_autocov
#>        t   lag  autocov
#>    <num> <num>    <num>
#> 1:  0.25     1 5.090136
#> 2:  0.50     1 4.937906
#> 3:  0.75     1 5.140885
#> 4:  0.25     2 2.836493
#> 5:  0.50     2 2.520016
#> 6:  0.75     2 2.168116

# Estimate empirical autocovariance with Cross-Validation bandwidth selection
dt_empirical_autocov_cv <- estimate_empirical_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), lag = c(1, 2), presmooth_bw = NULL,
  kernel_name = "epanechnikov")
dt_empirical_autocov_cv
#>        t   lag  autocov
#>    <num> <num>    <num>
#> 1:  0.25     1 5.379540
#> 2:  0.50     1 5.112767
#> 3:  0.75     1 5.357566
#> 4:  0.25     2 2.858931
#> 5:  0.50     2 2.700018
#> 6:  0.75     2 2.446969

```
