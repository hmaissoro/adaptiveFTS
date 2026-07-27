# Estimate the Mean Function

Estimates the mean function of the underlying process with the adaptive
estimator of Maissoro, Patilea and Vimond (2025), using at each point
the bandwidth that minimises the estimated risk.

## Usage

``` r
estimate_mean(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw = NULL,
  bw_grid = NULL,
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

  `vector (numeric)`. Points of \\\[0, 1\]\\ at which the risk is
  estimated.

- bw:

  `vector (numeric)`. Bandwidth to use at each point of `t`, recycled if
  a scalar. Default `NULL` selects it by minimising the risk estimated
  by
  [estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md).

- bw_grid:

  `vector (numeric)`. Candidate bandwidths, from which estimate_mean
  picks the risk-minimising one at each `t`. Default `NULL` builds the
  grid from the data; see Details.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per point of `t` and columns:

- `t`: the point at which the mean function is estimated.

- `optbw`: the bandwidth used at `t`.

- `Ht`: the estimated local exponent \\H_t\\.

- `Lt2`: the estimated squared Hölder constant \\L_t^2\\.

- `PN`: the number of curves contributing to the estimate at `t`.

- `muhat`: the estimated mean function.

## Details

Unless `bw` is supplied, the bandwidth is selected point by point by
minimising the risk of
[estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md)
over `bw_grid`, so neighbouring points may be smoothed differently
according to the local regularity of the process. Supplying `bw` skips
the risk estimation entirely, which is worth doing when the same
bandwidths are reused across many calls.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

## See also

[`estimate_mean_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).

## Examples

``` r
data("data_far")

dt_mean <- estimate_mean(
  data = data_far[data_far$id_curve <= 20, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.02, 0.15, length.out = 8),
  kernel_name = "epanechnikov")
dt_mean
#>        t optbw        Ht       Lt2    PN    muhat
#>    <num> <num>     <num>     <num> <num>    <num>
#> 1:  0.25  0.02 0.5116001  7.372248    20 243.5573
#> 2:  0.50  0.02 1.0000000 13.320572    19 241.4676
#> 3:  0.75  0.02 0.7944922 10.397411    18 241.3353

summary(dt_mean)
#> Adaptive mean function estimate
#>   Evaluation points  : 3 (t in [0.25, 0.75])
#>   Training curves    : 20
#>   Kernel             : epanechnikov
#>   Optimal bandwidth  : [0.02, 0.02]
#>   Curves used (PN)   : [18, 20]
#>   muhat              : [241, 244]
```
