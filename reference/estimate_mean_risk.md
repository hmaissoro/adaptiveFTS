# Estimate the Risk of the Mean Function Estimator

Estimates the risk \\R\_\mu(t;h)\\ of the adaptive mean function
estimator over a grid of candidate bandwidths, as described in Section
4.1 of Maissoro, Patilea and Vimond (2025). Minimising it over `h` at
each `t` is what
[estimate_mean](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md)
does to select its bandwidth.

## Usage

``` r
estimate_mean_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw_grid = NULL,
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL,
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

- bw_grid:

  `vector (numeric)`. Candidate bandwidths, from which
  [estimate_mean](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md)
  picks the risk-minimising one at each `t`. Default `NULL` builds the
  grid from the data; see Details.

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

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per (`t`, `h`) pair and columns:

- `t`: the point at which the risk is estimated.

- `h`: the candidate bandwidth.

- `PN`: the number of curves contributing to the estimate at `t`,
  \\P_N(t;h)\\.

- `locreg_bw`: the bandwidth used to estimate the local regularity.

- `Ht`: the estimated local exponent \\H_t\\.

- `Lt2`: the estimated squared Hölder constant \\L_t^2\\.

- `bias_term`, `variance_term`, `dependence_term`: the three components
  of the risk.

- `mean_risk`: the estimated risk.

## Details

The risk bound splits into three terms, returned separately so that the
selected bandwidth can be traced back to what drove it: a bias term
growing with \\h^{2H_t}\\ through the local regularity, a variance term
decreasing in \\h\\ through the number of usable points, and a
dependence term reflecting the serial dependence between curves. The
local regularity parameters are estimated internally at each `t`, so
`Ht` and `Lt2` are reported alongside the risk.

Left to `NULL`, `bw_grid` is a 20-point geometric grid running from
\\4(N\widehat\lambda)^{-0.9}\\ to \\4(N\widehat\lambda)^{-1/3}\\, where
\\N\\ is the number of curves and \\\widehat\lambda\\ the average number
of observation points per curve.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

## See also

[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md).

## Examples

``` r
data("data_far")

dt_mean_risk <- estimate_mean_risk(
  data = data_far[data_far$id_curve <= 20, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.02, 0.15, length.out = 8),
  kernel_name = "epanechnikov")

# The risk-minimising bandwidth at each t.
dt_mean_risk[, list(h = h[which.min(mean_risk)]), by = "t"]
#>        t     h
#>    <num> <num>
#> 1:  0.25  0.02
#> 2:  0.50  0.02
#> 3:  0.75  0.02
```
