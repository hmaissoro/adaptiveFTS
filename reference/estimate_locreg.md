# Estimate the Local Regularity Parameters

Estimates the local regularity parameters \\H_t\\ and \\L_t^2\\ of the
underlying process at each point of `t`, following Section 3 of
Maissoro, Patilea and Vimond (2025). \\H_t\\ is the local Hölder
exponent and \\L_t^2\\ the squared Hölder constant; both drive the
bandwidths of every adaptive estimator of the package.

## Usage

``` r
estimate_locreg(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = 1/2,
  Delta = NULL,
  presmooth_bw = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL,
  kernel_name = "epanechnikov",
  center = TRUE
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

  `vector (numeric)`. Points of \\\[0, 1\]\\ at which the local
  regularity parameters are estimated.

- Delta:

  `numeric (positive)`. Length of the neighbourhood around each point of
  `t` used to estimate the local regularity. Default `NULL` sets it from
  the data; see Details.

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth of the
  Nadaraya-Watson estimator used to presmooth each curve before the
  regularity is estimated. A scalar applies the same bandwidth to every
  curve; a vector must hold one bandwidth per curve, in the order the
  curves appear in `data`. Default `NULL` selects a single bandwidth by
  cross-validation over every curve, as
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md)
  does.

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

  `string`. Kernel of the presmoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

- center:

  `logical`. If `TRUE` (default), the curves are centred before the
  regularity is estimated.

## Value

A `data.table` with one row per point of `t` and columns:

- `t`: the point at which the local regularity is estimated.

- `locreg_bw`: the presmoothing bandwidth used.

- `Delta`: the length of the neighbourhood used around `t`.

- `Nused`: the number of curves that contributed a non-degenerate
  estimate at `t`.

- `Ht`: the estimated local exponent \\H_t\\.

- `Lt2`: the estimated squared Hölder constant \\L_t^2\\.

## Details

Each curve is presmoothed and evaluated at three points spread over a
neighbourhood of length `Delta` around `t`, and the regularity is read
off the ratio of the mean squared increments between those points.
`Delta` drives the bias-variance trade-off of that comparison: too small
and the three points carry the same information, too large and the local
regularity is averaged away. Left to `NULL` it is set from the average
number of observation points per curve, \\\widehat\lambda\\, as
\\\min\\\exp(-(\log\widehat\lambda)^{1/3}),\\ 0.2\\\\. Near the
boundaries of \\\[0, 1\]\\ the neighbourhood is shifted inwards rather
than truncated.

Curves whose presmoothed values fall outside the 2.5\\ any of the three
points are discarded, so `Nused` is smaller than the number of curves.
The exponent is clamped to \\\[0.1, 1\]\\: a boundary value usually
means the neighbourhood or the presmoothing bandwidth is unsuited to the
data rather than a genuinely extreme regularity.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md).

## Examples

``` r
data("data_far")

dt_locreg <- estimate_locreg(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = seq(0.2, 0.8, length.out = 5), Delta = NULL, presmooth_bw = NULL,
  kernel_name = "epanechnikov", center = TRUE)
dt_locreg
#>        t  locreg_bw     Delta Nused        Ht      Lt2
#>    <num>      <num>     <num> <num>     <num>    <num>
#> 1:  0.20 0.01740629 0.1923983   118 0.5645557 9.870244
#> 2:  0.35 0.01740629 0.1923983   127 0.5810445 7.120796
#> 3:  0.50 0.01740629 0.1923983   122 0.6187842 8.051270
#> 4:  0.65 0.01740629 0.1923983   120 0.5229759 4.141653
#> 5:  0.80 0.01740629 0.1923983   127 0.4184367 2.174434
```
