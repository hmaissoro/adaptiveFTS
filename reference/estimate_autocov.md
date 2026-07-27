# Estimate the Covariance or Autocovariance Function

Estimates the adaptive lag-\\\ell\\ autocovariance function for \\\ell =
0, 1, \ldots\\ (\\\ell = 0\\ being the covariance function), using at
each pair of points the bandwidths that minimise the estimated risk.

## Usage

``` r
estimate_autocov(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  bw_s = NULL,
  bw_t = NULL,
  bw_grid = NULL,
  common_bw = FALSE,
  center_curves = TRUE,
  correct_diagonal = TRUE,
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

  `integer (non-negative)`. Lag \\\ell\\ of the autocovariance;
  `lag = 0` gives the covariance function.

- bw_s:

  `vector (numeric)`. Bandwidth to use for `s` at each pair. Default
  `NULL` selects it by minimising the estimated risk.

- bw_t:

  `vector (numeric)`. Bandwidth to use for `t` at each pair. Default
  `NULL` selects it by minimising the estimated risk.

- bw_grid:

  `vector (numeric)`. Candidate bandwidths, from which estimate_autocov
  picks the risk-minimising one for each pair. Default `NULL` builds the
  grid from the data; see Details.

- common_bw:

  `logical`. If `TRUE`, a single bandwidth is selected for both
  arguments of the autocovariance; if `FALSE` (default), one bandwidth
  per argument. See Details.

- center_curves:

  `logical`. If `TRUE` (default), the curves are centred before
  smoothing.

- correct_diagonal:

  `logical`. If `TRUE` (default), the observation-noise variance is
  subtracted from the diagonal when `lag = 0`. See Details.

- kernel_name:

  `string`. Kernel of the smoothing estimator, one of "epanechnikov"
  (default), "biweight", "triweight", "tricube", "triangular" and
  "uniform".

## Value

A `data.table` with one row per pair (`s`, `t`) and columns:

- `s`, `t`: the arguments of the autocovariance function.

- `optbw_s`, `optbw_t`: the bandwidths used for `s` and for `t`.
  Identical when `common_bw = TRUE`.

- `Hs`, `Ls2`: the estimated local exponent \\H_s\\ and squared Hölder
  constant \\L_s^2\\ at `s`.

- `Ht`, `Lt2`: the same quantities at `t`.

- `PNs`, `muhat_s`: the number of curves used for the mean at `s` and
  the estimated mean there.

- `PNt`, `muhat_t`: the same quantities at `t`.

- `PNl`: the number of curves contributing to the estimate at (`s`,
  `t`), \\P\_{N,\ell}(s,t;h_s,h_t)\\.

- `autocov`: the estimated (auto)covariance.

## Details

Unless `bw_s` and `bw_t` are supplied, the bandwidths are selected pair
by pair by minimising the risk of
[estimate_autocov_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
over `bw_grid`. Set `common_bw = TRUE` for the one-bandwidth estimator
of Maissoro, Patilea and Vimond (2025) and `FALSE` (default) for the
two-bandwidth estimator of Maissoro, Patilea and Vimond (2026).

`center_curves` chooses how the mean is removed. With `TRUE` (default)
the curves are centred before smoothing, which estimates
\\\mathbb{E}(X_0(s) - \mu(s))(X\_{\ell}(t) - \mu(t))\\ in one pass. With
`FALSE` the two pieces \\\mathbb{E}X_0(s)X\_{\ell}(t)\\ and
\\\mu(s)\mu(t)\\ are estimated separately, the first with a bandwidth
from
[estimate_autocov_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
and the second with a bandwidth from
[estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md).
Both are centred estimates; they differ in which bandwidth is applied to
which piece.

At `lag = 0` the observation noise contributes to the estimate wherever
the two smoothing windows overlap, that is on and near the diagonal \\s
= t\\. `correct_diagonal = TRUE` (default) subtracts that contribution:
an estimate of \\\sigma(s)\sigma(t)\\ weighted by the overlap of the two
kernel weight vectors, which is largest at \\s = t\\ and decays as the
points move apart. It has no effect for `lag > 0`, where the noise of
two distinct curves is uncorrelated.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md),
[`estimate_facf()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_facf.md),
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md).

## Examples

``` r
data("data_far")
dt_small <- data_far[data_far$id_curve <= 20, ]
bwg <- seq(0.04, 0.15, length.out = 5)

# Lag-1 autocovariance.
dt_autocov <- estimate_autocov(
  data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5), t = c(1/4, 1/2), lag = 1, bw_grid = bwg,
  common_bw = FALSE, center_curves = TRUE, correct_diagonal = FALSE,
  kernel_name = "epanechnikov")
dt_autocov[, list(s, t, optbw_s, optbw_t, PNl, autocov)]
#>        s     t optbw_s optbw_t   PNl    autocov
#>    <num> <num>   <num>   <num> <num>      <num>
#> 1:   0.2  0.25    0.04    0.04    19  0.9357060
#> 2:   0.4  0.50    0.04    0.04    19 -0.1754777

# Covariance, with the noise variance removed from the diagonal.
dt_cov <- estimate_autocov(
  data = dt_small, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/4, 1/2), t = c(1/4, 1/2), lag = 0, bw_grid = bwg,
  common_bw = FALSE, center_curves = TRUE, correct_diagonal = TRUE,
  kernel_name = "epanechnikov")
dt_cov[, list(s, t, PNl, autocov)]
#>        s     t   PNl  autocov
#>    <num> <num> <num>    <num>
#> 1:  0.25  0.25    20 3.826240
#> 2:  0.50  0.50    20 5.198347
```
