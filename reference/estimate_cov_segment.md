# Estimate Covariance Segment Function for Functional Data

Estimates the covariance segment function \\\Gamma\_{N,0}(t,t;h_t,h_t)\\
for functional data using the Nadaraya–Watson estimator with a specified
kernel. This is part of the methodology described in Maissoro, Patilea
and Vimond (2026).

## Usage

``` r
estimate_cov_segment(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw = NULL,
  bw_grid = NULL,
  center_curves = TRUE,
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

  A numeric vector. Observation points where the mean function of the
  underlying process is estimated.

- bw:

  A numeric vector. Bandwidth to use at each point of `t`, recycled if a
  scalar. Default `NULL` selects it by minimising the risk estimated by
  [estimate_cov_segment_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md).

- bw_grid:

  A numeric vector. A bandwidth grid from which the best smoothing
  parameter is selected for each `t`. Default is `NULL`, in which case
  it is defined as an exponential grid of \\N \times \lambda\\.

- center_curves:

  Logical. If `TRUE` (default), the curves are centred before smoothing.
  It governs the moment and autocovariance estimators only: the local
  regularity step always centres the curves.

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

  Character string. Specifies the kernel function for estimation;
  default is `"epanechnikov"`. Supported kernels include:
  `"epanechnikov"`, `"biweight"`, `"triweight"`, `"tricube"`,
  `"triangular"`, and `"uniform"`.

## Value

A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
containing the following columns:

- `t` : The observation points at which the covariance segment function
  is estimated.

- `optbw` : The optimal bandwidth used to estimate covariance segment
  function at each `t`.

- `Ht` : Local exponent estimates for each `t`, corresponding to
  \\H_t\\.

- `Lt2` : Estimates of the Hölder constant for each `t`, corresponding
  to \\L_t^2\\.

- `PN` : The number of selected curves used in the estimation for each
  `t`.

- `cov_segment_hat` : Uncorrected covariance segment estimate.

- `covseg_correction` : Correction term based on measurement error
  variance.

- `cov_segment_hat_corrected` : Final corrected covariance segment
  estimate.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[estimate_cov_segment_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment_risk.md),
[estimate_locreg](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[estimate_sigma](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
[estimate_nw](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[estimate_empirical_autocov](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_autocov.md)

## Examples

``` r
data("data_far")

dt_cov_segment <- estimate_cov_segment(
  data = data_far[data_far$id_curve <= 20, ],
  idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.04, 0.15, length.out = 5),
  center_curves = TRUE, kernel_name = "epanechnikov")
dt_cov_segment
#>        t optbw        Ht       Lt2    PN cov_segment_hat covseg_correction
#>    <num> <num>     <num>     <num> <num>           <num>             <num>
#> 1:  0.25  0.04 0.5116001  7.372248    20        3.855708        0.02946757
#> 2:  0.50  0.04 1.0000000 13.320572    20        5.244993        0.04664560
#> 3:  0.75  0.04 0.7944922 10.397411    20        6.313028        0.12520059
#>    cov_segment_hat_corrected
#>                        <num>
#> 1:                  3.826240
#> 2:                  5.198347
#> 3:                  6.187827
```
