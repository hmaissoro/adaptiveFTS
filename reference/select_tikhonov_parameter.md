# Select the Tikhonov regularisation parameter

Selects the Tikhonov regularisation parameter \\\alpha\\ of the adaptive
BLUP. Currently only `method = "cv"` is implemented, a one-step-ahead
cross-validation.

## Usage

``` r
select_tikhonov_parameter(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  method = c("cv"),
  kernel_name = "epanechnikov",
  homoscedastic = TRUE,
  bw_subgrid_size = 10L,
  n_cv_curves = 30L,
  tikhonov_grid = NULL,
  bw_grid = NULL,
  density_bw = NULL,
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL
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

- method:

  Selection method. Currently only `"cv"` (one-step-ahead
  cross-validation) is available.

- kernel_name:

  Kernel name. Default `"epanechnikov"`.

- homoscedastic:

  If `TRUE` (default) a constant noise variance (median of the pointwise
  estimates) is used; otherwise the t-varying estimates.

- bw_subgrid_size:

  Number of points per axis of the coarse sub-grid on which the adaptive
  bandwidths are selected. Default `10`.

- n_cv_curves:

  Number of trailing curves held out for the Tikhonov cross-validation
  when `tikhonov` is `NULL`. Default `30`.

- tikhonov_grid:

  Candidate values. If `NULL`, a 25-point grid \\\\e^{-5}, \ldots,
  e^3\\\\ is used.

- bw_grid:

  Bandwidth grid for the adaptive mean/(auto)covariance risk. Default
  `NULL` sets a geometric grid from the data.

- density_bw:

  Optional fixed design-density bandwidth reused for every
  `estimate_density` call (independent design only). Default `NULL`
  selects it once via
  [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md).

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

## Value

A list with:

- `tikhonov_star`: the selected Tikhonov parameter.

- `tikhonov_grid`: the candidate grid.

- `cv_curve`: the mean cross-validation score per candidate.

- `cv_matrix`: the per-fold cross-validation scores (fold by candidate).

- `val_ids`: the ids of the validation curves.

## Details

Each of the last `n_cv_curves` curves is predicted from its immediate
predecessor and scored by the design-weighted squared prediction error
at its observation points, \\\sum_i \varrho\_{n,i}\\(Y\_{n,i} - \widehat
X_n(T\_{n,i};\alpha))^2\\, where \\\varrho\_{n,i}\\ is the design weight
of the held-out (target) curve. Two regimes:

- **Common design** — the operators are estimated once on the training
  block and only the conditioning values vary across the validation set.

- **Independent design** — a rolling origin: for each target curve the
  plug-in estimates are refreshed on all curves observed up to its
  predecessor, with the adaptive bandwidths held fixed at those selected
  once on the initial block (cached in the fit).

Only the \\M \times M\\ system depends on \\\alpha\\; it is solved for
the whole grid from a single eigendecomposition per fold.

## See also

[`blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/blup_fit.md),
[`predict.blup_fit()`](https://hmaissoro.github.io/adaptiveFTS/reference/predict.blup_fit.md).
