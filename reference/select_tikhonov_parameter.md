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
  n_subgrid_bw = 10L,
  n_cv_tikhonov = 30L,
  tikhonov_grid = NULL,
  bw_grid = NULL,
  density_bw = NULL
)
```

## Arguments

- data:

  A `data.table` (or `data.frame`), a `list` of `data.table` (or
  `data.frame`), or a `list` of `list`.

  - If `data.table`: It should contain the raw curve observations in at
    least three columns.

    - `idcol` : The name of the column containing the curve index in the
      sample. Each curve index is repeated according to the number of
      observation points.

    - `tcol` : The name of the column with observation points associated
      with each curve index.

    - `ycol` : The name of the column with observed values at each
      observation point for each curve index.

  - If `list` of `data.table`: In this case, each element in the `list`
    represents the observation data of a curve in the form of a
    `data.table` or `data.frame`. Each `data.table` contains at least
    two columns.

    - `tcol` : The name of the column with observation points for the
      curve.

    - `ycol` : The name of the column with observed values for the
      curve.

  - If `list` of `list`: In this case, `data` is a list where each
    element is the observation data of a curve, given as a `list` of two
    vectors.

    - `tcol` : The vector containing observation points for the curve.

    - `ycol` : The vector containing observed values for the curve.

- idcol:

  `character`. If `data` is given as a `data.table` or `data.frame`,
  this is the name of the column that holds the curve index. Each curve
  index is repeated according to the number of observation points. If
  `data` is a `list` of `data.table` (or `data.frame`) or a `list` of
  `list`, set `idcol = NULL`.

- tcol:

  `character`. The name of the column (or vector) containing the
  observation points for the curves.

- ycol:

  `character`. The name of the column with observed values for the
  curves.

- method:

  Selection method. Currently only `"cv"` (one-step-ahead
  cross-validation) is available.

- kernel_name:

  Kernel name. Default `"epanechnikov"`.

- homoscedastic:

  If `TRUE` (default) a constant noise variance (median of the pointwise
  estimates) is used; otherwise the t-varying estimates.

- n_subgrid_bw:

  Number of points per axis of the coarse sub-grid on which the adaptive
  bandwidths are selected. Default `10`.

- n_cv_tikhonov:

  Number of trailing curves used for the Tikhonov cross-validation when
  `tikhonov` is `NULL`. Default `30`.

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

## Value

A list with:

- `tikhonov_star`: the selected Tikhonov parameter.

- `tikhonov_grid`: the candidate grid.

- `cv_curve`: the mean cross-validation score per candidate.

- `cv_matrix`: the per-fold cross-validation scores (fold by candidate).

- `val_ids`: the ids of the validation curves.

## Details

Each of the last `n_cv_tikhonov` curves is predicted from its immediate
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
