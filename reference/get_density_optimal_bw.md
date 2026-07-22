# Select the design-density bandwidth on a subset of curves

Selects a single Parzen-Rosenblatt bandwidth for the design density by
least-squares cross-validation on a subset of curves, mirroring
[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Usage

``` r
get_density_optimal_bw(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  nsubset = NULL,
  bw_grid = NULL,
  kernel_name = "epanechnikov",
  lower = 0,
  upper = 1
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

- nsubset:

  `integer (positive)`. The number of curves to randomly and uniformly
  select for bandwidth selection. Default is `NULL`, in which case every
  curve is used.

- bw_grid:

  Numeric vector of candidate bandwidths passed to
  [`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md).
  Default is `NULL` (a fixed log-spaced grid). For reproducible Monte
  Carlo studies, pass a fixed grid.

- kernel_name:

  Kernel name: "epanechnikov" (default), "biweight", "triweight",
  "tricube", "triangular", or "uniform".

- lower, upper:

  Bounds of the domain \\I\\; default (0, 1\].

## Value

A `numeric` scalar: the median of the per-curve LSCV-optimal bandwidths
over the subset.

## Details

For each sampled curve the bandwidth minimising the LSCV score of
[`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md)
over `bw_grid` is computed, and the median of these per-curve optima is
returned. The adaptive BLUP calls this once and then reuses the returned
bandwidth for every curve (passing it as the `h` argument of
[`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md)),
so the design weights are consistent across curves and the
cross-validation is not repeated on every call.

## See also

[`estimate_density()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_density.md),
[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).
