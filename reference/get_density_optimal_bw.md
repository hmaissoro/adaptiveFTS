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
