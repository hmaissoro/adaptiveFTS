# Convert Raw Curve Observations to the Package Data Format

Reshapes raw curve observations into the three-column `data.table` that
every estimator of the package consumes, and checks that the result
meets the assumptions those estimators rely on. All exported estimators
call `format_data` on their `data` argument, so it rarely needs to be
called directly.

## Usage

``` r
format_data(data, idcol = NULL, tcol = "tobs", ycol = "X")
```

## Arguments

- data:

  Raw curve observations, as a `data.table` (or `data.frame`) in long
  format, or as a `list` with one element per curve. See `format_data`
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

## Value

A `data.table` with three columns, sorted by `id_curve` then `tobs`:

- `id_curve`: the curve index, renumbered \\1, \ldots, N\\ in order of
  first appearance in `data`.

- `tobs`: the observation points of each curve.

- `X`: the values observed at each `tobs`.

## Details

Three input layouts are accepted.

- A `data.table` (or `data.frame`) in long format, with one row per
  observation point and at least the columns `idcol`, `tcol` and `ycol`.
  The curve index is repeated once per observation point of that curve.

- A `list` with one element per curve, each element a `data.table` (or
  `data.frame`) holding at least the columns `tcol` and `ycol`. Set
  `idcol = NULL`: the curve index is the position in the list.

- A `list` with one element per curve, each element a `list` of two
  vectors of equal length named `tcol` and `ycol`. Set `idcol = NULL`.

Curves are renumbered \\1, \ldots, N\\ in order of first appearance in
`data`. That order defines the order of the series: the lag-\\\ell\\
estimators pair curve \\n\\ with curve \\n + \ell\\, so the curves must
arrive in chronological order. Within a curve, rows need neither be
contiguous nor sorted; the returned table is always sorted by
`id_curve`, then by `tobs`.

The observation points must lie in \\\[0, 1\]\\, the domain the
estimators assume for the curves; rescale them beforehand if they are
recorded on another scale. Missing values are rejected rather than
dropped, so that the number of observation points per curve is the one
the caller intends.

## Examples

``` r
data("data_far")

# Long format: one row per observation point.
dt <- format_data(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X")
head(dt)
#>    id_curve       tobs        X
#>       <int>      <num>    <num>
#> 1:        1 0.02173214 242.0116
#> 2:        1 0.02273069 241.9948
#> 3:        1 0.05885535 239.9626
#> 4:        1 0.08130179 240.0331
#> 5:        1 0.09080029 240.4602
#> 6:        1 0.09862983 240.7074

# One list element per curve: the curve index is the position in the list.
curves <- split(dt, by = "id_curve", keep.by = FALSE)
head(format_data(data = curves, tcol = "tobs", ycol = "X"))
#>    id_curve       tobs        X
#>       <int>      <num>    <num>
#> 1:        1 0.02173214 242.0116
#> 2:        1 0.02273069 241.9948
#> 3:        1 0.05885535 239.9626
#> 4:        1 0.08130179 240.0331
#> 5:        1 0.09080029 240.4602
#> 6:        1 0.09862983 240.7074
```
