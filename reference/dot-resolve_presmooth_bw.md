# Validate or select the presmoothing bandwidth of the empirical estimators

Validate or select the presmoothing bandwidth of the empirical
estimators

## Usage

``` r
.resolve_presmooth_bw(presmooth_bw, data, N, kernel_name)
```

## Arguments

- presmooth_bw:

  The user-supplied bandwidth: a scalar, a vector of length `N`, or
  `NULL` to select it by cross-validation.

- data:

  A formatted `data.table`, as returned by
  [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md).

- N:

  `integer(1)`. The number of curves in `data`.

- kernel_name:

  `character(1)`. The kernel used for the presmoothing.

## Value

A `numeric` vector of length `N`, one bandwidth per curve.
