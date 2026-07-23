# Retrieve the metadata attached to an adaptive-estimator object

Retrieve the metadata attached to an adaptive-estimator object

## Usage

``` r
.est_meta(x, field = NULL)
```

## Arguments

- x:

  An `adaptiveFTS_est` object.

- field:

  `character(1)` or `NULL`. If `NULL` (default), the whole metadata
  list; otherwise the named element (`NULL` when absent).

## Value

The metadata list or the requested element.
