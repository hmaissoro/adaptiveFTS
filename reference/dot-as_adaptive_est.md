# Tag an estimator result with its adaptive-estimator S3 class

Tag an estimator result with its adaptive-estimator S3 class

## Usage

``` r
.as_adaptive_est(dt, subclass, meta = list())
```

## Arguments

- dt:

  A
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  (modified in place by reference).

- subclass:

  `character(1)`. The per-estimator subclass, e.g. `"mean_est"`.

- meta:

  `list`. Context stored in the `adaptive_meta` attribute and read back
  by the `summary`/`plot` methods.

## Value

`dt`, invisibly re-classed and carrying the `adaptive_meta` attribute.
