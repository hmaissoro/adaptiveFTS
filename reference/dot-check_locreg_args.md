# Validate the local regularity arguments forwarded to estimate_locreg_cpp

The adaptive estimators reach the local regularity step from inside C++,
so these four arguments are only checked here, once, before they are
handed down. All four are optional; `NULL` means "use the C++ default",
which is what every estimator did before they became settable.

## Usage

``` r
.check_locreg_args(
  presmooth_bw = NULL,
  Delta = NULL,
  presmooth_bw_grid = NULL,
  presmooth_nsubset = NULL
)
```

## Arguments

- presmooth_bw:

  The presmoothing bandwidth: a positive scalar, a vector of one
  bandwidth per curve, or `NULL` to select it by cross-validation.

- Delta:

  The neighbourhood length: a scalar in (0, 1), or `NULL`.

- presmooth_bw_grid:

  The candidate grid of that cross-validation: a numeric vector of at
  least two positive values, or `NULL`.

- presmooth_nsubset:

  The number of curves used by that cross-validation: a positive
  integer, or `NULL`.

## Value

`NULL`, invisibly. Called for the side effect of erroring out.
