# Project a matrix onto the PSD cone (C++ core)

Symmetrises `M` and floors its eigenvalues at zero. Called by
[`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md);
not intended to be used directly.

## Usage

``` r
psd_project_cpp(M)
```

## Arguments

- M:

  A numeric matrix.

## Value

The nearest positive-semidefinite matrix.
