# Least-squares cross-validation score of the design-density bandwidth

Least-squares cross-validation score of the design-density bandwidth

## Usage

``` r
.density_cv_score(hval, x, kern, lower, upper)
```

## Arguments

- hval:

  Bandwidth.

- x:

  Observation times of a single curve.

- kern:

  Kernel function.

- lower, upper:

  Bounds of the domain.

## Value

The LSCV score at `hval`.
