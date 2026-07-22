# Leave-one-out Parzen-Rosenblatt estimate at the observation points

Leave-one-out Parzen-Rosenblatt estimate at the observation points

## Usage

``` r
.density_loo_estimate(x, hval, kern)
```

## Arguments

- x:

  Observation times of a single curve.

- hval:

  Bandwidth.

- kern:

  Kernel function.

## Value

The vector \\\widehat g\_{h}^{(-i)}(x_i)\\, in the order of `x`.
