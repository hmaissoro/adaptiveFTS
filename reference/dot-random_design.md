# Generate a random design of the place where the observation is made.

Generate a random design of the place where the observation is made.

## Usage

``` r
.random_design(N, lambda, M_distribution = rpois, t_distribution = runif, ...)
```

## Arguments

- N:

  `integer`. Number of curves.

- lambda:

  `integer`. Mean of the number of observations per curve.

- M_distribution:

  `function`. Distribution of the number of observation points per
  curve. The first argument of the function must correspond to `N` and
  the second to `lambda`. Default `M_distribution = rpois`.

- t_distribution:

  `function`. Distribution of the observation point in the domain.
  Currently only `runif` is accepted.

- ...:

  Additional argument of `t_distribution`.

## Value

A `data.table` containing 3 column :

- id_curve : Index of the curve. It goes from 1 to N.

- Mn : Number of sampled observation location.

- Tn : Sampled observation location.
