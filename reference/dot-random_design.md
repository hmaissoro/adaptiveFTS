# Generate a random design of the place where the observation is made.

Generate a random design of the place where the observation is made.

## Usage

``` r
.random_design(N, lambda, Mdistribution = rpois, tdistribution = runif, ...)
```

## Arguments

- N:

  `integer`. Number of curves.

- lambda:

  `integer`. Mean of the number of observations per curve.

- Mdistribution:

  `function`. Distribution of the number of observation points per
  curve. The first argument of the function must correspond to `N` and
  the second to `lambda`. Default `Mdistribution = rpois`.

- tdistribution:

  `function`. Distribution of the observation point in the domain.
  Currently only `runif` is accepted.

- ...:

  Additional argument of `tdistribution`.

## Value

A `data.table` containing 3 column :

- id_curve : Index of the curve. It goes from 1 to N.

- Mn : Number of sampled observation location.

- Tn : Sampled observation location.
