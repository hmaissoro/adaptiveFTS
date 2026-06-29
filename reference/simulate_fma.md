# Functional Moving Average process of order 1 (FMA(1)) simulation

Functional Moving Average process of order 1 (FMA(1)) simulation

## Usage

``` r
simulate_fma(
  N = 2L,
  lambda = 70L,
  tdesign = "random",
  Mdistribution = rpois,
  tdistribution = runif,
  tcommon = seq(0.2, 0.8, len = 50),
  hurst_fun = hurst_logistic,
  L = 4,
  fma_kernel = function(s, t) 9/4 * exp(-(t + 2 * s)^2),
  fma_mean = function(t) 4 * sin(1.5 * pi * t),
  int_grid = 100L,
  burnin = 100L,
  remove_burnin = TRUE
)
```

## Arguments

- N:

  `integer`. Number of curves.

- lambda:

  `integer`. Mean of the number of observations per curve.

- tdesign:

  `character`. Type of the design. It is either 'random' or 'common'.

- Mdistribution:

  `function`. Distribution of the number of observation points per
  curve. The first argument of the function must correspond to `N` and
  the second to `lambda`. Default `Mdistribution = rpois`.

- tdistribution:

  `function (or NULL)`. Observation point distribution if
  `tdesign = 'random'` and `NULL` otherwise.

- tcommon:

  `vector (float)`. Observation point vector if `tdesign = 'common'`. If
  `tdesign = 'random'` and if we want to run some tests at a particular
  observation position, this can also be specified.

- hurst_fun:

  `function`. Hurst function. It can be
  [`hurst_arctan`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
  [`hurst_linear`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md),
  [`hurst_logistic`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md).

- L:

  `float (positive)`. Hölder constant.

- fma_kernel:

  `function`. Kernel function of the operator of the FMA(1).

- fma_mean:

  `function`. Mean function of the FMA(1).

- int_grid:

  `integer`. Length of the grid used to approximate the integral.

- burnin:

  `integer`. Burnin period of the FMA(1).

- remove_burnin:

  `boolean`. If `TRUE`, burnin period is removed.

## Value

A `data.table` containing 3 column :

- id_curve : Index of the curve. It goes from 1 to N.

- tobs : Sampled observation points, for each `id_curve`.

- ttag : Tag on the observations points, for each `id_curve`. It is
  either `tcommon` for common design grid or `tcommon` pour random
  design.

- fma_mean : The mean of the process evaluate at `tobs`, for each
  `id_curve`.

- X : The process observed at tobs, for each `id_curve`.

## Examples

``` r

if (FALSE) { # \dontrun{
dt_fma <- simulate_fma(N = 2L, lambda = 70L,
                       tdesign = "random",
                       Mdistribution = rpois,
                       tdistribution = runif,
                       tcommon = seq(0.2, 0.8, len = 50),
                       hurst_fun = hurst_logistic,
                       L = 4,
                       fma_kernel = function(s,t) 9/4 * exp(- (t + 2 * s) ** 2),
                       fma_mean = function(t) 4 * sin(1.5 * pi * t),
                       int_grid = 100L,
                       burnin = 100L,
                       remove_burnin = TRUE)
# plot simulated curve
library(ggplot2)

ggplot(data = dt_fma[ttag == "trandom", .("id_curve" = as.factor(id_curve), tobs, X)],
       mapping = aes(x = tobs, y = X, group = id_curve, color = id_curve)) +
  geom_line() +
  scale_colour_grey() +
  theme_minimal()

} # }
```
