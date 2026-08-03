# Simulate a Functional Moving Average Process of Order One

Simulates `N` curves of a FMA(1) process, \\X_n = \varepsilon_n +
\Psi(\varepsilon\_{n-1})\\, where \\\Psi\\ is the integral operator with
kernel `fma_kernel` and the innovations are multifractional Brownian
motions whose roughness follows `hurst_fun`. It is the moving-average
counterpart of
[`simulate_far`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md)
and shares its arguments and output.

## Usage

``` r
simulate_fma(
  N = 2L,
  lambda = 70L,
  design = "random",
  M_distribution = rpois,
  t_distribution = runif,
  t_common = seq(0.2, 0.8, len = 50),
  hurst_fun = hurst_logistic,
  L = 4,
  intercept_var = 0,
  fma_kernel = function(s, t) 9/4 * exp(-(t + 2 * s)^2),
  fma_mean = function(t) 4 * sin(1.5 * pi * t),
  n_int_grid = 100L,
  n_burnin = 100L,
  remove_burnin = TRUE
)
```

## Arguments

- N:

  `integer`. Number of curves.

- lambda:

  `integer`. Mean of the number of observations per curve.

- design:

  `character`. Type of the design. It is either 'random' or 'common'.

- M_distribution:

  `function`. Distribution of the number of observation points per
  curve. The first argument of the function must correspond to `N` and
  the second to `lambda`. Default `M_distribution = rpois`.

- t_distribution:

  `function (or NULL)`. Observation point distribution if
  `design = 'random'` and `NULL` otherwise.

- t_common:

  `vector (float)`. Observation point vector if `design = 'common'`. If
  `design = 'random'` and if we want to run some tests at a particular
  observation position, this can also be specified.

- hurst_fun:

  `function`. Hurst function. It can be
  [`hurst_arctan`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_arctan.md),
  [`hurst_linear`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_linear.md),
  [`hurst_logistic`](https://hmaissoro.github.io/adaptiveFTS/reference/hurst_logistic.md).

- L:

  `float (positive)`. Hölder constant.

- intercept_var:

  `float (non-negative)`. Variance of a per-curve random intercept added
  to each innovation, expressed relative to the scale of the innovation,
  so that the intercept has variance `L * intercept_var`. It displaces
  each innovation on the ordinate axis without changing its local
  regularity. Passed to
  [`simulate_mfBm`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md),
  whose Details section describes it. Default is `intercept_var = 0`.

- fma_kernel:

  `function`. Kernel function of the operator of the FMA(1).

- fma_mean:

  `function`. Mean function of the FMA(1).

- n_int_grid:

  `integer`. Length of the grid used to approximate the integral.

- n_burnin:

  `integer`. Burnin period of the FMA(1).

- remove_burnin:

  `boolean`. If `TRUE`, n_burnin period is removed.

## Value

A `data.table` containing 3 column :

- id_curve : Index of the curve. It goes from 1 to N.

- tobs : Sampled observation points, for each `id_curve`.

- ttag : Tag on the observation points, for each `id_curve`. It is
  either `tcommon` for the common design grid or `trandom` for the
  random design.

- fma_mean : The mean of the process evaluate at `tobs`, for each
  `id_curve`.

- X : The process observed at tobs, for each `id_curve`.

## Examples

``` r

dt_fma <- simulate_fma(N = 2L, lambda = 70L,
                       design = "random",
                       M_distribution = rpois,
                       t_distribution = runif,
                       t_common = seq(0.2, 0.8, len = 50),
                       hurst_fun = hurst_logistic,
                       L = 4,
                       fma_kernel = function(s,t) 9/4 * exp(- (t + 2 * s) ** 2),
                       fma_mean = function(t) 4 * sin(1.5 * pi * t),
                       n_int_grid = 100L,
                       n_burnin = 100L,
                       remove_burnin = TRUE)
# plot simulated curve
library(ggplot2)

ggplot(data = dt_fma[ttag == "trandom", .("id_curve" = as.factor(id_curve), tobs, X)],
       mapping = aes(x = tobs, y = X, group = id_curve, color = id_curve)) +
  geom_line() +
  scale_colour_grey() +
  theme_minimal()



```
