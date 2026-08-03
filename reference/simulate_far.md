# Simulate a Functional Autoregressive Process of Order One

Simulates `N` curves of a FAR(1) process, \\X_n = \Psi(X\_{n-1}) +
\varepsilon_n\\, where \\\Psi\\ is the integral operator with kernel
`far_kernel` and the innovations are multifractional Brownian motions
whose roughness follows `hurst_fun`. Each curve is observed at random or
common design points, and the returned sample is the one this package's
estimators consume.

## Usage

``` r
simulate_far(
  N = 2L,
  lambda = 70L,
  design = "random",
  M_distribution = rpois,
  t_distribution = runif,
  t_common = seq(0.2, 0.8, len = 50),
  hurst_fun = hurst_logistic,
  L2 = 4,
  intercept_var = 0,
  far_kernel = function(s, t) 9/4 * exp(-(t + 2 * s)^2),
  far_mean = function(t) 4 * sin(1.5 * pi * t),
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

- L2:

  `float (positive)`. Squared Hölder constant \\L_t^2\\ of the
  innovation, the quantity the package's estimators report in their
  `Lt2`/`Ls2` columns. The Hölder constant itself is \\L_t =
  \sqrt{L2}\\: the innovation increments satisfy \\E\[(\varepsilon(t +
  \delta) - \varepsilon(t))^2\] = L_t^2 \\ \delta^{2 H_t}\\.

- intercept_var:

  `float (non-negative)`. Variance of a per-curve random intercept added
  to each innovation, expressed relative to the scale of the innovation,
  so that the intercept has variance `L2 * intercept_var`. It displaces
  each innovation on the ordinate axis without changing its local
  regularity. Passed to
  [`simulate_mfBm`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md),
  whose Details section describes it. Default is `intercept_var = 0`.

- far_kernel:

  `function`. Kernel function of the operator of the FAR(1).

- far_mean:

  `function`. Mean function of the FAR(1).

- n_int_grid:

  `integer`. Length of the grid used to approximate the integral.

- n_burnin:

  `integer`. Burnin period of the FAR(1).

- remove_burnin:

  `boolean`. If `TRUE`, n_burnin period is removed.

## Value

A `data.table` containing 3 column :

- id_curve : Index of the curve. It goes from 1 to N.

- tobs : Sampled observation points, for each `id_curve`.

- ttag : Tag on the observation points, for each `id_curve`. It is
  either `tcommon` for the common design grid or `trandom` for the
  random design.

- far_mean : The mean of the process evaluate at `tobs`, for each
  `id_curve`.

- X : The process observed at tobs, for each `id_curve`.

## Details

The process is built on a regular integration grid of `n_int_grid`
points and iterated for `n_burnin` steps before the `N` curves are kept,
so that the returned sample is (close to) stationary; set
`remove_burnin = FALSE` to keep the burn-in curves as well. Each curve
is then observed at `M` points, with `M` drawn from `M_distribution`: at
random locations drawn from `t_distribution` when `design = "random"`,
or at the shared grid `t_common` when `design = "common"`.

## Examples

``` r

dt_far <- simulate_far(N = 2L, lambda = 70L,
                       design = "random",
                       M_distribution = rpois,
                       t_distribution = runif,
                       t_common = seq(0.2, 0.8, len = 50),
                       hurst_fun = hurst_logistic,
                       L2 = 4,
                       far_kernel = function(s,t) 9/4 * exp(- (t + 2 * s) ** 2),
                       far_mean = function(t) 4 * sin(1.5 * pi * t),
                       n_int_grid = 100L,
                       n_burnin = 100L,
                       remove_burnin = TRUE)

# Give each innovation a random intercept, so the curves do not share an origin.
dt_far_shifted <- simulate_far(N = 3L, lambda = 40L, design = "random",
                               t_common = NULL, L2 = 4, intercept_var = 0.05,
                               n_int_grid = 60L, n_burnin = 40L)


```
