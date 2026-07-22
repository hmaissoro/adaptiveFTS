# Leave-one-out Parzen-Rosenblatt density estimator

Computes the leave-one-out Parzen-Rosenblatt estimator of the design
density evaluated at the observation points themselves, as required by
the independent-design weights of the adaptive BLUP.

## Usage

``` r
estimate_density(
  x,
  h = NULL,
  bw_grid = NULL,
  kernel_name = "epanechnikov",
  lower = 0,
  upper = 1
)
```

## Arguments

- x:

  Numeric vector of observation times of a *single* curve
  \\T\_{n_0,\cdot}\\ (not pooled across curves).

- h:

  Optional numeric scalar bandwidth. If supplied, the least-squares
  cross-validation is skipped and this bandwidth is used directly. This
  is the path used by the adaptive BLUP, where a single bandwidth
  selected once on a subset of curves (see
  [`get_density_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_density_optimal_bw.md))
  is reused for every curve. Default is `NULL`, in which case the
  bandwidth is selected by LSCV.

- bw_grid:

  Numeric vector of candidate bandwidths for the LSCV search (used only
  when `h` is `NULL`). If `NULL` (default), a fixed log-spaced grid
  close to `exp(seq(log(0.01), log(0.3), length.out = 30))` is used. For
  Monte Carlo studies, pass a fixed grid so the cross-validation is
  comparable across replications.

- kernel_name:

  Kernel name: "epanechnikov" (default), "biweight", "triweight",
  "tricube", "triangular", or "uniform".

- lower, upper:

  Bounds of the domain \\I\\; default (0, 1\].

## Value

A list with:

- `h_star`: the bandwidth used.

- `bw_grid`: the candidate bandwidth grid (as supplied or defaulted).

- `cv_curve`: the LSCV score per candidate bandwidth, or `NULL` when `h`
  is supplied.

- `kernel_name`: the kernel used.

- `estimate`: the vector \\\widehat g\_{n_0,i}(T\_{n_0,i})\\, one value
  per observation point, in the order of `x`.

## Details

For a single curve with observation times \\T\_{n_0,i}\\, \\1 \le i \le
M\_{n_0}\\, the estimator at the observation points is \$\$\widehat
g\_{n_0,i}(T\_{n_0,i}) = \frac{1}{(M\_{n_0}-1)\\h} \sum\_{j \ne i}
K\\\left(\frac{T\_{n_0,i} - T\_{n_0,j}}{h}\right),\$\$ feeding the
independent-design weights \\\varrho\_{n_0,i} = \\M\_{n_0}\\\widehat
g\_{n_0,i}(T\_{n_0,i})\\^{-1}\\. The bandwidth \\h\\ is either supplied
directly (via `h`) or selected by least-squares (Rudemo-Bowman)
cross-validation over `bw_grid`, which minimises an unbiased estimate
(up to an \\h\\-independent constant) of the density MISE: \$\$CV(h) =
\int \widehat g_h(x)^2 dx - \frac{2}{M} \sum_i \widehat
g_h^{(-i)}(T\_{n_0,i}).\$\$
