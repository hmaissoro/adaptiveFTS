# Local Regularity Parameters Estimation

This function estimates the local regularity parameters \\H_t\\ and
\\L_t^2\\ defined in Section 3 of Maissoro et al. (2024) .

## Usage

``` r
estimate_locreg(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = 1/2,
  Delta = NULL,
  h = NULL,
  kernel_name = "epanechnikov",
  center = TRUE
)
```

## Arguments

- data:

  A `data.table` (or `data.frame`), a `list` of `data.table` (or
  `data.frame`), or a `list` of `list`.

  - If `data.table`: It should contain the raw curve observations in at
    least three columns.

    - `idcol` : The name of the column containing the curve index in the
      sample. Each curve index is repeated according to the number of
      observation points.

    - `tcol` : The name of the column with observation points associated
      with each curve index.

    - `ycol` : The name of the column with observed values at each
      observation point for each curve index.

  - If `list` of `data.table`: In this case, each element in the `list`
    represents the observation data of a curve in the form of a
    `data.table` or `data.frame`. Each `data.table` contains at least
    two columns.

    - `tcol` : The name of the column with observation points for the
      curve.

    - `ycol` : The name of the column with observed values for the
      curve.

  - If `list` of `list`: In this case, `data` is a list where each
    element is the observation data of a curve, given as a `list` of two
    vectors.

    - `tcol` : The vector containing observation points for the curve.

    - `ycol` : The vector containing observed values for the curve.

- idcol:

  `character`. If `data` is given as a `data.table` or `data.frame`,
  this is the name of the column that holds the curve index. Each curve
  index is repeated according to the number of observation points. If
  `data` is a `list` of `data.table` (or `data.frame`) or a `list` of
  `list`, set `idcol = NULL`.

- tcol:

  `character`. The name of the column (or vector) containing the
  observation points for the curves.

- ycol:

  `character`. The name of the column with observed values for the
  curves.

- t:

  `vector (numeric)`. Observation points at which to estimate the local
  regularity parameters of the underlying process.

- Delta:

  `numeric (positive)`. The length of the neighborhood around each point
  in `t` for local regularity estimation. Default is `Delta = NULL`, in
  which case it will be estimated from the data.

- h:

  `numeric (positive vector or scalar)`. Bandwidth parameter for the
  Nadaraya-Watson estimator used in local regularity estimation. Default
  is `h = NULL`, allowing it to be determined by cross-validation on a
  subset of curves. If `h` is a scalar, the same bandwidth is applied
  across all curves. If it is a vector, it must match the number of
  curves in `data`, with each element corresponding to a curve.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

- center:

  `logical`. If `TRUE`, centers the curves before regularity estimation.

## Value

A `data.table` with columns:

- `t` : The observation points around which local regularity is
  estimated.

- `locreg_bw` : The presmoothing bandwidth used for estimation.

- `Delta` : The neighborhood length around each `t` used in local
  regularity estimation.

- `Nused` : The number of curves contributing non-degenerate estimates
  around each `t`.

- `Ht` : Local exponent estimates, denoted by \\H_t\\.

- `Lt` : Hölder constant estimates, corresponding to \\L_t^2\\.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.

## See also

[`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[`estimate_nw_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw_bw.md),
[`simulate_far()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_far.md),
etc.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Load data
 data("data_far")

 # Define observation points for local regularity estimation
 t0 <- seq(0.2, 0.8, length.out = 8)

 # Estimate local regularity parameters
 dt_locreg <- estimate_locreg(data = data_far,
                              idcol = "id_curve",
                              tcol = "tobs",
                              ycol = "X",
                              t = t0,
                              Delta = NULL,
                              h = NULL,
                              kernel_name = "epanechnikov",
                              center = TRUE)

 dt_locreg
} # }
```
