# Estimate empirical \\p\\-th order moment of \\X(t)\\.

This function estimates the \\p\\-th order moment of \\X(t)\\, used in
the empirical study section of the papers Maissoro et al. (2024) and
Maissoro et al. (2025) .

## Usage

``` r
estimate_empirical_mom(
  data,
  idcol = NULL,
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  mom_order = 1,
  h = NULL,
  center = TRUE,
  kernel_name = "epanechnikov"
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

  `vector (numeric)`. Observation points at which the \\p\\-th order
  moment of \\X(t)\\ is estimated. Each element should be a value
  between 0 and 1.

- mom_order:

  `numeric (positive scalar)`. The order of the moment to be computed
  (e.g., 1 for mean, 2 for variance).

- h:

  `numeric (positive vector or scalar)`. The smoothing bandwidth
  parameter. Default `h = NULL`, in which case the bandwidth will be
  estimated by Cross-Validation on a subset of curves. If `h` is a
  scalar, then all curves will be smoothed with the same bandwidth. If
  `h` is a vector, its length must equal the number of curves in `data`,
  with each element corresponding to a curve in the same order as in
  `data`.

- center:

  `logical`. If `TRUE`, then the \\p\\-th order moment of the centered
  \\X(t)\\ is estimated. Default is `TRUE`.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` with three columns: `t`, `mom_order`, and `mom_estimate`
corresponding to the estimated \\p\\-th order moment of \\X(t)\\ at each
time point specified in `t`.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.  
  
Maissoro H, Patilea V, Vimond M (2025). “Adaptive prediction for
Functional Times Series.” *arXiv preprint arXiv:2501.xxxxx*.

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data("data_far")  # Replace with actual data containing observed curves

# Define parameters
observation_points <- c(0.25, 0.5, 0.75)  # Points at which to estimate moments
moment_order <- 2                         # Example: 2nd order moment (variance)
bandwidth <- 0.1                          # Smoothing parameter; can be NULL for CV-estimated

# Estimate the 2nd order moment (variance) at specified observation points
moment_estimates <- estimate_empirical_mom(
  data = data_far,
  t = observation_points,
  mom_order = moment_order,
  h = NULL,
  center = TRUE,
  kernel_name = "epanechnikov"
)

# View the result
print(moment_estimates)
} # }
```
