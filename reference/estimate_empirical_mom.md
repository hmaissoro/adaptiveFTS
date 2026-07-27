# Estimate empirical \\p\\-th order moment of \\X(t)\\.

This function estimates the \\p\\-th order moment of \\X(t)\\, used in
the empirical study section of the papers Maissoro, Patilea and Vimond
(2025) and Maissoro, Patilea and Vimond (2026).

## Usage

``` r
estimate_empirical_mom(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  mom_order = 1,
  presmooth_bw = NULL,
  center = TRUE,
  kernel_name = "epanechnikov"
)
```

## Arguments

- data:

  Raw curve observations, as a `data.table` (or `data.frame`) in long
  format, or as a `list` with one element per curve. See
  [`format_data`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md)
  for the accepted layouts and for the `id_curve` / `tobs` / `X` columns
  they are converted to.

- idcol:

  `character(1)` or `NULL`. Name of the column holding the curve index
  when `data` is a single table. Must be `NULL` when `data` is a list of
  curves.

- tcol:

  `character(1)`. Name of the column (or vector) holding the observation
  points of the curves.

- ycol:

  `character(1)`. Name of the column (or vector) holding the values
  observed at those points.

- t:

  `vector (numeric)`. Observation points at which the \\p\\-th order
  moment of \\X(t)\\ is estimated. Each element should be a value
  between 0 and 1.

- mom_order:

  `numeric (positive scalar)`. The order of the moment to be computed
  (e.g., 1 for mean, 2 for variance).

- presmooth_bw:

  `numeric (positive vector or scalar)`. Bandwidth used to presmooth
  each curve before the estimation. A scalar applies the same bandwidth
  to every curve; a vector must hold one bandwidth per curve, in the
  order the curves appear in `data`. Default `NULL` selects it by
  cross-validation, see
  [get_nw_optimal_bw](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

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

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

Maissoro, H., Patilea, V. and Vimond, M. (2026). Adaptive Prediction for
Functional Time Series. *arXiv preprint* arXiv:2609.xxxxx.

## See also

[`get_nw_optimal_bw()`](https://hmaissoro.github.io/adaptiveFTS/reference/get_nw_optimal_bw.md).

## Examples

``` r
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
  presmooth_bw = NULL,
  center = TRUE,
  kernel_name = "epanechnikov"
)

# View the result
print(moment_estimates)
#>        t mom_order mom_estimate
#>    <num>     <num>        <num>
#> 1:  0.25         2     9.133629
#> 2:  0.50         2     9.740608
#> 3:  0.75         2    10.465018

```
