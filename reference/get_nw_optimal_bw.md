# Estimate Optimal Bandwidth for Nadaraya-Watson Estimator on a Subset of Curves

This function estimates the median optimal bandwidth for the
Nadaraya-Watson kernel estimator using cross-validation on a subset of
curves.

## Usage

``` r
get_nw_optimal_bw(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  bw_grid = NULL,
  nsubset = NULL,
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

- bw_grid:

  `vector (numeric)`. A grid of candidate bandwidth values for
  cross-validation. Default is `bw_grid = NULL`, which sets an
  exponential grid based on the average number of observation points per
  curve.

- nsubset:

  `integer (positive)`. The number of curves to randomly and uniformly
  select for bandwidth optimization. Default is `nsubset = NULL`, in
  which case an optimal bandwidth is calculated for each curve.

- kernel_name:

  `string`. A string specifying the name of the kernel function to use,
  with "epanechnikov" as the default. Supported kernels: "epanechnikov",
  "biweight", "triweight", "tricube", "triangular", and "uniform".

## Value

A `numeric` scalar representing the estimated optimal bandwidth as the
median of the best bandwidths from the subset of curves.

## Details

This function performs cross-validation to determine the optimal
bandwidth for each curve in a specified subset. It returns the median of
these best bandwidths as the final estimate, providing a representative
bandwidth that can generalize across curves.

## Examples

``` r
# Load the dataset
data(data_far)

# Estimate the optimal bandwidth on a subset of 30 curves
hbest <- get_nw_optimal_bw(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
                           bw_grid = NULL, nsubset = 30, kernel_name = "epanechnikov")
# Display the optimal bandwidth
hbest
#> [1] 0.02581033


```
