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
if (FALSE) { # \dontrun{
# Load the dataset
data(data_far)

# Estimate the optimal bandwidth on a subset of 30 curves
hbest <- get_nw_optimal_bw(data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
                           bw_grid = NULL, nsubset = 30, kernel_name = "epanechnikov")
# Display the optimal bandwidth
hbest
} # }

```
