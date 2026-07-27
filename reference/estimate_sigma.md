# Estimate the the standard deviation of the observation error

This function estimates the the standard deviation of the observation
error using the estimator proposed by Maissoro, Patilea and Vimond
(2025).

## Usage

``` r
estimate_sigma(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4)
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

  `vector (numeric)`. Observation points at which we want to estimate
  the standard deviation of the error.

## Value

A data.table with two columns: `t` and `sig` corresponding to the
estimated standard deviation.

## References

Maissoro, H., Patilea, V. and Vimond, M. (2025). Adaptive Estimation for
Weakly Dependent Functional Time Series. *Journal of Time Series
Analysis*. [doi:10.1111/jtsa.70006](https://doi.org/10.1111/jtsa.70006)

## Examples

``` r
# Load data
data("data_far")

# Estimate the standar-deviation of the error term
estimate_sigma(data = data_far, t = c(1/4, 1/2, 3/4))
#>        t       sig
#>    <num>     <num>
#> 1:  0.25 0.4468025
#> 2:  0.50 0.3524165
#> 3:  0.75 0.4441715



```
