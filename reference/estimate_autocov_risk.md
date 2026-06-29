# Estimate the Risk of the Covariance or Autocovariance Function

This function estimates the risk function of the adaptive lag-\\\ell\\
autocovariance function estimator, where \\\ell\\ = 0, 1, ..., using one
bandwidth parameter as proposed in Maissoro et al. (2024) or two
bandwidth parameters as proposed in Maissoro et al. (2025) .

## Usage

``` r
estimate_autocov_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  bw_grid = NULL,
  use_same_bw = FALSE,
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

- s:

  `vector (numeric)`. The first argument of the autocovariance function,
  corresponding to observation points `s` in the pair (`s`, `t`). Must
  be of the same length as `t`.

- t:

  `vector (numeric)`. The second argument of the autocovariance
  function, corresponding to observation points `t` in the pair (`s`,
  `t`). Must be of the same length as `s`.

- lag:

  `integer (positive integer)`. The lag of the autocovariance.

- bw_grid:

  `vector (numeric)`. Bandwidth grid for selecting the optimal smoothing
  parameter for each pair (`s`, `t`). Defaults to `NULL`, which
  generates an exponential grid of \\N \lambda\\.

- use_same_bw:

  `logical`. Indicates whether the same bandwidth should be used for
  both `s` and `t`. Defaults to `FALSE`.

- center:

  `logical`. If `TRUE`, centers the data before estimation. Default is
  `TRUE`.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` containing the following columns:

- `s` : The first argument of the autocovariance function.

- `t` : The second argument of the autocovariance function.

- `hs` : Candidate bandwidth for `s`. If `use_same_bw = TRUE`, `hs` and
  `ht` will contain the same values.

- `ht` : Candidate bandwidth for `t`.

- `PNl` : Number of curves used in the autocovariance estimation at
  (`s`, `t`). Corresponds to \\P\_{N,\ell}(s,t;h_s, h_t)\\.

- `locreg_bw` : Bandwidth used for estimating local regularity
  parameters.

- `Hs` : Local exponent estimates for `s`, denoted as \\H_s\\.

- `Ls` : Estimates of the Hölder constant for `s`, corresponding to
  \\L_s^2\\.

- `Ht` : Local exponent estimates for `t`, denoted as \\H_t\\.

- `Lt` : Estimates of the Hölder constant for `t`, corresponding to
  \\L_t^2\\.

- `bias_term` : Bias term of the risk function.

- `variance_term` : Variance term of the risk function.

- `dependence_term` : Dependence term of the risk function.

- `autocov_risk` : Estimated risk of the covariance/autocovariance
  function.

## References

Maissoro H, Patilea V, Vimond M (2024). “Adaptive estimation for Weakly
Dependent Functional Times Series.” *arXiv preprint arXiv:2403.13706*.  
  
Maissoro H, Patilea V, Vimond M (2025). “Adaptive prediction for
Functional Times Series.” *arXiv preprint arXiv:2501.xxxxx*.

## See also

[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_sigma()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_sigma.md),
[`estimate_nw()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_nw.md),
[`estimate_empirical_XsXt_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_empirical_XsXt_autocov.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load data
data("data_far")

# Estimate risk function with same bandwidth for s and t
dt_autocov_risk <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
  bw_grid = NULL, use_same_bw = TRUE, center = TRUE,
  kernel_name = "epanechnikov"
)

# Visualize mean risk function for different (s, t) pairs
dt_dcast <- data.table::dcast(data = dt_autocov_risk,
                              formula = hs ~ s + t,
                              value.var = "autocov_risk")
manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.2, 0.25)" = `0.2_0.25`)],
                      main = "lag = 1 - (s, t) = (0.2, 0.25)",
                      xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.4, 0.5)" = `0.4_0.5`)],
                      main = "lag = 1 - (s, t) = (0.4, 0.5)",
                      xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.8, 0.75)" = `0.8_0.75`)],
                      main = "lag = 1 - (s, t) = (0.8, 0.75)",
                      xlab = "h", ylab = "Risk Function")
  ),
  nrow = 3
)

# Estimate risk function with separate bandwidths for s and t
dt_autocov_risk_2bw <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
  bw_grid = NULL, use_same_bw = FALSE, center = TRUE,
  kernel_name = "epanechnikov"
)

# Display selected columns of the results
dt_autocov_risk_2bw[, .(s, t, lag, hs, ht, autocov_risk)]
} # }
```
