# Estimate the Covariance or Autocovariance Function

This function estimates the adaptive lag-\\\ell\\ autocovariance
function, where \\\ell\\ = 0, 1, ..., using one bandwidth parameter as
proposed in Maissoro et al. (2024) or two bandwidth parameters as
proposed in Maissoro et al. (2025) .

## Usage

``` r
estimate_autocov(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  lag = 1,
  optbw_s = NULL,
  optbw_t = NULL,
  bw_grid = NULL,
  use_same_bw = FALSE,
  center = TRUE,
  correct_diagonal = TRUE,
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

- optbw_s:

  `vector (numeric)`. The optimal bandwidths for `s`. Default is `NULL`.

- optbw_t:

  `vector (numeric)`. The optimal bandwidths for `t`. Default is `NULL`.

- bw_grid:

  `vector (numeric)`. Bandwidth grid for selecting the optimal smoothing
  parameter for each pair (`s`, `t`). Defaults to `NULL`, which
  generates an exponential grid of \\N \lambda\\.

- use_same_bw:

  `logical`. Indicates whether the same bandwidth should be used for
  both `s` and `t`. Defaults to `FALSE`.

- center:

  `logical (TRUE or FALSE)`. Default `center = TRUE` and so the curves
  are centred when the autocovariance is estimated:
  \\\mathbb{E}(X_0(s) - \mu(s))(X\_{\ell}(t) - \mu(t))\\. Otherwise, the
  two parts \\\mathbb{E}X_0(s)X\_{\ell}(t)\\ and \\\mu(s)\mu(t)\\ will
  be estimated separately. The first part with a bandwidth obtained with
  [estimate_autocov_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md)
  and the second part with a bandwidth obtained with
  [estimate_mean_risk](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean_risk.md).

- correct_diagonal:

  `logical (TRUE or FALSE)`. Indicates whether the diagonal of the
  covariance should be corrected when `lag=0`.

- kernel_name:

  `string`. Specifies the kernel function for estimation; default is
  "epanechnikov". Supported kernels include: "epanechnikov", "biweight",
  "triweight", "tricube", "triangular", and "uniform".

## Value

A `data.table` containing the following columns.

- s : The first argument of the autocovariance function.

- t : The second argument of the autocovariance function.

- optbw_s : The optimal bandwidth for the first argument of the
  autocovariance function. If `use_same_bw = TRUE`, the same bandwidth
  candidate is used for `s` and for `t`, so the 3rd and 4th columns
  contain the same values.

- optbw_t : The optimal bandwidth for the second argument of the
  autocovariance function.

- Hs : The estimates of the local exponent for each `s`. It corresponds
  to \\H_s\\.

- Ls : The estimates of the Hölder constant for each `s`. It corresponds
  to \\L_s^2\\.

- Ht : The estimates of the local exponent for each `t`. It corresponds
  to \\H_t\\.

- Lt : The estimates of the Hölder constant for each `t`. It corresponds
  to \\L_t^2\\.

- PNs : The number of curves used to estimate the mean at `s`. It
  corresponds to \\P_N(s;h)\\.

- muhat_s : The estimates of the mean at `s`. It corresponds to
  \\\widehat{\mu}\_N(s;h)\\.

- PNt : The number of curves used to estimate the mean at `t`. It
  corresponds to \\P_N(t;h)\\.

- muhat_t : The estimates of the mean at `t`. It corresponds to
  \\\widehat{\mu}\_N(t;h)\\.

- PNl : The number of curves used to estimate the autocovariance at
  `(s,t)`. It corresponds to \\P\_{N,\ell}(s,t;h_s, h_t)\\.

- autocov : The estimates of the covariance/autocovariance.

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
[`estimate_autocov_risk()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov_risk.md).

## Examples

``` r
if (FALSE) { # \dontrun{
#' # Load data
data("data_far")

# Estimate adaptive lag-1 autocovariance
dt_autocov <- estimate_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 1,
  optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
  use_same_bw = FALSE, center = TRUE, correct_diagonal = FALSE,
  kernel_name = "epanechnikov")
dt_autocov[, .(s, t, lag, PNl, autocov)]

# Estimate adaptive covariance

dt_cov <- estimate_autocov(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4), lag = 0,
  optbw_s = NULL, optbw_t = NULL, bw_grid = NULL,
  use_same_bw = FALSE, center = TRUE, correct_diagonal = TRUE,
  kernel_name = "epanechnikov")
dt_cov[, .(s, t, lag, PNl, autocov)]
} # }
```
