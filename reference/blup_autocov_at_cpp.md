# (Auto)covariance block at new locations using cached bandwidths (C++ core)

Reuses the adaptive (auto)covariance bandwidths cached in a `blup_fit`
object. Called by
[`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md);
not intended to be used directly.

## Usage

``` r
blup_autocov_at_cpp(data, opt_bw, s, t, lag, correct_diagonal, kernel_name)
```

## Arguments

- data:

  A DataFrame with columns `id_curve`, `tobs`, `X`.

- opt_bw:

  Cached (auto)covariance bandwidth matrix (`s`, `t`, `bw_s`, `bw_t`).

- s, t:

  Evaluation locations (rows indexed by `s`, columns by `t`).

- lag:

  0 for the covariance, 1 for the lag-1 autocovariance.

- correct_diagonal:

  Whether to correct the covariance diagonal.

- kernel_name:

  Kernel name.

## Value

A `length(s)` by `length(t)` matrix of \\\hat c\_{lag}(s_i, t_j)\\.
