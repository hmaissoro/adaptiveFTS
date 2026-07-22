# Mean at new locations using cached bandwidths (C++ core)

Reuses the adaptive mean bandwidths cached in a `blup_fit` object.
Called by
[`select_tikhonov_parameter()`](https://hmaissoro.github.io/adaptiveFTS/reference/select_tikhonov_parameter.md);
not intended to be used directly.

## Usage

``` r
blup_mean_at_cpp(data, opt_mean, t, kernel_name)
```

## Arguments

- data:

  A DataFrame with columns `id_curve`, `tobs`, `X`.

- opt_mean:

  Cached mean adaptive-bandwidth matrix (`t`, `optbw`).

- t:

  Evaluation locations (assumed sorted).

- kernel_name:

  Kernel name.

## Value

The mean estimates at `t`.
