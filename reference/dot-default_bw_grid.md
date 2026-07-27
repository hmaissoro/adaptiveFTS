# Default candidate bandwidth grid of the adaptive risk functions

A 20-point geometric grid running from \\4(N\lambda)^{-0.9}\\ to
\\4(N\lambda)^{-1/3}\\, where \\N\\ is the number of curves and
\\\lambda\\ the average number of observation points per curve.

## Usage

``` r
.default_bw_grid(data)
```

## Arguments

- data:

  A formatted `data.table`, as returned by
  [`format_data()`](https://hmaissoro.github.io/adaptiveFTS/reference/format_data.md).

## Value

A `numeric` vector of candidate bandwidths.
