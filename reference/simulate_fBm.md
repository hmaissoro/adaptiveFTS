# Draw a fractional Brownian motion sample path.

Draw a fractional Brownian motion sample path.

## Usage

``` r
simulate_fBm(
  t = seq(0.2, 0.8, len = 20),
  hurst = 0.6,
  L = 1,
  intercept_var = 0,
  tied = TRUE
)
```

## Arguments

- t:

  `vector (float)`. Grid of points between 0 and 1 where we want to
  generate the sample path.

- hurst:

  `float (positive)`. The Hurst exponent scalar value between 0 and 1.

- L:

  `float (positive)`. Hölder constant.

- intercept_var:

  `float (non-negative)`. Variance of a per-curve random intercept added
  to the sample path, expressed relative to the scale of the process, so
  that the intercept has variance `L * intercept_var`. It displaces the
  path on the ordinate axis without changing its local regularity.
  Default is `intercept_var = 0`, which adds no intercept. Ignored when
  `tied = TRUE`. See the Details section.

- tied:

  `boolean`. If `TRUE`, the sample path is tied-down.

## Value

A `data.table` containing 2 column : `t` and `fBm`, the sample path.

## Details

Let \\\xi\\ denote the standardised fractional Brownian motion with
exponent `hurst`, that is the centred Gaussian process with covariance
\\(u^{2H} + v^{2H} - \|u - v\|^{2H}) / 2\\. Its variance is
\\Var(\xi(t)) = t^{2H}\\, so that \\Var(\xi(1)) = 1\\, and its
increments satisfy \\E\[(\xi(t + \delta) - \xi(t))^2\] = \delta^{2H}\\.
The returned sample path is \\\sqrt{L} \\ \xi(t)\\ when
`intercept_var = 0` and `tied = FALSE`, and \\\sqrt{L} \\ (\xi(t) + Z)\\
otherwise, where \\Z\\ is a centred Gaussian variable of variance
`intercept_var`, drawn independently of \\\xi\\ and constant in `t`. See
the Details section of
[`simulate_mfBm`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)
for the role of the intercept.

This is the same process as
[`simulate_mfBm`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md)
given a Hurst function constant at `hurst`.

## See also

[`simulate_mfBm()`](https://hmaissoro.github.io/adaptiveFTS/reference/simulate_mfBm.md).

## Examples

``` r

t0 <- seq(0.2, 0.8, len = 20)
dt_fBm <- simulate_fBm(t = t0, hurst = 0.6, L = 1, tied = TRUE)
plot(x = dt_fBm$t, y = dt_fBm$fBm, type = "l", col = "red")

```
