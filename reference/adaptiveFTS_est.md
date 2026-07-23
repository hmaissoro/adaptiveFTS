# Adaptive functional time series estimator objects

The adaptive estimators of adaptiveFTS (for example
[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md)
and their risk counterparts) return a
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
carrying an extra S3 class so that dedicated
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[ggplot2::autoplot()](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods are available. The object is still a genuine `data.table`: every
column and every data.table operation behaves exactly as before.

## Details

Each result inherits from the shared parent class `adaptiveFTS_est` and
from a per-estimator subclass (`mean_est`, `mean_risk`, `autocov_est`,
`autocov_risk`, `locreg_est`, `cov_segment_est`, `cov_segment_risk`,
`fts_acf`). A short `adaptive_meta` attribute stores context used by the
methods (design type, kernel, number of curves, lag, ...).

## Subsetting caveat

A data.table `[` subset returns a plain `data.table` — the leading
`adaptiveFTS_est` class is dropped, so
[`summary()`](https://rdrr.io/r/base/summary.html) /
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) should be
called on the estimator's fresh result. Re-tagging after a subset is
possible but not generally needed.

## See also

[`estimate_mean()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_mean.md),
[`estimate_autocov()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_autocov.md),
[`estimate_locreg()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_locreg.md),
[`estimate_cov_segment()`](https://hmaissoro.github.io/adaptiveFTS/reference/estimate_cov_segment.md).
