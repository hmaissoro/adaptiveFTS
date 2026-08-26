# Save a plot as a (standalone) TikZ/LaTeX figure

Renders a plot to a `.tex` file through
[`tikzDevice::tikz()`](https://rdrr.io/pkg/tikzDevice/man/tikz.html)
and, optionally, compiles it to PDF with `pdflatex`. This is the package
version of the `plot_save_tkz()`/`save_tikz()` helper that is otherwise
copied between analysis projects.

## Usage

``` r
save_plot_tikz(
  plot,
  file,
  width = 7,
  height = 5,
  standalone = TRUE,
  packages = NULL,
  compile = FALSE,
  clean = TRUE,
  ...
)
```

## Arguments

- plot:

  A printable plot object — typically a
  [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html),
  a `ggpubr::ggarrange` arrangement, or any object whose
  [`print()`](https://rdrr.io/r/base/print.html) method draws to the
  current graphics device.

- file:

  `character(1)`. Output path. A `.tex` extension is appended when
  missing, and the parent directory is created if needed.

- width, height:

  `numeric(1)`. Device size in inches. Defaults `7` and `5`.

- standalone:

  `logical(1)`. Passed to
  [`tikzDevice::tikz()`](https://rdrr.io/pkg/tikzDevice/man/tikz.html)
  as `standAlone`. `TRUE` (default) writes a self-contained document
  that `pdflatex` can compile directly.

- packages:

  `character` or `NULL`. Extra LaTeX packages forwarded to
  [`tikzDevice::tikz()`](https://rdrr.io/pkg/tikzDevice/man/tikz.html)
  via its `package` argument, e.g.
  `c(getOption("tikzLatexPackages"), "\\usepackage{amsmath}")`.

- compile:

  `logical(1)`. If `TRUE`, run `pdflatex` on the generated file
  (requires `standalone = TRUE` and `pdflatex` on the `PATH`). Default
  `FALSE`.

- clean:

  `logical(1)`. If `TRUE` (default) and `compile` ran, remove the `.aux`
  and `.log` by-products.

- ...:

  Further arguments passed to
  [`tikzDevice::tikz()`](https://rdrr.io/pkg/tikzDevice/man/tikz.html).

## Value

The path to the written `.tex` file, invisibly.

## See also

[`tikzDevice::tikz()`](https://rdrr.io/pkg/tikzDevice/man/tikz.html),
[adaptiveFTS-autoplot](https://hmaissoro.github.io/adaptiveFTS/reference/adaptiveFTS-autoplot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("tikzDevice", quietly = TRUE)) {
  data("data_far")
  g <- ggplot2::autoplot(estimate_mean(data_far, t = seq(0.1, 0.9, 0.1)))
  save_plot_tikz(g, file = file.path(tempdir(), "mean.tex"),
                 width = 7, height = 5)
}
} # }
```
