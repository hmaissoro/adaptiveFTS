#' Save a plot as a (standalone) TikZ/LaTeX figure
#'
#' Renders a plot to a `.tex` file through [tikzDevice::tikz()] and, optionally,
#' compiles it to PDF with `pdflatex`. This is the package version of the
#' `plot_save_tkz()`/`save_tikz()` helper that is otherwise copied between
#' analysis projects.
#'
#' @param plot A printable plot object — typically a [ggplot2::ggplot], a
#'   `ggpubr::ggarrange` arrangement, or any object whose `print()` method draws
#'   to the current graphics device.
#' @param file \code{character(1)}. Output path. A `.tex` extension is appended
#'   when missing, and the parent directory is created if needed.
#' @param width,height \code{numeric(1)}. Device size in inches. Defaults `7`
#'   and `5`.
#' @param standalone \code{logical(1)}. Passed to [tikzDevice::tikz()] as
#'   `standAlone`. `TRUE` (default) writes a self-contained document that
#'   `pdflatex` can compile directly.
#' @param packages \code{character} or `NULL`. Extra LaTeX packages forwarded to
#'   [tikzDevice::tikz()] via its `package` argument, e.g.
#'   `c(getOption("tikzLatexPackages"), "\\usepackage{amsmath}")`.
#' @param compile \code{logical(1)}. If `TRUE`, run `pdflatex` on the generated
#'   file (requires `standalone = TRUE` and `pdflatex` on the `PATH`). Default
#'   `FALSE`.
#' @param clean \code{logical(1)}. If `TRUE` (default) and `compile` ran, remove
#'   the `.aux` and `.log` by-products.
#' @param ... Further arguments passed to [tikzDevice::tikz()].
#'
#' @return The path to the written `.tex` file, invisibly.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("tikzDevice", quietly = TRUE)) {
#'   data("data_far")
#'   g <- ggplot2::autoplot(estimate_mean(data_far, t = seq(0.1, 0.9, 0.1)))
#'   save_plot_tikz(g, file = file.path(tempdir(), "mean.tex"),
#'                  width = 7, height = 5)
#' }
#' }
#'
#' @seealso [tikzDevice::tikz()], [adaptiveFTS-autoplot].
#' @export
save_plot_tikz <- function(plot, file, width = 7, height = 5,
                           standalone = TRUE, packages = NULL,
                           compile = FALSE, clean = TRUE, ...) {
  if (!requireNamespace("tikzDevice", quietly = TRUE))
    stop("Package 'tikzDevice' is required by save_plot_tikz(); ",
         "please install it.", call. = FALSE)
  if (!grepl("\\.tex$", file)) file <- paste0(file, ".tex")
  dir_out <- dirname(file)
  if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

  tikzDevice::tikz(file = file, width = width, height = height,
                   standAlone = standalone, package = packages, ...)
  # Close the device even if printing errors.
  tryCatch(print(plot), finally = grDevices::dev.off())

  if (compile) {
    if (!standalone)
      stop("'compile = TRUE' requires 'standalone = TRUE'.", call. = FALSE)
    if (nzchar(Sys.which("pdflatex"))) {
      status <- system2(
        "pdflatex",
        args = c("-interaction=nonstopmode",
                 paste0("-output-directory=", shQuote(dir_out)),
                 shQuote(file)),
        stdout = FALSE, stderr = FALSE)
      if (!identical(as.integer(status), 0L))
        warning("pdflatex returned a non-zero status (", status, ").")
      if (clean) {
        base <- tools::file_path_sans_ext(basename(file))
        junk <- file.path(dir_out, paste0(base, c(".aux", ".log")))
        file.remove(junk[file.exists(junk)])
      }
    } else {
      warning("'pdflatex' not found on the PATH; skipping PDF compilation.")
    }
  }
  invisible(file)
}
