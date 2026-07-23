#' Estimate the functional autocorrelation function (FACF)
#'
#' Computes the functional autocorrelation \eqn{\widehat\rho_\ell} for lags
#' \eqn{\ell = 1, \dots, } `lag.max`, using the adaptive lag-\eqn{\ell}
#' autocovariance estimator [estimate_autocov()]:
#' \deqn{\widehat\rho_\ell = \lVert \widehat\Gamma_{N,\ell} \rVert
#'   \left[ \int \widehat\Gamma_{N,0}(t, t)\, dt \right]^{-1},}
#' where \eqn{\lVert\cdot\rVert} is the \eqn{\mathbb{L}^2} norm on
#' \eqn{[0,1]^2}. This is the functional analogue of the autocorrelation function
#' of a scalar time series and is useful for detecting serial dependence /
#' non-stationarity in a functional time series.
#'
#' @details
#' The lag-0 variance function \eqn{\widehat\Gamma_{N,0}(t,t)} (the denominator)
#' and each lag-\eqn{\ell} autocovariance surface (the numerator) are estimated
#' adaptively with [estimate_autocov()], which handles both the common and the
#' independent observation design. The integrals use the trapezoidal rule on the
#' evaluation grid `t`. Computing a full \eqn{|t| \times |t|} surface for each lag
#' can be costly; keep `n_grid` modest, pass a coarser `t`, or supply `bw_grid`.
#'
#' @inheritParams estimate_autocov
#' @param lag.max \code{integer(1)}. Largest lag. Default `5`. Clamped to
#'   `N - 1` (with a warning) when it reaches the number of curves.
#' @param t \code{numeric} or `NULL`. Evaluation grid in \eqn{[0,1]}. When `NULL`
#'   (default) it is chosen design-aware: the shared observation grid under the
#'   common design (subsampled to at most `n_grid` points), or a regular grid of
#'   `n_grid` points under the independent design.
#' @param n_grid \code{integer(1)}. Size of the default evaluation grid when `t`
#'   is `NULL`. Default `25`.
#'
#' @return An object of class `fts_acf` (a classed [data.table::data.table]; see
#'   [adaptiveFTS_est]) with columns:
#'   \itemize{
#'     \item `lag`: the lag \eqn{\ell}.
#'     \item `norm`: the \eqn{\mathbb{L}^2} norm \eqn{\lVert\widehat\Gamma_{N,\ell}\rVert}.
#'     \item `facf`: the functional autocorrelation \eqn{\widehat\rho_\ell}.
#'   }
#'
#' @references
#' Horváth, L., Rice, G. and Whipple, S. (2016). Adaptive bandwidth selection
#' in the long run covariance estimator of functional time series.
#' \emph{Computational Statistics and Data Analysis}, 100, 676--693.
#' \doi{10.1016/j.csda.2014.06.008}
#'
#' @seealso [estimate_autocov()], [autoplot.fts_acf()], [summary.fts_acf()].
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' data("data_far")
#' facf <- estimate_facf(data = data_far, lag.max = 5, n_grid = 20)
#' summary(facf)
#' if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::autoplot(facf)
#' }
estimate_facf <- function(data, idcol = "id_curve", tcol = "tobs", ycol = "X",
                          lag.max = 5L, t = NULL, n_grid = 25L,
                          bw_grid = NULL, common_bw = FALSE,
                          center_curves = TRUE, kernel_name = "epanechnikov") {
  data <- format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)
  kernel_name <- match.arg(
    arg = kernel_name,
    choices = c("epanechnikov", "biweight", "triweight", "tricube", "triangular", "uniform"))

  N <- data[, length(unique(id_curve))]
  lag.max <- as.integer(lag.max)
  if (length(lag.max) != 1L || is.na(lag.max) || lag.max < 1L)
    stop("'lag.max' must be a positive integer.")
  if (lag.max >= N) {
    lag.max <- N - 1L
    warning("'lag.max' >= the number of curves; using lag.max = ", lag.max, ".")
  }

  common <- .is_common_design(data = data, idcol = "id_curve", tcol = "tobs")
  if (is.null(t)) {
    if (common) {
      tt <- data[, sort(unique(tobs))]
      if (length(tt) > n_grid)
        tt <- tt[round(seq(1, length(tt), length.out = n_grid))]
      t <- tt
    } else {
      t <- seq(0.05, 0.95, length.out = n_grid)
    }
  }
  if (!(methods::is(t, "numeric") && all(data.table::between(t, 0, 1))))
    stop("'t' must be a numeric vector with values between 0 and 1.")
  t <- sort(unique(t))
  if (length(t) < 2L) stop("The evaluation grid 't' needs at least two points.")

  # Denominator: integral of the lag-0 variance function Gamma_0(t, t).
  dt_cov0 <- estimate_autocov(
    data = data, s = t, t = t, lag = 0L, bw_grid = bw_grid,
    common_bw = common_bw, center_curves = center_curves, correct_diagonal = FALSE,
    kernel_name = kernel_name)
  vec_cov0 <- dt_cov0[order(s), autocov]
  vec_cov0[is.nan(vec_cov0)] <- 0
  denom <- .trapz(t, vec_cov0)
  if (!is.finite(denom) || denom <= 0)
    stop("The lag-0 variance integral is non-positive; cannot normalise the FACF.")

  # Numerators: L^2 norm of each lag-l autocovariance surface over the t x t grid.
  grid <- data.table::CJ(s = t, t = t)
  norms <- vapply(seq_len(lag.max), function(l) {
    dt_l <- estimate_autocov(
      data = data, s = grid$s, t = grid$t, lag = l, bw_grid = bw_grid,
      common_bw = common_bw, center_curves = center_curves, correct_diagonal = FALSE,
      kernel_name = kernel_name)
    dt_l[is.nan(autocov), autocov := 0]
    m <- as.matrix(
      data.table::dcast(dt_l[order(s, t)], s ~ t, value.var = "autocov")[, -1L])
    inner <- apply(m ^ 2, 2L, function(z) .trapz(t, z))  # integrate over s
    sqrt(.trapz(t, inner))                                # then over t
  }, numeric(1))

  dt_facf <- data.table::data.table(
    lag = seq_len(lag.max), norm = norms, facf = norms / denom)
  .as_adaptive_est(
    dt_facf, "fts_acf",
    meta = list(kernel = kernel_name, N = N,
                design = if (common) "common" else "independent",
                n_grid = length(t), denom = denom))
}

#' @rdname adaptiveFTS-summary
#' @export
summary.fts_acf <- function(object, ...) {
  imax <- which.max(abs(object$facf))
  cat("Functional autocorrelation function (FACF)\n")
  cat(sprintf("  Design             : %s\n", .est_meta(object, "design")))
  cat(sprintf("  Evaluation grid    : %s points\n", .est_meta(object, "n_grid")))
  cat(sprintf("  Training curves    : %s\n", .est_meta(object, "N")))
  cat(sprintf("  Kernel             : %s\n", .est_meta(object, "kernel")))
  cat(sprintf("  Lags               : 1 to %d\n", max(object$lag)))
  cat(sprintf("  rho                : %s\n", .fmt_range(object$facf)))
  cat(sprintf("  Largest |rho|      : %s (lag %d)\n",
              .fmt_num(object$facf[imax]), object$lag[imax]))
  invisible(object)
}

#' @rdname adaptiveFTS-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.fts_acf <- function(object, ...) {
  .require_ggplot2()
  ggplot2::ggplot(object, ggplot2::aes(x = lag, y = facf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0),
                          linewidth = 0.6, colour = .col_main) +
    ggplot2::geom_point(size = 1.8, colour = .col_main) +
    ggplot2::scale_x_continuous(breaks = object$lag) +
    ggplot2::labs(title = "Functional autocorrelation (FACF)",
                  x = "lag", y = expression(hat(rho)[l])) +
    .gg_base()
}

#' @rdname adaptiveFTS-autoplot
#' @export
plot.fts_acf <- function(x, ...) {
  p <- autoplot.fts_acf(x, ...); print(p); invisible(p)
}
