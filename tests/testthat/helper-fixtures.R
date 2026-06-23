# Deterministic fixtures shared by the test suite (Phase 3) and by the
# reference-capture scripts used to validate refactors (Phases 1-2).
#
# Every fixture is fully reproducible: it is either a fixed slice of the
# packaged `data_far` dataset (no RNG) or a `set.seed()`-seeded simulation.
# Keep these small so estimator calls stay fast and CRAN-legal.

# A small, deterministic slice of the packaged FAR dataset.
# No RNG involved, so the result is identical on every platform/session.
fixture_data_far <- function(n_curves = 15L) {
  e <- new.env()
  utils::data("data_far", package = "adaptiveFTS", envir = e)
  dt <- data.table::as.data.table(e$data_far)
  dt[dt[["id_curve"]] <= n_curves, ]
}

# A small, seeded FAR simulation for edge cases (sparse / irregular designs).
# Reproducible given the same `seed` because simulate_far() draws from R's RNG.
fixture_far_small <- function(seed = 42L, N = 20L, lambda = 30L) {
  set.seed(seed)
  simulate_far(
    N = N, lambda = lambda,
    tdesign = "random",
    Mdistribution = stats::rpois,
    tdistribution = stats::runif,
    tcommon = NULL,
    hurst_fun = hurst_logistic,
    L = 4,
    far_kernel = function(s, t) 9 / 4 * exp(-(t + 2 * s) ** 2),
    far_mean = function(t) 4 * sin(1.5 * pi * t),
    int_grid = 100L,
    burnin = 100L,
    remove_burnin = TRUE
  )
}

# A common evaluation grid used across reference captures.
fixture_tgrid <- function() c(0.2, 0.4, 0.5, 0.6, 0.8)
