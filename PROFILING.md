# Phase 2.0 — Profiling baseline

Establishes, with measurements, where runtime actually goes so the Phase-2
sub-phases are ordered by real impact (not assumption). All numbers are from the
post-Phase-1 numerical state.

## Method
- Platform: R 4.2.1 (ucrt) + rtools42, Windows, 8 cores.
- Fixture: deterministic 12-curve slice of `data_far` (`fixture_data_far(12)`),
  ~1080 rows, irregular per-curve designs.
- Timing: base `system.time` elapsed (single call; the dominant items are seconds
  so single-call timing is sufficient to rank). `bench`/`microbenchmark` are not
  installed; base R was used to avoid adding a dev dependency.
- Each exported `*_cpp` function and its R wrapper was exercised on the fixture
  and its output saved to `tests/testthat/_refs/*.rds`. These doubles as the
  Phase-3 regression anchors and the Phase-2 numerical-comparison references.

## Headline result
`predict_curve` / `estimate_curve_cpp` (the full adaptive BLUP pipeline) dominate
everything else by ~100×. On the 12-curve fixture:

| Function | elapsed |
|---|---|
| `predict_curve` (wrapper) | ~4.3 s |
| `estimate_curve_cpp` | ~4.2 s |
| `estimate_autocov_risk_cpp` (3 pairs) | 0.04 s |
| `estimate_autocov_cpp` (3 pairs) | 0.03 s |
| `estimate_mean_*`, `estimate_locreg_*` | 0.02 s |
| `estimate_cov_segment_*` | 0.01 s |
| kernels, matrix helpers, `estimate_nw*`, `estimate_sigma_cpp` | ~0 s |

## Where the pipeline time goes
`estimate_curve_cpp` (src/08_estimate_curve_cpp.cpp:555-558) builds a **fixed
19×19 = 361-point grid** and calls `estimate_autocov_risk_cpp` over **all 361
(s,t) pairs, twice** (lag 0 and lag 1), with the default bandwidth grid and
`use_same_bw = FALSE` (the O(bw²) two-bandwidth loop). The subsequent
`estimate_autocov_cpp` calls (lines 579-593) pass *precomputed* bandwidths, so
they skip the risk surface and are comparatively cheap.

Direct measurement of those two calls in isolation (same 361-grid, default bw,
2-bandwidth):

| Step | elapsed | share of `estimate_curve_cpp` |
|---|---|---|
| `estimate_autocov_risk_cpp` lag 0 (361 pairs) | ~5.6 s | |
| `estimate_autocov_risk_cpp` lag 1 (361 pairs) | ~5.4 s | |
| `estimate_mean_risk_cpp` (40 pts) | 0.04 s | |
| **two risk calls combined** | **~10.9 s** | **~78%** |

(Absolute seconds vary with machine load; the **~78% share** is the stable, and
decisive, figure.)

## Conclusion — sub-phase ordering
1. **2.1 `estimate_autocov_risk_cpp` (two-bandwidth risk) — TOP priority.** It is
   ~78% of the BLUP pipeline. Optimising it speeds up `estimate_curve_cpp` /
   `predict_curve` directly. Targets: cache per-bandwidth weight vectors and
   per-point local-regularity outside the O(bw²) loop; avoid arma copies; exploit
   symmetry; verify OpenMP has no hidden copies / reduction-order issues.
2. **2.6/diagonal correction in `estimate_autocov_cpp`** — the remaining ~22% of
   the pipeline (estimation + diagonal correction + BLUP solve). Second.
3. Everything else (`estimate_mean*`, `locreg`, `cov_segment`, `nw`, constants,
   helpers) is negligible on realistic inputs — optimise only opportunistically
   while in the relevant file, and only when numerically safe.

Note: the Rubin–Panaretos R wrappers (`estimate_mean_bw_rp`,
`estimate_autocov_bw_rp`, `estimate_autocov_rp`) are separately pathologically
slow (O(pairs) × per-call `format_data` + repeated `gc()`); they are not on the
adaptive BLUP path but are R-wrapper hotspots to address when their file is
touched.

## References / regression anchors
`tests/testthat/_refs/*.rds` — one per function (41 functions) plus
`_timings.rds`. Every Phase-2 change is validated with `all.equal()` against the
matching reference (fallback abs/rel tol 1e-10).

---

# Phase 6: porting the restructuring of the fast estimators (2026-11-06)

Branch `perf/optimal-core`. Reference: M. Vimond's `estimate_mean_fast()` /
`estimate_autocov_fast()` (FTS repository, `scripts/R/CHANGES.md` sections
16-18), which made the pure-R estimators 19-88x faster with identical results by
reorganising the computation: curves split once, per-curve quantities on the
distinct points, dependence terms from matrices, one smoothing pass, risk and
estimate in one pass, BLUP operators estimated directly at the points used.

## What remained in our C++ core
The same patterns, in C++ form: `arma::find()` scans of whole tables inside
nested loops (per bandwidth, point, lag or pair), `std::function` kernels
allocating a vector per call, sub-estimators re-converting the data and
re-running the presmoothing cross-validation, and per-scalar R closures (density
LSCV, Rubin-Panaretos `CJ()` products of all the observation pairs).

## What was done
- `src/00_core.h`: curves split once (`core::Curves`), scalar kernels,
  compact-support windows on the sorted points (the original 0/1 predicates are
  re-evaluated, so the indicators are unchanged).
- `src/00_internal.h`: `*_impl` entry points shared by all the estimators.
- Mean / covariance-segment risk: t x lag autocovariance matrix built once;
  `rho_l` from integer co-occurrence counts of the active curves.
- (Auto)covariance risk: point constants once on the distinct points of s and t;
  per-curve cubes from the windows; sums over the curves active on both sides
  (bit sets); fused argmin (the former finite-index argmin quirk is reproduced).
- (Auto)covariance estimate: one smoothing pass per (point, bandwidth) tuple for
  both the means and the cross products; diagonal correction over the
  intersection of the windows.
- Local regularity: one scalar quantile instead of a vector of identical ones.
- Deterministic OpenMP (each output computed by one thread, same summation
  order); no exception inside the parallel regions.
- R level: density LSCV in C++; Rubin-Panaretos double sums factorised (with an
  exact emulation of the former recycling when the mean vector is misaligned,
  see below); simulators split the curves once; `gc()` calls removed.
- BLUP: bandwidths selected directly at the points (no sub-grid). This is the
  only change of results, by design.

## Validation
A/B harness (scratch scripts): 63 cases, i.e. every exported function plus the
`_cpp` primitives, on the fixtures, an N = 150 / lambda = 40 random design, a
common design and edge cases (duplicated pairs, scalar and per-pair bandwidths,
tiny bandwidths, every kernel). Every non-BLUP output agrees with `main` within
1e-10 (largest relative difference 2.4e-13; the simulators and several
primitives are bit-identical). The test suite (402 expectations) passes.

## Timings (N = 150, lambda = 40, Windows, 8 threads)
| Call | main | branch | speed-up |
|---|---|---|---|
| `blup_fit()`, defaults | 33.2 s | 4.8 s | 6.9x |
| `select_tikhonov_parameter()` (bandwidths now re-selected in each fold) | 10.6 s | 4.6 s | 2.3x |
| `get_density_optimal_bw()` | 16.6 s | 0.17 s | 97x |
| `simulate_far()` | 13.4 s | 8.5 s | 1.6x |
| `estimate_facf()` | 16.6 s | 0.28 s | 59x |
| `estimate_mean_risk_cpp()`, 80 t | 3.6 s | 0.02 s | 178x |
| `estimate_mean()`, 80 t | 3.6 s | 0.03 s | 121x |
| `estimate_autocov_cpp()`, lag 0, bw NULL, 35 x 35 | 3.6 s | 0.04 s | 90x |
| `estimate_cov_segment()`, 80 t | 2.4 s | 0.06 s | 40x |
| `predict_curve()`, 40 curves | 1.1 s | 0.20 s | 5.7x |
| `estimate_autocov_rp()`, 20 curves, 4 pairs | 8.4 s | 0.01 s | ~840x |
| `estimate_mean_bw_rp()`, 20 curves, 5 bandwidths | 6.8 s | 0.30 s | 23x |
| `estimate_autocov_bw_rp()`, 4 curves x 4 points, 2 folds, 2 bandwidths | 440 s | 0.1 s | ~4000x |

`simulate_far()` is now dominated by `MASS::mvrnorm()` and the mfBm covariance,
which must stay as they are to keep the same random draws.

## Pre-existing behaviours kept (for a separate decision)
- `estimate_mean_cpp`, `estimate_cov_segment_cpp`, `estimate_autocov_cpp`: the
  optimal bandwidth is taken at the position of the minimum *among the finite
  risks*, applied to the list of *all* the candidates, so a non-finite risk
  before the minimum shifts the choice.
- `estimate_autocov_bw_rp()`: the mean estimated on the full sample is indexed
  with the renumbered ids of the training curves, and the cartesian products
  sorted the times, values and means independently; with unequal curve sizes
  the means are recycled or NA. Reproduced exactly (`rp_q_recycled_cpp`), which
  keeps that path quadratic in the number of observations.
- The two-bandwidth branch of the autocovariance risk squares `pn_t` in its
  variance term (noted earlier).

## Round 3: kernels that do less work (2026-11-06)

Measured at N = 300, lambda = 100 (one Tikhonov fold: 1.8 s before this round),
each idea validated on the A/B harness (all non-BLUP outputs within 1e-10 of
main; the BLUP identical to the previous round).

- **Nested supports.** pi_n(x; h) = 1{some |T - x| / h <= 1} only switches from
  0 to 1 as h grows (IEEE division is monotone), so along the sorted bandwidths
  each curve enters a support once. In the two-bandwidth risk, P_N, B_1, C_1
  (s-side quantities over the t-active curves) and B_2 become cumulative sums
  along the entry order: one pass over the curves per bandwidth instead of one
  dot product per bandwidth couple. Only C_3 = sum max W_s max W_t stays a
  (vectorised) dot product per couple: 5 dot products -> about 1. The nesting is
  checked per point; a dense fallback remains.
- **Dense vectorised dot products** (`omp simd`): inactive curves contribute exact
  zeros (pi_n = 0 implies b_n = max W = 0), so the sums run contiguously over all
  the curves instead of following bit sets.
- **Precomputed D terms**: D(u; h) = num / P_N^3 once per (point, bandwidth);
  P_N^3 of a count is an exact product; bandwidth powers once per pair.
- **Shared powers**: |(T - x)/h|^{2H} = |T - x|^{2H} h^{-2H}, so each observation's
  power is computed once for the 15 bandwidths (mean, covariance-segment and
  (auto)covariance risks).
- **Shared point constants** between C0 and C1 of the BLUP (local regularity,
  noise, moments, D numerators, per-curve cubes on T_d U t_pred).
- **Lag sums**: NaN-aware (values set to 0 plus a finiteness indicator, same
  mean-of-finite-products semantics), four lags per pass (independent
  accumulators, each lag summed in its own order).
- **Popcount lag counts** in the mean risk: c_l = popcount(pi AND (pi >> l)).
- Parallel D numerators and smoothing matrices.

| Call | main | round 2 | round 3 |
|---|---|---|---|
| `blup_fit()` defaults, N = 150, lambda = 40 | 33.2 s | 4.8 s | 2.6 s |
| `select_tikhonov_parameter()`, N = 150 | 10.6 s | 4.6 s | 2.3 s |
| `blup_fit()`, N = 300, lambda = 100 | (not run) | 58 s | 25.6 s |
| one Tikhonov fold, N = 300 | | 1.8 s | 0.8 s |
| `estimate_mean()`, 80 t, N = 150 | 3.6 s | 0.03 s | 0.01 s |
| `estimate_facf()`, N = 150 | 16.6 s | 0.28 s | 0.19 s |
| `estimate_cov_segment()`, 80 t | 2.4 s | 0.06 s | 0.03 s |

What remains in a fold at N = 300 is close to the work the estimator needs: the
X(s)X(t) autocovariances (all the lags of every pair, O(N^2) per pair) and the
C_3 products (B^2 N per pair).

## Round 4: thread control and the X(s)X(t) lag autocovariances as matrix products (2026-11-07)

- **`set_threads()` / `get_threads()`** (and `ADAPTIVEFTS_NUM_THREADS` at load):
  the OpenMP thread count of the compiled estimators. Results do not depend on
  it (tested: 1, 2 and the default give `identical()` outputs). For Monte Carlo
  studies, one worker per processor with `set_threads(1)` avoids
  oversubscription.
- **Matrix-product form of the X(s)X(t) dependence quantities** (risk of the
  (auto)covariance, dense grids of pairs). With P_k = a_k(s) b_{k+l}(t),
  sum (P_k - g)(P_{k+L} - g) = T1 - g (T2 + T3) + g^2 n_f, and
  T1 = sum [a_k a_{k+L}] [b_{k+l} b_{k+l+L}] is, for each lag, a product of an
  s-only and a t-only factor: one matrix product per lag for the whole grid
  (4 x 4 blocks of vectorised dot products). T2, T3 and n_f come from prefix sums
  of P, corrected on the few positions where a smoothed value is NaN (vectorised,
  all the lags at once). Pairs with very many NaN, duplicated pairs, sparse pair
  sets and other lag sets keep the direct computation.

One Tikhonov fold at N = 300, lambda = 100 with the CV-selected presmoothing
bandwidth: 0.87 s -> 0.7 s; `select_tikhonov_parameter()` there: 25.9 s -> 23.2 s.
At N = 150: one fold 0.10 s -> 0.08 s.
