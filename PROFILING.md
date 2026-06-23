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
