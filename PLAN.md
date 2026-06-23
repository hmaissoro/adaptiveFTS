# adaptiveFTS — Cleanup, Performance, Tests, CRAN & CI Plan

## Context

`adaptiveFTS` implements the estimators from Maissoro, Patilea & Vimond,
"Adaptive prediction for functional time series" — adaptive BLUP for stationary
functional time series observed with error at irregular, random design points.
The statistical core is in `src/` (Rcpp + RcppArmadillo); the R layer in `R/`
orchestrates, validates inputs, and does cross-validation.

The package works but is not release-grade: it pulls in `caret` for one helper,
carries an undeclared `parallel` dependency, has a runtime bug that disables
`get_nw_optimal_bw()`, has **zero automated tests**, and **no CI**. The goal is a
CRAN-acceptable, well-tested, faster package whose numerical results are provably
unchanged (beyond floating-point tolerance) by every refactor.

**Top priority: numerical correctness.** No refactor may change estimator output
beyond fp tolerance. Every change is gated on the numerical-comparison rule below.

### Environment confirmed
- R 4.2.1 at `C:\Program Files\R\R-4.2.1` (not on PATH; invoke directly).
- rtools42 present (matches R 4.2.x) → C++ can be compiled & checked locally.
- LinkingTo RcppArmadillo; OpenMP enabled via `$(SHLIB_OPENMP_CXXFLAGS)`.

### Numerical-comparison rule (applies to EVERY refactor)
1. Before changing a function, capture its output on fixed **seeded** inputs and
   save as a reference `.rds` under `tests/testthat/_refs/` (or a scratch dir
   pre-Phase-3).
2. After the change, re-run on the SAME inputs and compare with `all.equal()`.
3. Only when `all.equal()` is impractical (structurally different objects), fall
   back to absolute/relative tolerance `1e-10`.
4. Accept the change ONLY if the comparison passes; otherwise revert and report.
5. Every PR that could touch statistical behavior shows the before/after diff.

### Workflow (all phases)
- One task = one branch off up-to-date `main` = one PR. Suggested branch names:
  `chore/baseline-fixes`, `feat/remove-caret`, `perf/autocov-rcpp`,
  `test/exported-coverage`, `cran/check-clean`, `ci/actions-matrix`.
- Commit per logical step; clear messages. **Never commit to `main`.**
- Never change statistical behavior without showing before/after in the PR.

---

## Baseline inventory (from exploration)

**Deps:** Imports = caret, data.table, fastmatrix, MASS, methods, Rcpp, Rdpack;
LinkingTo = Rcpp, RcppArmadillo; Suggests = knitr, rmarkdown. `parallel` is used
but **undeclared**. `stats` is used (`rpois`,`runif`) but **not in Imports**.
- `caret` → only `caret::createFolds` at `R/04_estimate_mean.R:430` and
  `R/05_estimate_autocovariance.R:652` (K-fold CV index creation).
- `fastmatrix` → only `kronecker.prod` in `estimate_mean_rp` (2 calls). Keep
  (lightweight, specialised) but evaluate a base-R `kronecker()` swap.
- `parallel` → `mclapply(..., mc.cores = 20)` at `R/05_...:655,712` (serial on
  Windows; hardcoded 20).

**Pre-existing bug:** `R/02_smoothing.R:243` calls
`get_nw_optimal_bw_cpp(..., kernel_name = )` with no value → `get_nw_optimal_bw()`
errors; this cascades into `R/07_estimate_constants.R` callers.

**src/ hotspots (ranked):** diagonal correction in `estimate_autocov_cpp`
(O(n_curve·n²·n_points) nested weight eval, `05_...:725-740`); two-bandwidth
`estimate_autocov_risk_cpp` (O(bw²·n·n_curve)); `estimate_numerator_dependence_term_DD_cpp`;
`estimate_locreg_cpp` (3 NW calls per curve per t). OpenMP already on bandwidth
loops in mean/autocov/cov-segment risk.

**Tests:** none. **CI:** none. **Examples:** `predict_curve` has none;
`estimate_cov_segment[_risk]` say "Example coming soon"; `format_data` none; many
others `\dontrun`. `data-raw/data_far.R` is a `usethis` stub, not real gen code.

---

## Phase 0 — Baseline & enabling fixes  (branch: `chore/baseline-fixes`)

Goal: establish a reproducible baseline and unblock broken paths before any
refactor, so later phases have something to compare against.

Tasks:
1. **Commit `PLAN.md`** to repo root (this document).
2. **Capture the CRAN baseline.** Run `R CMD build` then
   `R CMD check --as-cran` with R 4.2.1; record verbatim every ERROR/WARNING/NOTE
   into `PLAN.md` (a "Baseline check" appendix). This is the 0-of-everything
   target for Phase 4.
3. **Fix the `kernel_name = ` bug** at `R/02_smoothing.R:243` → pass
   `kernel_name = kernel_name`. This *enables* behavior rather than changing it
   (the function currently cannot run). Validate by calling `get_nw_optimal_bw()`
   on a seeded slice of `data_far` and confirming it returns a finite scalar;
   record the value as a regression reference.
4. **Declare hidden deps:** add `stats` and (temporarily) `parallel` to Imports
   so the baseline check is honest; `parallel` is removed again in Phase 1.
5. **Create a seeded fixture script** `tests/testthat/helper-fixtures.R` (used
   from Phase 3 on, but authored now) that builds a small deterministic dataset
   (e.g. `set.seed(42); simulate_far(...)` subset, or a fixed slice of
   `data_far`) reused by all reference captures.

Files: `R/02_smoothing.R`, `DESCRIPTION`, new `PLAN.md`, new fixture helper.
Validation: package builds; `get_nw_optimal_bw()` runs; baseline check recorded.
Risks: R 4.2.1 is older than CRAN's current release — some NOTES may be
version-specific; we'll re-check under newer R in CI (Phase 5).

---

## Phase 1 — Remove unnecessary dependencies  (branch: `feat/remove-caret`)

Goal: drop `caret` entirely and tidy the dependency list, with zero change to
seeded results.

Tasks:
1. **Re-implement `createFolds` faithfully (exact reproduction).** Add an
   internal `.create_folds(y, k, list = TRUE)` in a new `R/utils_folds.R` that
   replicates caret's stratification + assignment algorithm so that, **under the
   same `set.seed`, the returned fold lists are identical to caret's**. caret's
   algorithm for numeric `y`: cut `y` into ~`groups = min(5, ...)` percentile
   bins, then within each bin assign a shuffled, recycled `1:k` labelling via
   `sample()`. Replicate this exactly (note: here `y = unique(id_curve)`, integer
   1:N, so binning is on integers).
   - Validation: for several seeds and `(N, k)` combos,
     `identical(.create_folds(...), caret::createFolds(...))` must be TRUE. Then
     swap the two call sites (`R/04...:430`, `R/05...:652`) and confirm seeded
     `estimate_mean_bw_rp` / `estimate_autocov_bw_rp` outputs match the Phase-0
     reference via `all.equal()`.
2. **Replace `parallel::mclapply`** (`R/05...:655,712`) with plain `lapply`
   (deterministic, cross-platform). Optionally expose a `ncores`/backend arg
   defaulting to serial; if kept, gate behind a declared backend. Remove
   `parallel` from Imports (added in Phase 0).
   - Validation: seeded `estimate_autocov_bw_rp` output `all.equal()` to the
     pre-change reference (serial vs the old Windows-serial mclapply are already
     identical; on Linux confirm too).
3. **Evaluate `fastmatrix`.** Check whether `fastmatrix::kronecker.prod` can be
   replaced by base `kronecker()`/`%x%` with identical results
   (`all.equal`). If yes, drop `fastmatrix`; if perf-sensitive, keep and justify.
4. **Remove `caret` from DESCRIPTION Imports**; run `R CMD check` to confirm no
   residual references (also grep for `caret`, `parallel`).

Files: new `R/utils_folds.R`, `R/04_estimate_mean.R`, `R/05_estimate_autocovariance.R`,
`R/04...rp` (kronecker), `DESCRIPTION`, `NAMESPACE` (if import tags change).
Validation: dependency-removal grep clean; all affected seeded estimators
`all.equal()` to Phase-0 references; `R CMD check` no new notes.
Risks: caret's `createFolds` internals must be matched precisely — the validation
gate (identical fold lists) catches any drift; if exact match proves infeasible
for some edge config, surface it rather than silently diverging.

---

## Phase 2 — Performance: ALL of src/ + the R wrappers (split into sub-phases)

Goal: cut runtime with no change to results, across **every** `src/` module
(not just autocov) **and** the R wrapper functions in `R/` that call them.
Profile first, then work module-by-module; each sub-phase is its own branch/PR,
each function independently validated by the numerical rule and reverted if it
fails. Report before/after `bench::mark` medians per change.

Coverage map — every compiled unit and its R wrapper(s) is in scope:

| src file | C++ functions | R wrapper(s) |
|---|---|---|
| `02_smoothing_rcpp.cpp` | 6 kernels, `select_kernel`, `estimate_nw_cpp`, `estimate_nw_bw_cpp`, `get_nw_optimal_bw_cpp` | `R/02_smoothing.R`: `estimate_nw`, `estimate_nw_bw`, `get_nw_optimal_bw`, kernel fns |
| `03_estimate_locreg_rcpp.cpp` | `estimate_locreg_cpp` | `R/03_estimate_regularity.R`: `estimate_locreg` |
| `04_estimate_mean_cpp.cpp` | `estimate_mean_risk_cpp`, `estimate_mean_cpp` | `R/04_estimate_mean.R`: `estimate_mean`, `estimate_mean_risk`, `estimate_mean_rp`, `estimate_mean_bw_rp` |
| `05_estimate_autocov_cpp.cpp` | `estimate_autocov_risk_cpp`, `estimate_autocov_cpp`, `get_upper_tri_couple`, `sort_by_columns` | `R/05_estimate_autocovariance.R`: `estimate_autocov`, `estimate_autocov_risk`, `estimate_autocov_rp`, `estimate_autocov_bw_rp` |
| `06_estimate_cov_segment_cpp.cpp` | `estimate_cov_segment_risk_cpp`, `estimate_cov_segment_cpp` | `R/06_estimate_covariance_segment.R`: `estimate_cov_segment`, `estimate_cov_segment_risk` |
| `07_estimate_constants_cpp.cpp` | `estimate_sigma_cpp`, `estimate_empirical_mom_cpp`, `estimate_empirical_autocov_cpp`, `estimate_empirical_XsXt_autocov_cpp`, `estimate_numerator_DD_single_t`, `estimate_numerator_dependence_term_DD_cpp` | `R/07_estimate_constants.R`: `estimate_sigma`, `estimate_empirical_mom`, `estimate_empirical_autocov`, `estimate_empirical_XsXt_autocov` |
| `08_estimate_curve_cpp.cpp` | `estimate_curve_cpp`, `reshape_matrix`, `combine_matrices`, `build_grid`, `get_best/nearest_*_bw`, `ensure_positive_definite` | `R/08_predict_curve.R`: `predict_curve` |

### Phase 2.0 — Profiling baseline  (branch: `perf/profile-baseline`)
Use `profvis`/`bench` on representative seeded workflows that exercise each
module (full `predict_curve`/`estimate_curve`, plus each estimator with no
precomputed bandwidth). Record flame graphs and top-N per workflow. This produces
the empirical hotspot ranking that orders the sub-phases below; **do not edit C++
before this confirms where time actually goes.** Capture seeded reference outputs
for every exported `*_cpp` function and every R wrapper here (shared with Phase 3
refs). Files: scratch profiling scripts only. Validation: profile artifacts saved.

### Sub-phases (ordered by 2.0's measured ranking; representative scope each)
Each does, for its module: hoist invariants out of loops; pass arma objects by
`const&`; `reserve()/set_size()` before fill loops; replace R-ish element loops
with bit-identical arma vectorised ops; remove redundant recomputation; exploit
operator symmetry/structure; review OpenMP for false sharing and reduction-order
safety; and tighten the **R wrapper** (avoid `data.table` copies, redundant
`format_data`/`gc()` calls, `expand.grid` blowups, repeated re-smoothing).

- **2.1 `perf/autocov-rcpp`** — `05_estimate_autocov_cpp.cpp` (expected #1):
  diagonal correction `src/05...:725-740` (precompute per-curve kernel weights,
  upper-triangle-only via existing `get_upper_tri_couple`, kill copies);
  two-bandwidth `estimate_autocov_risk_cpp` `src/05...:289-414` (cache per-bw
  weights & per-point regularity outside the O(bw²) loop). + `R/05` wrappers
  (`estimate_autocov_bw_rp` `expand.grid`/CV path).
- **2.2 `perf/constants-rcpp`** — `07_estimate_constants_cpp.cpp`, esp.
  `estimate_numerator_dependence_term_DD_cpp` / `..._single_t` (DD term feeds
  multiple risk fns) and `estimate_empirical_XsXt_autocov_cpp`. + `R/07` wrappers.
- **2.3 `perf/locreg-rcpp`** — `03_estimate_locreg_rcpp.cpp`: reuse NW smooths
  across the 3 evaluation offsets where mathematically identical; cut per-t
  quantile re-filtering. + `R/03` wrapper.
- **2.4 `perf/mean-cov-segment-rcpp`** — `04_estimate_mean_cpp.cpp` &
  `06_estimate_cov_segment_cpp.cpp` (shared risk-term structure; OpenMP critical
  sections). + `R/04`, `R/06` wrappers.
- **2.5 `perf/smoothing-rcpp`** — `02_smoothing_rcpp.cpp`: `estimate_nw_cpp`
  weight-matrix construction, `estimate_nw_bw_cpp` CV loop, `get_nw_optimal_bw_cpp`
  subset loop; vectorise kernels. + `R/02` wrappers.
- **2.6 `perf/curve-pipeline-rcpp`** — `08_estimate_curve_cpp.cpp`: orchestration
  & helpers (`ensure_positive_definite` eig cost, `reshape_matrix`,
  `get_nearest_best_*_bw`). + `R/08` `predict_curve` wrapper. Push any remaining
  hot R-level loops into Rcpp **only where 2.0 justifies it** (no speculative
  ports).

Files per sub-phase: the named `src/*.cpp` + its `.h`, the matching `R/*.R`
wrapper, and `src/RcppExports.cpp` / `R/RcppExports.R` (regenerated only if a
signature changes — avoid signature changes).
Validation (every sub-phase, every function): seeded reference from 2.0; after
change `all.equal()` (fallback abs/rel tol `1e-10`); validate OpenMP paths with
threads=1 **and** threads>1; report before/after `bench::mark` medians.
Risks: floating-point reassociation from reordering sums or enabling threads can
exceed tolerance — keep reduction order, document any order-sensitive sums, and
revert any change that fails the gate.

---

## Phase 3 — Test suite  (branch: `test/exported-coverage`)

Goal: a `testthat` (edition 3) suite where **every exported function has ≥1
test**, plus regression tests that lock current seeded numerics.

Tasks:
1. Scaffold `tests/testthat.R` + `tests/testthat/`; add `testthat` to Suggests;
   `Config/testthat/edition: 3`.
2. **Coverage tests** for all 38 exports (kernels, `format_data`, simulators,
   Hurst fns, `estimate_locreg/mean/autocov/cov_segment` + risk/rp/bw variants,
   `estimate_sigma`, empirical-moment/autocov fns, NW smoothing, `predict_curve`,
   real-data fns, exported `.Qpq_fun`/`.Spq_fun`).
3. **Regression tests** using the Phase-0/seeded references: snapshot estimator
   outputs (`expect_equal` with explicit tolerance) so later phases can't silently
   change results. Store refs under `tests/testthat/_refs/`.
4. **Edge cases:** sparse curves (small `M_n`), varying per-curve `M_n`, boundary
   domain points (t at 0/1), single-curve / minimal-N inputs, lag-0 vs lag>0.
5. **Rcpp-path tests:** call the `*_cpp` functions directly (they're exported to
   R) to cover C++ branches (`use_same_bw` T/F, `correct_diagonal` T/F,
   `center` T/F, each kernel).
6. Keep each test fast (small data) so example/test runtime stays CRAN-legal.

Files: new `tests/` tree, `DESCRIPTION` (Suggests + Config).
Validation: `devtools::test()` green; `covr::package_coverage()` shows every
export hit at least once.
Risks: long-running estimators — use minimal fixtures and precomputed bandwidths
to keep tests sub-second where possible.

---

## Phase 4 — CRAN-acceptable submission  (branch: `cran/check-clean`)

Goal: `R CMD check --as-cran` → 0 errors / 0 warnings / 0 notes.

Tasks (driven by the Phase-0 baseline list):
1. **Docs/examples:** add runnable `@examples` for `predict_curve`, `format_data`,
   `estimate_cov_segment[_risk]` (replace "coming soon"); make examples fast or
   `\donttest` (not `\dontrun`) so they execute under check; ensure each export
   has full roxygen (`@param`,`@return`,`@examples`). Regenerate `man/`.
2. **DESCRIPTION:** real ORCID or drop the `YOUR-ORCID-ID` placeholder; valid
   `Authors@R` (add Patilea, Vimond if appropriate); tighten `Description`
   prose (no starting article, single paragraph, methods + `<doi:...>`/`<arXiv:...>`
   refs); confirm `License: AGPL (>= 3)` + `LICENSE.md` consistency.
3. **Rcpp registration & Makevars:** confirm `useDynLib(.registration=TRUE)`,
   regenerated `RcppExports`, `CXX_STD`/OpenMP flags portable; no compiler
   warnings under `-Wall`.
4. **Hygiene:** no `options()`/`par()` mutation without restore (none found —
   keep it that way), no writing outside `tempdir()`, no `<<-` to globalenv;
   `.Rbuildignore` covers `inst/*_inst*.R` dev scripts if they trip examples.
5. **Vignettes:** ensure they build without heavy/Suggested-but-undeclared pkgs
   (dygraphs, plotly, manipulateWidget, ggplot2 appear in vignettes) — either add
   to Suggests with `eval=FALSE`/guarded chunks, or trim. Spelling
   (`devtools::spell_check()`), URL/DOI validation (`urlchecker`).
6. **`data-raw/data_far.R`:** replace the `usethis` stub with the real, seeded
   generation script that reproduces `data/data_far.rda`.

Files: all `R/*.R` roxygen blocks, `man/*`, `DESCRIPTION`, `vignettes/*.Rmd`,
`src/Makevars*`, `data-raw/data_far.R`, `.Rbuildignore`.
Validation: `R CMD check --as-cran` clean under R 4.2.1 locally; also note any
items that only newer R flags (defer their fix to CI confirmation).
Risks: example runtime limits — gate expensive examples with `\donttest` and tiny
inputs; vignette deps are the most likely NOTE source.

---

## Phase 5 — CI/CD  (branch: `ci/actions-matrix`)

Goal: GitHub Actions with a full 3-OS check matrix + quality gates.

Tasks:
1. **`R-CMD-check.yaml`** (r-lib/actions/check-r-package): matrix over
   **ubuntu-latest, macos-latest, AND windows-latest**, each compiling the C++,
   with R **release + oldrel** (and optionally devel on Ubuntu). Use
   `r-lib/actions/setup-r` + `setup-r-dependencies`.
2. **`test-coverage.yaml`** (covr → Codecov), add codecov badge.
3. **`lint.yaml`** (lintr; add `.lintr` config tuned to the code style).
4. **`pkgdown.yaml`** build + deploy to `gh-pages`; add `_pkgdown.yml`.
5. README badges (check, coverage, CRAN-status placeholder).

Files: `.github/workflows/{R-CMD-check,test-coverage,lint,pkgdown}.yaml`,
`.lintr`, `_pkgdown.yml`, `codecov.yml`, `README.md`.
Validation: all matrix jobs green on the PR; coverage uploads; pkgdown builds.
Risks: Windows + OpenMP + RcppArmadillo build flakiness — pin Rtools, ensure
`Makevars.win` flags are CI-correct; macOS OpenMP may need a guard.

---

## Cross-phase awareness (do NOT act on now)
A later Python port is anticipated. Where cheap, keep the numerical core in `src/`
cleanly separated and well-specified (clear function contracts, no R-only
coupling in the math) so porting is straightforward — but no Python work and no
scope changes for it now.

## Decisions locked
- caret folds → **reproduce caret's algorithm exactly** (seeded fold lists must be
  identical); validate before swapping call sites.
- `parallel::mclapply` → **replace with portable serial `lapply`** (optional
  opt-in backend), drop undeclared `parallel`.
- `R/02_smoothing.R:243` bug → **fix in Phase 0** as a behavior-enabling prereq.

## Overall verification strategy
- Every estimator change gated by `all.equal()` (fallback `1e-10`) vs a seeded
  pre-change reference; failing changes are reverted and reported.
- End state: `R CMD check --as-cran` 0/0/0 on 3 OSes × {release, oldrel}; full
  testthat suite green; every export covered; CI/coverage/lint/pkgdown passing.

---

## Appendix A — Baseline `R CMD check --as-cran` (to be filled in Phase 0)

_Pending: output of `R CMD check --as-cran` under R 4.2.1 will be recorded here._
