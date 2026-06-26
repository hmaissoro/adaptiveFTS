# Access internal (non-exported) C++ wrappers and read reference fixtures.
g <- function(nm) getFromNamespace(nm, "adaptiveFTS")

read_ref <- function(nm) readRDS(test_path("_refs", paste0(nm, ".rds")))

# Standard reference inputs (must match scripts that generated tests/testthat/_refs).
ref_inputs <- function() {
  dt <- fixture_data_far(12L)
  list(
    dt   = dt,
    dtf  = as.data.frame(dt),
    N    = length(unique(dt$id_curve)),
    tt   = c(0.25, 0.5, 0.75),
    ss   = c(0.25, 0.5, 0.75),
    ttac = c(0.30, 0.50, 0.70),
    bwg  = seq(0.05, 0.20, length.out = 6),
    uvec = seq(-1.2, 1.2, by = 0.1),
    one  = dt[dt[["id_curve"]] == 1L][order(get("tobs"))],
    vg   = c(.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95)
  )
}
