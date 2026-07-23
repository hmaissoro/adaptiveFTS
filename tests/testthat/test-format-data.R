# The data contract every estimator relies on: three columns, curves renumbered
# in order of first appearance, sorted output, and rejection of input the
# estimators cannot honour.

# data.table functions are accessed via data.table::

long_fixture <- function(n_curves = 4L) {
  dt <- fixture_data_far(n_curves)
  format_data(dt, idcol = "id_curve", tcol = "tobs", ycol = "X")
}

test_that("the long format yields the documented columns, order and numbering", {
  out <- long_fixture()
  expect_true(data.table::is.data.table(out))
  expect_identical(names(out), c("id_curve", "tobs", "X"))
  expect_identical(sort(unique(out[["id_curve"]])), 1:4)
  expect_false(is.unsorted(out[["id_curve"]]))
  expect_true(all(out[, !is.unsorted(tobs), by = "id_curve"][["V1"]]))
})

test_that("the three input layouts give the same table", {
  from_long <- long_fixture()
  curves <- split(from_long, by = "id_curve", keep.by = FALSE)

  from_tables <- format_data(curves, tcol = "tobs", ycol = "X")
  expect_equal(from_tables, from_long)

  from_vectors <- format_data(
    lapply(curves, function(d) list(tobs = d[["tobs"]], X = d[["X"]])),
    tcol = "tobs", ycol = "X")
  expect_equal(from_vectors, from_long)
})

test_that("curves whose rows are not contiguous are still grouped correctly", {
  # Regression: the curve index used to be rebuilt from run lengths counted by
  # value but assigned in row order, which silently scattered the observation
  # points of a curve across its neighbours.
  reference <- long_fixture()
  raw <- fixture_data_far(4L)
  interleaved <- raw[c(seq(1L, nrow(raw), by = 2L), seq(2L, nrow(raw), by = 2L))]

  out <- format_data(interleaved, idcol = "id_curve", tcol = "tobs", ycol = "X")
  expect_equal(out, reference)
})

test_that("curves are numbered by first appearance, not by sorted index value", {
  raw <- fixture_data_far(12L)
  raw[["id_curve"]] <- paste0("curve_", raw[["id_curve"]])  # "curve_10" < "curve_2"
  out <- format_data(raw, idcol = "id_curve", tcol = "tobs", ycol = "X")
  expect_equal(out, long_fixture(12L))
})

test_that("the input is left untouched", {
  raw <- fixture_data_far(3L)
  before <- data.table::copy(raw)
  invisible(format_data(raw, idcol = "id_curve", tcol = "tobs", ycol = "X"))
  expect_equal(raw, before)
})

test_that("unusable input is rejected with a message naming the problem", {
  raw <- fixture_data_far(3L)
  curves <- split(long_fixture(3L), by = "id_curve", keep.by = FALSE)

  expect_error(format_data(42), "must be a data.table")
  expect_error(format_data(raw, tcol = "tobs", ycol = "X"), "'idcol' must name")
  expect_error(format_data(raw, idcol = "absent"), "no column named")
  expect_error(format_data(curves, idcol = "id_curve"), "must be NULL")
  expect_error(format_data(list(), tcol = "tobs", ycol = "X"), "empty list")
  expect_error(
    format_data(list(list(tobs = c(0.1, 0.2), X = 1)), tcol = "tobs", ycol = "X"),
    "different lengths")
  expect_error(
    format_data(data.frame(i = 1L, tobs = 1.5, X = 1), idcol = "i"),
    "must lie in \\[0, 1\\]")
  expect_error(
    format_data(data.frame(i = 1:2, tobs = c(0.1, NA), X = c(1, 2)), idcol = "i"),
    "missing values")
  expect_error(
    format_data(data.frame(i = 1:2, tobs = c("a", "b"), X = c(1, 2)), idcol = "i"),
    "must be numeric")
})

test_that("repeated observation points warn instead of passing silently", {
  expect_warning(
    format_data(data.frame(i = c(1L, 1L), tobs = c(0.5, 0.5), X = c(1, 2)),
                idcol = "i"),
    "repeated observation points")
})
