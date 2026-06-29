# Sample for Functional Autoregressive Process of Order 1

This dataset contains a set of 150 curves, each with an average of 90
observation points and some added noise.

## Usage

``` r
data_far
```

## Format

An object of class `data.table` (inherits from `data.frame`) with 13346
rows and 3 columns.

## Details

A `data.table` with three columns:

- `id_curve` ::

  The index of each curve.

- `tobs` ::

  The observation points for each `id_curve`.

- `X` ::

  The observed values of the curve at each `tobs` point.
