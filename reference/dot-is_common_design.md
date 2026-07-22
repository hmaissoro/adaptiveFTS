# Check whether all curves share the same observation design

Check whether all curves share the same observation design

## Usage

``` r
.is_common_design(data, idcol = "id_curve", tcol = "tobs")
```

## Arguments

- data:

  A prepared functional data.table.

- idcol, tcol:

  Identifier and observation-time column names.

## Value

`TRUE` if every curve shares the same sorted design, else `FALSE`.
