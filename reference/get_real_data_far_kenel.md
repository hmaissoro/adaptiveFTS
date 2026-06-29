# FAR kernel learned from the voltage curves of the electricity

For more details see the vignette:
`vignette("hybrid-simulation-setup", package = "adaptiveFTS")`

## Usage

``` r
get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
```

## Arguments

- s:

  `numeric (positive)`. A vector or scalar value(s) between 0 and 1.

- t:

  `numeric (positive)`. A vector or scalar value(s) between 0 and 1.

- operator_norm:

  `numeric (positive)`. A scalar corresponding to the norm of the
  integral operator associated with this kernel function.

## Value

A vector (or scalar) of `numeric` values corresponding to the value of
the kernel function evaluated at (`s`, `t`).

## Examples

``` r

# get the value of the kernel at (s,t) = (0.2, 0.3)
kerval <- get_real_data_far_kenel(s = 0.2, t = 0.3, operator_norm = 0.5)
kerval
#> [1] 0.4533972
```
