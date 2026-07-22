# Select the kernel function used by the density estimator

Maps a kernel name to the corresponding (vectorised) kernel evaluated by
the package's compiled kernels. Only the kernels shared with the rest of
the package are supported.

## Usage

``` r
.select_density_kernel(kernel_name)
```

## Arguments

- kernel_name:

  A string giving the kernel name.

## Value

A function taking a numeric vector `u` and returning the kernel values
`K(u)`.
