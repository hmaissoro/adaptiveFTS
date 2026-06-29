# **adaptiveFTS**

[![Version](https://img.shields.io/badge/version-0.1.1-blue)](https://github.com/hmaissoro/adaptiveFTS)  
[![License](https://img.shields.io/badge/license-AGPL%20%3E%3D%203-lightgrey)](https://hmaissoro.github.io/adaptiveFTS/LICENSE)
[![R-CMD-check](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/test-coverage.yaml)
[![lint](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/lint.yaml/badge.svg)](https://github.com/hmaissoro/adaptiveFTS/actions/workflows/lint.yaml)
[![codecov](https://codecov.io/gh/hmaissoro/adaptiveFTS/branch/main/graph/badge.svg)](https://app.codecov.io/gh/hmaissoro/adaptiveFTS)

------------------------------------------------------------------------

## **Description**

`adaptiveFTS` provides tools for adaptive estimation procedures for
weakly dependent functional time series. Developed by Maissoro et
al. (2024, 2025), it includes:

- Estimators for local regularity parameters
- Mean function estimation
- Autocovariance function estimation
- An adaptive Best Linear Unbiased Predictor (BLUP)

This package leverages the computational efficiency of Rcpp and
RcppArmadillo for high-performance functional data analysis.

------------------------------------------------------------------------

## **Table of Contents**

1.  [Features](#features)  
2.  [Installation](#installation)  
3.  [Usage](#usage)  
4.  [Examples](#examples)  
5.  [Development](#development)  
6.  [Contributing](#contributing)  
7.  [License](#license)

------------------------------------------------------------------------

## **Features**

- Adaptive estimators for functional time series analysis.
- Efficient computation using Rcpp and RcppArmadillo.
- Implements methodologies from cutting-edge research (2024, 2025).
- Tools for predicting functional time series via BLUP.

------------------------------------------------------------------------

## **Installation**

### **From CRAN**

Once released on CRAN, install the stable version with:

``` r

install.packages("adaptiveFTS")
```

### **From GitHub**

To install the development version of **adaptiveFTS**, use the
`devtools` package:

``` r

# Install devtools if not already installed
install.packages("devtools")

# Install adaptiveFTS from GitHub
devtools::install_github("hmaissoro/adaptiveFTS")
```

------------------------------------------------------------------------

## **Usage**

Load the package and start using its functions:

``` r

# Load adaptiveFTS
library(adaptiveFTS)

# Example of mean function estimation

## Load data
data("data_far")

## Estimate risk function for the mean
dt_mean_risk <- estimate_mean_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL,
  kernel_name = "epanechnikov"
)
```

(Refer to the [examples section](#examples) for detailed code snippets
and use cases.)

------------------------------------------------------------------------

## **Examples**

Example of mean function estimation

``` r

# Load data
data("data_far")

# Estimate risk function for the mean
dt_mean_risk <- estimate_mean_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = NULL,
  kernel_name = "epanechnikov"
)

# Visualize mean risk at various observation points
dt_dcast <- data.table::dcast(data = dt_mean_risk, formula = h ~ t, value.var = "mean_risk")
manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.25" = `0.25`)],
      main = "t = 0.25", xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.5" = `0.5`)],
      main = "t = 0.5", xlab = "h", ylab = "Risk Function"),
    dygraphs::dygraph(
      data = dt_dcast[, list(h, "t = 0.75" = `0.75`)],
      main = "t = 0.75", xlab = "h", ylab = "Risk Function")
  ),
  nrow = 3
)

# Estimate mean function with optimal bandwidths
dt_mean <- estimate_mean(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  t = c(1/4, 1/2, 3/4), bw_grid = seq(0.005, 0.15, len = 45),
  kernel_name = "epanechnikov"
)

# Display rounded estimates of the mean function
DT::datatable(data = dt_mean[, lapply(.SD, function(X) round(X, 3))])
```

------------------------------------------------------------------------

## **Development**

A `Makefile` wraps the common development and CRAN-submission tasks. Run
them from the package root (R, `Rscript` and a C++ toolchain must be on
the `PATH`):

``` bash
make help        # list all targets
make document    # regenerate Rd files and NAMESPACE (roxygen2)
make test        # run the testthat suite
make check       # build + R CMD check --as-cran on the tarball
make coverage    # test-coverage report (covr)
make pkgdown     # build the documentation website
make cran        # pre-submission gate: document + build + check --as-cran
```

If R is not on the `PATH` (e.g. on Windows), override the interpreters:

``` bash
make check R="/c/Program Files/R/R-4.2.1/bin/x64/R.exe" \
           RSCRIPT="/c/Program Files/R/R-4.2.1/bin/x64/Rscript.exe"
```

Continuous integration (GitHub Actions) runs `R CMD check` on Linux,
macOS and Windows (R-release and R-oldrel, plus R-devel on Linux), test
coverage, lint, and the pkgdown site on every push and pull request.

------------------------------------------------------------------------

## **Contributing**

We welcome contributions! Follow these steps:

1.  Fork the repository:
    [adaptiveFTS](https://github.com/hmaissoro/adaptiveFTS)

2.  Clone your fork:

    ``` bash
    git clone https://github.com/your-username/adaptiveFTS.git
    ```

3.  Create a feature branch:

    ``` bash
    git checkout -b feature-name
    ```

4.  Commit your changes:

    ``` bash
    git commit -m "Add feature description"
    ```

5.  Push the branch:

    ``` bash
    git push origin feature-name
    ```

6.  Submit a pull request.

------------------------------------------------------------------------

## **License**

This package is licensed under the **AGPL (\>= 3)** license. See the
[LICENSE](https://hmaissoro.github.io/adaptiveFTS/LICENSE) file for more
details.

------------------------------------------------------------------------

## **Acknowledgements**

- Built using **Rcpp** and **RcppArmadillo** for efficient computation.
- Inspired by methodologies developed in Maissoro et al. (2024, 2025).

------------------------------------------------------------------------

## **References**

1.  **Hassan Maissoro, Valentin Patilea, and Myriam Vimond.** *Adaptive
    Estimation for Weakly Dependent Functional Time Series.* 2024.
    Available at [arXiv:2403.13706](https://arxiv.org/abs/2403.13706).

2.  **Hassan Maissoro, Valentin Patilea, and Myriam Vimond.** *Adaptive
    prediction for Functional Times Series.* 2024. [Work in
    progress](https://hassan.maissoro.com/assets/pdf/2024-adaptive-estimation-for-functional-time-series.pdf).

------------------------------------------------------------------------
