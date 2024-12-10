
# **adaptiveFTS**

[![Version](https://img.shields.io/badge/version-0.1.1-blue)](https://github.com/hmaissoro/adaptiveFTS)  
[![License](https://img.shields.io/badge/license-AGPL%20%3E%3D%203-lightgrey)](LICENSE)

---

## **Description**

`adaptiveFTS` provides tools for adaptive estimation procedures for weakly dependent functional time series. Developed by Maissoro et al. (2024, 2025), it includes:

- Estimators for local regularity parameters
- Mean function estimation
- Autocovariance function estimation
- An adaptive Best Linear Unbiased Predictor (BLUP)

This package leverages the computational efficiency of Rcpp and RcppArmadillo for high-performance functional data analysis.

---

## **Table of Contents**

1. [Features](#features)  
2. [Installation](#installation)  
3. [Usage](#usage)  
4. [Examples](#examples)  
5. [Contributing](#contributing)  
6. [License](#license)  

---

## **Features**

- Adaptive estimators for functional time series analysis.
- Efficient computation using Rcpp and RcppArmadillo.
- Implements methodologies from cutting-edge research (2024, 2025).
- Tools for predicting functional time series via BLUP.

---

## **Installation**

### **From GitHub**
To install the development version of **adaptiveFTS**, use the `devtools` package:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install adaptiveFTS from GitHub
devtools::install_github("hmaissoro/adaptiveFTS")
```

---

## **Usage**

Load the package and start using its functions:

```R
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

(Refer to the [examples section](#examples) for detailed code snippets and use cases.)

---

## **Examples**

Example of mean function estimation

```R
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

---

## **Contributing**

We welcome contributions! Follow these steps:

1. Fork the repository: [adaptiveFTS](https://github.com/hmaissoro/adaptiveFTS)
2. Clone your fork:
   ```bash
   git clone https://github.com/your-username/adaptiveFTS.git
   ```
3. Create a feature branch:
   ```bash
   git checkout -b feature-name
   ```
4. Commit your changes:
   ```bash
   git commit -m "Add feature description"
   ```
5. Push the branch:
   ```bash
   git push origin feature-name
   ```
6. Submit a pull request.

---

## **License**

This package is licensed under the **AGPL (>= 3)** license. See the [LICENSE](LICENSE) file for more details.

---

## **Acknowledgements**

- Built using **Rcpp** and **RcppArmadillo** for efficient computation.
- Inspired by methodologies developed in Maissoro et al. (2024, 2025).

---

## **References**

**Hassan Maissoro, Valentin Patilea, and Myriam Vimond.**  
*Adaptive Estimation for Weakly Dependent Functional Time Series.*  
2024. Available at [arXiv:2403.13706](https://arxiv.org/abs/2403.13706).

**Hassan Maissoro, Valentin Patilea, and Myriam Vimond.**  
*Adaptive prediction for Functional Times Series.*  
2024. [Work in progress](https://hassan.maissoro.com/assets/pdf/2024-adaptive-estimation-for-functional-time-series.pdf).



---
