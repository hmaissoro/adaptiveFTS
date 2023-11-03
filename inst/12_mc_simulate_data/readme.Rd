In this folder, we build a unique simulation setup to generate

- FAR(1) process using fBm and mfBm white noise and
- FMA(1) process using fBm and mfBm white noise,

to estimate

- The local regularity paramètres
- The mean function
- The lag-1 autocovariance function.

The data generation strategy is :

- to generate a zero-mean process
- and add to this process the mean function.

For the mean and autocovariance functions, we also estimate the MSE for each h in the bandwidth grid.

