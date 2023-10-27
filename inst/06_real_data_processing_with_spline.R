library(data.table)
library(magrittr)
library(ggplot2)
library(latex2exp)
library(splines2)

# Import data ----
dt_raw <- fread("../electricity_consumption_data/household_power_consumption.txt", sep = ";")

# Information about the data
names(dt_raw)
dt_raw[, length(unique(Date))]

# One observation per minute, expect "16/12/2006" and "26/11/2010" where some minutes are not observed
dt_raw[, .N, by = "Date"][, summary(N)]
dt_raw[! Date %in% c("16/12/2006", "26/11/2010"), .N, by = "Date"]

# Remove "16/12/2006" and "26/11/2010"
dt_raw <- dt_raw[! Date %in% c("16/12/2006", "26/11/2010")]
dt_raw[, date := as.Date(Date, format = "%d/%m/%Y")]

# Normalize time
dt_raw[, t := (1:1440) / 1440, by = "Date"]

# Define colour by season
yy <- dt_raw[, unique(year(date))]
dt_color <- data.table("date" = dt_raw[, unique(date)],
                       "season" = "", "color" = "")
for(yyi in  yy){
  date_winter <- c(seq(lubridate::as_date(paste0(yyi, "-12-22")), lubridate::as_date(paste0(yyi, "-12-31")), by = 1),
                   seq(lubridate::as_date(paste0(yyi, "-01-01")), lubridate::as_date(paste0(yyi, "-03-19")), by = 1))
  date_spring <- seq(lubridate::as_date(paste0(yyi, "-03-20")), lubridate::as_date(paste0(yyi, "-06-20")), by = 1)
  date_summer <- seq(lubridate::as_date(paste0(yyi, "-06-21")), lubridate::as_date(paste0(yyi, "-09-22")), by = 1)
  date_autumn <- seq(lubridate::as_date(paste0(yyi, "-09-23")), lubridate::as_date(paste0(yyi, "-12-21")), by = 1)

  dt_color[date %in% date_winter, c("Season", "color") := .("Winter", "#0095EF")]
  dt_color[date %in% date_spring, c("Season", "color") := .("Spring", "#6A38B3")]
  dt_color[date %in% date_summer, c("Season", "color") := .("Summer", "#FE433C")]
  dt_color[date %in% date_autumn, c("Season", "color") := .("Autumn", "#F39C12")]
}

# Add colour and extract only Voltage curves
dt <- data.table::merge.data.table(x = dt_raw[, .(t, date, Voltage)], y = dt_color, by = "date")
rm(dt_raw)
gc()

# Convert voltage curve as nunmeric
dt[, Voltage := as.numeric(Voltage)]
date_na <- dt[is.na(Voltage), unique(date)]
dt <- dt[! date %in% date_na]

# Plot curves ----
figures_path <- "../../../report/learning-smmoothness/Learning-smoothness/figures/"
theme_set(theme_minimal())

## All curves
ggplot(dt, aes(x = t, y = Voltage, color = Season, group = date)) +
  geom_line() +
  ylim(220, 255) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "real_data_all_curves.png"), units = "px", dpi = 300)

ggplot(dt, aes(x = t, y = Voltage, color = Season, group = date, linetype = Season)) +
  geom_line() +
  ylim(220, 260) +
  scale_linetype_manual(values = c("Autumn" = "solid", "Spring" = "dotted", "Summer" = "longdash", "Winter" = "twodash")) +
  scale_color_grey() +
  theme(legend.position = "top")

## Autumn 2009 to Winter 2010 : 167 curves
beginning_autumn <- lubridate::as_date("2009-09-23")
end_winter <- lubridate::as_date("2010-03-19")

dt_slice <- dt[between(date, beginning_autumn, end_winter)]
ggplot(dt_slice, aes(x = t, y = Voltage, color = Season, group = date)) +
  geom_line() +
  ylim(220, 255) +
  scale_color_grey() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width= unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "real_data_selected_curves.png"), units = "px", dpi = 300)


# Mean function estimation ----
# Import data
dt <- fread(file = "../electricity_consumption_data/household_voltage_autumn2009_winter2010.csv")

# Estimation of the Fourier basis coefficients
# Indeed \mu(t) is expressed in fourier basis
# Number of basis elements : 101

## Empirical mean function
dt_mu <- dt[order(tobs), .("mu" = mean(voltage, na.rm = TRUE)), by = "tobs"]
dygraphs::dygraph(dt_mu)

## Regression to estimate coefficients

### co-variables
tobs <- dt[, unique(sort(tobs))]
K <- 5

cos_mat <- outer(X = tobs, Y = 1:K, function(t, k) sqrt(2) * cos(2 * pi * k * t))
colnames(cos_mat) <- paste0("cos", 1:K)

sin_mat <- outer(X = tobs, Y = 1:K, function(t, k) sqrt(2) * sin(2 * pi * k * t))
colnames(sin_mat) <- paste0("sin", 1:K)

mat_covariable <- cbind(cos_mat, sin_mat)

### LASSO regression
cv_model <- glmnet::cv.glmnet(x = mat_covariable, y = dt_mu[, mu],
                              intercept = TRUE, alpha = 1, nfolds = 20)
plot(cv_model)
best_lambda <- cv_model$lambda.min

mu_model <- glmnet::glmnet(x = mat_covariable, y = dt_mu[, mu],
                           intercept = TRUE, alpha = 1, lambda = best_lambda)
mu_coef <- coef(mu_model)

mu <- predict(mu_model, mat_covariable)
dygraphs::dygraph(data = data.table::data.table(tobs, mu))

## mean function construction

get_real_data_mean <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    2.424666e+02, 2.045053e-01, 1.428404e-01, 3.570635e-01,
    1.090163e-01, 2.118234e-01, 1.301123e-01, 9.508016e-02,
    1.379302e-01, 4.820783e-03, -6.806055e-02, -2.388833e-02,
    -7.403338e-02, -2.313821e-02, -5.963765e-02, -2.788351e-02,
    6.017267e-02, 4.237483e-03, 2.406135e-02, 9.984997e-03,
    2.732683e-02, -3.947718e-02, -2.744963e-03, -5.624731e-03,
    -7.232138e-02, 1.228444e-02, -2.983196e-02, -1.373429e-02,
    -5.086314e-03, -5.206650e-03, 1.983353e-02, -1.671532e-02,
    1.694785e-02, 1.663588e-02, -8.924121e-03, -1.470650e-02,
    -1.568260e-02, -3.873785e-03, -1.642147e-02, 0.000000e+00,
    1.706651e-04, 3.662220e-03, 3.599005e-03, 2.163418e-02,
    2.079180e-02, -5.805163e-03, -6.198625e-03, -3.051126e-03,
    -3.050078e-02, -2.183053e-02, -2.030866e-02, 6.524725e-01,
    1.448886e+00, -3.113890e-01, -4.050726e-01, 1.478745e-01,
    -1.578363e-01, -2.346072e-01, -5.265881e-02, -5.181928e-02,
    -8.382788e-02, -5.055031e-02, 1.471149e-01, -8.402212e-03,
    -4.316013e-02, 7.528717e-02, 3.718024e-02, -4.602782e-03,
    -4.930040e-02, -7.104138e-03, -3.485272e-02, -5.034491e-02,
    2.230170e-02, 5.058664e-02, -2.996308e-02, 0.000000e+00,
    1.773518e-02, 1.664768e-03, -5.118570e-04, 2.536951e-02,
    1.103531e-02, -3.781447e-02, -9.837124e-03, 3.219296e-03,
    -1.163841e-02, -1.604513e-02, -8.183994e-03, 3.309498e-02,
    7.700235e-03, 1.578432e-02, 5.755486e-03, -3.571603e-03,
    -1.118589e-03, -1.883942e-03, -5.265843e-03, -2.892450e-02,
    1.032219e-02, 1.451413e-02, 6.348425e-04, 1.621365e-02,
    1.322576e-02
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

dt_smooth <- data.table("t" = tobs, "mean" = get_real_data_mean(t = tobs))

### Plot for paper
figures_path <- "../../../report/learning-smmoothness/Learning-smoothness/figures/"

ggplot(dt_mu, aes(x = tobs, y = mu)) +
  geom_line() +
  ylim(238.5, 246) +
  xlab(label = "t") +
  ylab(label = expression(mu(t))) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
ggsave(filename = file.path(figures_path, "empirical_mean.png"), units = "px", dpi = 300)

ggplot(dt_smooth, aes(x = t, y = mean)) +
  geom_line() +
  ylim(238.5, 246) +
  ylab(label = expression(mu(t))) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
ggsave(filename = file.path(figures_path, "smooth_mean.png"), units = "px", dpi = 300)

# Empirical covariance estimation ----
## E[X(s)X(t)]
dt_x <- data.table::dcast(data = dt, formula = tobs ~ date, value.var = "voltage")
dt_x <- dt_x[order(tobs)]
mat_x <- as.matrix(dt_x[, .SD, .SDcols = ! 'tobs'])
rownames(mat_x) <- dt_x[, tobs]
mat_XsXt <- (1 / ncol(mat_x)) * (mat_x %*% t(mat_x))

## mu(s)mu(t)
dt_mu <- dt_mu[order(tobs)]
mat_mu <- matrix(data = dt_mu[, mu], ncol = 1)
rownames(mat_mu) <- dt_mu[order(tobs), tobs]
mat_musmut <- mat_mu %*% t(mat_mu)

## C(s,t)
mat_cov <- mat_XsXt - mat_musmut

rm(mat_XsXt, mat_x, dt_x)

## Plot covariance function
ggrid <- expand.grid(s = tobs, t = tobs)
dt_cov <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "cov" = c(mat_cov)
)

ggplot(dt_cov, aes(x = s, y = t, z = cov)) +
  geom_contour_filled(breaks = c(-2, 0, 3, 6, 9, 11, 12)) +
  xlab(label = "s") +
  ylab(label = "t") +
  labs(fill = expression(C(s,t))) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
# , legend.key.width= unit(0.8, 'cm')) + guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(figures_path, "empirical_cov.png"), units = "px", dpi = 300)

# Empirical lag-1 auto-covariance ----

## E[X_n(s)X_{n+1}(t)]
d <- dt[, unique(date)]
dt_xn <- dt[! date == d[length(d)]]
dt_xn_plus_1 <- dt[! date == d[1]]
dt_xn <- dt_xn[order(tobs)]
dt_xn_plus_1 <- dt_xn_plus_1[order(tobs)]

dt_xn <- data.table::dcast(data = dt_xn, formula = tobs ~ date, value.var = "voltage")
dt_xn_plus_1 <- data.table::dcast(data = dt_xn_plus_1, formula = tobs ~ date, value.var = "voltage")
mat_xn <- as.matrix(dt_xn[, .SD, .SDcols = ! 'tobs'])
mat_xn_plus_1 <- as.matrix(dt_xn_plus_1[, .SD, .SDcols = ! 'tobs'])
rownames(mat_xn) <- dt_xn[, tobs]
rownames(mat_xn_plus_1) <- dt_xn_plus_1[, tobs]

mat_Xns_Xn_plus_1_t <- (1 / ncol(mat_xn)) * (mat_xn %*% t(mat_xn_plus_1))

## C_1(s,t)
mat_autocov <- mat_Xns_Xn_plus_1_t - mat_musmut

rm(mat_musmut, mat_xn, mat_xn_plus_1, dt_xn_plus_1, dt_xn, mat_Xns_Xn_plus_1_t, d, mat_mu)
gc()

## plot C_1(s,t)
ggrid <- expand.grid(s = tobs, t = tobs)
dt_autocov <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "autocov" = c(mat_autocov)
)

ggplot(dt_autocov, aes(x = s, y = t, z = autocov)) +
  geom_contour_filled(breaks = c(-15, -10, 0, 10, 20, 25, 30)) +
  xlab(label = "s") +
  ylab(label = "t") +
  labs(fill = latex2exp::TeX("$C_1(s,t)$")) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
ggsave(filename = file.path(figures_path, "empirical_lag1_autocov.png"), units = "px", dpi = 300)

# Operator kernel estimation ----
## LASSO regression
### fourier basis
tobs <- dt[, unique(sort(tobs))]
L <- 2

cos_mat <- outer(X = tobs, Y = 1:L, function(t, l) sqrt(2) * cos(2 * pi * l * t))
colnames(cos_mat) <- paste0("cos", 1:L)

sin_mat <- outer(X = tobs, Y = 1:L, function(t, l) sqrt(2) * sin(2 * pi * l * t))
colnames(sin_mat) <- paste0("sin", 1:L)

theta <- cbind(1, cos_mat, sin_mat)
eta <- cbind(1, cos_mat, sin_mat)
dim(eta)

### Compute Z_l(s) = \int_{0}^{1} c(s,u) theta(u) du
### Riemann approximation
ZZ <- (1 / length(tobs)) * mat_cov %*% theta

### Compute the co-variable matrix
### Note that the following is faster then : FF <- kronecker(ZZ, eta)
FF <- fastmatrix::kronecker.prod(eta, ZZ)

### Compute the dependent variable vector
### Note that c(matrix) = (matrix[, 1], matrix[, 2], ..., matrix[, ncol(matrix)])
### Thus by doing...
### CC1 <- c(mat_autocov)
### we obtain the vector
# vector(c_1(s_1, t_1),...,c_1(s_1440, t_1),
#        c_1(s_1, t_2),...,c_1(s_1440, t_2)
#        ,...,
#        c_1(s_1, t_1440),...,c_1(s_1440, t_1440))
#
# which is compatible with the above kronecker product
CC1 <- c(t(mat_autocov))

rm(mat_autocov, mat_cov, ZZ, theta, eta, cos_mat, sin_mat)

### Lasso regression
beta_cv_model <- glmnet::cv.glmnet(x = FF, y = CC1,
                                   intercept = FALSE, alpha = 1, nfolds = 10)
plot(beta_cv_model)
beta_lambda <- beta_cv_model$lambda.min

beta_model <- glmnet::glmnet(x = FF, y = CC1,
                             intercept = FALSE, alpha = 1, lambda = beta_lambda)
beta_coef <- coef(beta_model)

CC1_prev <- predict(beta_model, FF)

### Build kernel function
get_real_data_far_kenel <- function(s = 0.2, t = 0.3, operator_norm = 0.5){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    2.57941060039562, -4.65662130992959, 1.5225055060594, 4.3037464414024,
    0.13756488746035, -0.290043632023058, 0.307363812510174, -0.0185231460161384,
    -0.0294750471444012, 0.178548954594415, 0.149939543585357, 0.166076151730632,
    0.251806576768987, -0.0816563105543888, -0.0543167925659293, 0.03780271855891,
    -0.131443317726857, -0.125582864547152, 0.144243986892005, 0.134101454330207,
    0.508908339757622, -0.101104935110456, 0.30474744066116, 0.103156715437213,
    -0.0352343962613528)
  # Transform to (K, L) matrix
  basis_coef_mat <- t(matrix(data = basis_coef, ncol = 5))

  ker_values <- mapply(function(s,t, coef_mat){
    # \eta(s)
    etas <- c(1, sqrt(2) * cos(2 * pi * 1:2 * s), sqrt(2) * sin(2 * pi * 1:2 * s))

    # \theta(t)
    thetat <- c(1, sqrt(2) * cos(2 * pi * 1:2 * t), sqrt(2) * sin(2 * pi * 1:2 * t))

    # Basis function
    ker_val <- matrix(etas, nrow = 1) %*% coef_mat %*% matrix(thetat, ncol = 1)
    return(c(ker_val))
  }, s = s, t = t, MoreArgs = list(coef_mat = basis_coef_mat))

  # Normalize values using operator norm
  op_norm <- 4.588783
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef) ; gc()

  return(ker_values)
}

### Plot of the kernel function

ggrid <- expand.grid(s = (1:1440) / 1440, t = (1:1440) / 1440)
dt_kernel <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "Kernel_value" = get_real_data_far_kenel(s = ggrid$s, t = ggrid$t, operator_norm = 4.588783)
)

#### Contour plot
ggplot(dt_kernel, aes(x = s, y = t, z = Kernel_value)) +
  geom_contour_filled(breaks = c(-35, -25, -15, -5, 5, 15, 25)) +
  xlab(label = "s") +
  ylab(label = "t") +
  labs(fill = latex2exp::TeX("$\\psi(s,t)$  ")) +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
ggsave(filename = file.path(figures_path, "far_kernel.png"), units = "px", dpi = 300)

### Surface plot
svec <- seq(0.01, 0.99, len = 200)
tvec <- seq(0.01, 0.99, len = 200)
ggrid <- expand.grid(s = svec, t = tvec)
dt_kernel <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "Kernel_value" = get_real_data_far_kenel(s = ggrid$s, t = ggrid$t, operator_norm = 4.588783)
)
ker_mat <- matrix(dt_kernel$Kernel_value, ncol = length(tvec))
plot3D::persp3D(x = seq(0.01, 0.99, len = 200),
                y = seq(0.01, 0.99, len = 200),
                z = ker_mat,
                colvar = ker_mat,
                breaks = c(-35, -25, -15, -5, 5, 15, 25),
                col = RColorBrewer::brewer.pal(n = 8, "Greys")[8:3],
                xlab = "s", ylab = "t", zlab = "Θ(s,t)",
                ticktype = 'detailed', nticks = 5)

