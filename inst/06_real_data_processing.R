library(data.table)
library(magrittr)
library(ggplot2)
library(latex2exp)

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
rm(dt_raw) ; gc()

# Convert voltage curve as nunmeric
dt[, Voltage := as.numeric(Voltage)]
date_na <- dt[is.na(Voltage), unique(date)]
dt <- dt[! date %in% date_na]
dt[, length(unique(date))]

# Plot curves ----
figures_path <- "./inst/12_mc_simulate_data/graphs/paper_graphs/"
theme_set(theme_minimal())

## All curves
g_real_data_all_curves <- ggplot(dt, aes(x = t, y = Voltage, color = Season, group = date)) +
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

# Save and clean
ggsave(plot = g_real_data_all_curves, filename = file.path(figures_path, "real_data_all_curves.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
rm(g_real_data_all_curves) ; gc()


# Mean function estimation ----
# Import data
dt <- dt[, .(date, "tobs" = t, "voltage" = Voltage)]

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
paste0(mu_coef, collapse = ", ")

mu <- predict(mu_model, mat_covariable)
dygraphs::dygraph(data = data.table::data.table(tobs, mu))

## mean function construction

get_real_data_mean <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:5, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    240.851203112216, 0.509378236915314, 0.0666785737279956, 0.402943145860831,
    0.161933581079031, 0.112863126651063, 0.420525704902966, 1.00346816098248,
    -0.242895339672357, -0.259141006436404, 0.00114630490804474
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
g_empirical_mean <- ggplot(dt_mu, aes(x = tobs, y = mu)) +
  geom_line() +
  ylim(238, 244) +
  xlab(label = "t") +
  ylab(label = expression(mu(t))) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
ggsave(plot = g_empirical_mean, filename = file.path(figures_path, "empirical_mean.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

g_smooth_mean <- ggplot(dt_smooth, aes(x = t, y = mean)) +
  geom_line() +
  ylim(238, 244) +
  ylab(label = expression(mu(t))) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
ggsave(plot = g_smooth_mean, filename = file.path(figures_path, "smooth_mean.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

rm(g_smooth_mean, g_empirical_mean) ; gc()

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

rm(mat_XsXt, mat_x, dt_x) ; gc()

## Plot covariance function
ggrid <- expand.grid(s = tobs, t = tobs)
dt_cov <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "cov" = c(mat_cov)
)

g_empirical_cov <- ggplot(dt_cov, aes(x = s, y = t, z = cov)) +
  geom_contour_filled(bins = 6) +
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
ggsave(plot = g_empirical_cov, filename = file.path(figures_path, "empirical_cov.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
rm(g_empirical_cov) ; gc()

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

rm(mat_musmut, mat_xn, mat_xn_plus_1, dt_xn_plus_1, dt_xn, mat_Xns_Xn_plus_1_t, d, mat_mu) ; gc()

## plot C_1(s,t)
ggrid <- expand.grid(s = tobs, t = tobs)
dt_autocov <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "autocov" = c(mat_autocov)
)

g_empirical_lag1_autocov <- ggplot(dt_autocov, aes(x = s, y = t, z = autocov)) +
  geom_contour_filled(bins = 6) +
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
ggsave(plot = g_empirical_lag1_autocov, filename = file.path(figures_path, "empirical_lag1_autocov.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
rm(g_empirical_lag1_autocov) ; gc()

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
paste0(beta_coef, collapse = ", ")

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
    0.887265486496153, -0.158284777828367, -0.433123270896265, -0.383368407909871,
    0.145655492369033, -0.00932791858596785, 0.25405721049976, 0.0360507006945973,
    0.0389539855934984, 0, -0.0133553863644848, 0.0177582032888235, 0.189761421268642,
    0.0195864450427664, 0.0887495150023169, 0, 0.0347257788913602, 0, 0.298938773778208,
    0.360062724244617, 0.00694075505838772, 0.0383993219719295, 0.0889742879270508,
    0.108124616829882, 0.597015339786177
  )
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
  op_norm <- 0.9128311
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef)
  gc()

  return(ker_values)
}

### Estimate the operator norm
svec <- (1:1500) / 1500
tvec <- (1:1500) / 1500

dt_grid <- expand.grid(s = svec, t = tvec)

dt_far_ker <- get_real_data_far_kenel(s = dt_grid$s, t = dt_grid$t, operator_norm = 0.7)
mat_far_ker <- matrix(dt_far_ker, ncol = length(tvec))

plotly::plot_ly(x = svec, y = tvec, z = mat_far_ker) %>%
  plotly::add_surface() %>%
  plotly::layout(title = "Operator kernel")

## Calcul de la norm
op_norm <- max(apply(X = mat_far_ker, MARGIN = 1, FUN = function(r){
  pracma::trapz(x = (1:1500) / 1500, y = r)
}))
op_norm

### Plot of the kernel function
ggrid <- expand.grid(s = (1:1440) / 1440, t = (1:1440) / 1440)
dt_kernel <- data.table::data.table(
  "s" = ggrid$s,
  "t" = ggrid$t,
  "Kernel_value" = get_real_data_far_kenel(s = ggrid$s, t = ggrid$t, operator_norm = 0.9128311)
)

#### Contour plot
g_far_kernel <- ggplot(dt_kernel, aes(x = s, y = t, z = Kernel_value)) +
  geom_contour_filled(bins = 6) +
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
ggsave(plot = g_far_kernel, filename = file.path(figures_path, "far_kernel.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

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

# Plot hurst function ----
# Local exponent function for the simulation
Hlogistic <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

dt_local_exponent <- data.table::data.table("t" = (1:1440) / 1440, "Ht" = Hlogistic((1:1440) / 1440))
g_local_exponent <- ggplot(data = dt_local_exponent, aes(x = t, y = Ht)) +
  geom_line(size = 0.9)  +
  ylim(0.35, 0.65) +
  xlab(label = "t") +
  labs(fill = latex2exp::TeX("$H_t$  ")) +
  scale_color_grey() +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = , r = 10, b = 0, l = 0)),
        axis.text.x =  element_text(size = 16),
        axis.text.y =  element_text(size = 16))
g_local_exponent
ggsave(plot = g_local_exponent, filename = file.path(figures_path, "local_exponent.png"),
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")

