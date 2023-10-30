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

# Local regularity estimation ----
# Import data
dt <- fread(file = "../electricity_consumption_data/household_voltage_autumn2009_winter2010.csv")
## Global parameters
N <- length(dt[, unique(date)])
lambda <- mean(dt[, .N, by = "date"][, N])
date_vec <- dt[, unique(date)]
delta <- exp(-log(lambda) ** 1/3)

## Estimate bandwidth
K <- 100
b0 <- 1 / lambda
bK <- lambda ** (- 1 / 3)
a <- exp((log(bK) - log(b0)) / K)
bw_grid <- b0 * a ** (seq_len(K))
rm(b0, bK, a, K) ; gc()

bw <- unlist(lapply(date_vec, function(di, data, bw_grid){
  estimate_nw_bw(
    y = data[date == di, voltage],
    t = data[date == di, tobs],
    bw_grid = bw_grid, smooth_ker = epanechnikov)
}, data = dt, bw_grid))

## Estimate the local regularity parameters
dt_locreg <- estimate_locreg(
  data = dt, idcol = "date", tcol = "tobs", ycol = "voltage",
  t = seq(0.1, 0.9, len = 100), Delta = delta,
  h = bw, smooth_ker = epanechnikov)

dygraphs::dygraph(data = dt_locreg[, .(t, Ht)])


# Mean function estimation ----
# Import data
dt <- fread(file = "../electricity_consumption_data/household_voltage_autumn2009_winter2010.csv")

# Estimation of the Fourier basis coefficients
# Indeed \mu(t) is expressed in N cubic spline basis
# Number of basis elements : 10

## Empirical mean function
dt_mu <- dt[order(tobs), .("mu" = mean(voltage, na.rm = TRUE)), by = "tobs"]
dygraphs::dygraph(dt_mu)

## Regression to estimate coefficients

### co-variables
tobs <- dt[, unique(sort(tobs))]
K <- 10

## Define basis
mat_covariable <- splines2::nsp(x = tobs, df = (K + 1 + 1), intercept = TRUE)

### LASSO regression
cv_model <- glmnet::cv.glmnet(x = mat_covariable, y = dt_mu[, mu],
                              intercept = TRUE, alpha = 1, nfolds = 20)
plot(cv_model)
best_lambda <- cv_model$lambda.min

mu_model <- glmnet::glmnet(x = mat_covariable, y = dt_mu[, mu],
                           intercept = FALSE, alpha = 1, lambda = best_lambda)
mu_coef <- coef(mu_model)

mu <- predict(mu_model, mat_covariable)
dygraphs::dygraph(data = data.table::data.table(tobs, mu))

## mean function construction

get_real_data_mean <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  eta <- splines2::nsp(x = t, df = 12, intercept = TRUE)

  # Basis coeffients
  basis_coef <- matrix(
    data = c(731.389571901671, -0.52006982827689, 246.098166747946,
             243.059070562643, 239.632261351878, 242.246665268361,
             241.476817297158, 246.281122481508, 241.410779068207,
             238.415397856336, -5.9613322334532, 734.400332002359),
    ncol = 1)

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(eta, basis_coef) ; gc()

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
L <- 5

theta <- splines2::nsp(x = tobs, df = (L + 1 + 1), intercept = TRUE)
eta <- splines2::nsp(x = tobs, df = (L + 1 + 1), intercept = TRUE)
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

rm(mat_autocov, mat_cov, ZZ, theta, eta)

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
get_real_data_far_kenel <- function(s = 0.2, t = 0.3, operator_norm = 0.7){
  # Basis coefficient
  # For each fixed {\eta_k(s), k = 1,...,K} and {\theta_l(t), l = 1,...,L}, we have
  # c(b_{11}, b_{12}, ..., b_{1L},
  #   b_{21}, b_{22}, ..., b_{2L},
  #   ...,
  #   b_{K1}, b_{K2}, ..., b_{KL})
  basis_coef <- c(
    -152.339167113767, 109.891801325509, 8.01218483618224, 40.7773902997773,
    9.29538508447265, -130.869348227434, 114.349788382373, -86.6452258918943,
    66.6620550406495, 2.85516854777547, 5.4726154735248, 8.16996388562305,
    -56.2974864050636, 82.3309005212671, -104.821148135886, 70.806466322378,
    3.20308498330693, 6.34208585976186, 12.8240898172477, -68.4244432976903,
    57.5035523478083, -79.009963744684, 55.5058015942419, 7.42288493639332,
    20.1537864363113, 2.35376692860042, -61.7830117509254, 77.2298006322408,
    -80.679220971239, 59.4135944681524, 4.1358055271817, 18.1150446807183,
    5.8204767292337, -60.5952551453644, 60.9830028643384, -137.595821369872,
    83.2753892033171, -7.64204301659673, -4.06678434607622, 11.8936539033712,
    -66.1910800695624, 85.8780261168699, -142.247945431835, 121.819129805756,
    9.45741883242566, 39.8753096159282, 6.85395396362438, -132.236594114549,
    128.212492255768)

  # Transform to (K, L) matrix
  basis_coef_mat <- matrix(data = basis_coef, ncol = 1)

  # \eta(s)
  etas <- splines2::nsp(x = s, df = 5 + 1 + 1, intercept = TRUE)

  # \theta(t)
  thetat <- splines2::nsp(x = t, df = 5 + 1 + 1, intercept = TRUE)

  mat_basis <- lapply(1:ncol(etas), function(id_eta_col, eta_mat, theta_mat){
    eta_mat[, id_eta_col] * theta_mat
  }, eta_mat = etas, theta_mat = thetat)
  mat_basis <- do.call(cbind, mat_basis)

  # Basis function
  ker_values <-  mat_basis %*% basis_coef_mat
  ker_values <-  c(ker_values)

  # Normalize values using operator norm
  op_norm <- 3.446037
  op_scale <- operator_norm / op_norm
  ker_values <- ker_values * op_scale

  # clean
  rm(basis_coef_mat, basis_coef, etas, thetat, mat_basis) ; gc()

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

# Estimate the norm ----
svec <- (1:1500) / 1500
tvec <- (1:1500) / 1500

dt_grid <- expand.grid(
  s = svec, t = tvec
)

dt_far_ker <- get_real_data_far_kenel(
  s = dt_grid$s, t = dt_grid$t, operator_norm = 0.7
)
mat_far_ker <- matrix(dt_far_ker, ncol = length(tvec))

plot_ly(x = svec, y = tvec, z = mat_far_ker) %>%
  add_surface() %>%
  layout(title = "Operator kernel")

## Calcul de la norm
op_norm <- max(1 / length(tvec) * mat_far_ker %*% matrix(data = 1, nrow = length(tvec)))
op_norm

op_norm <- max(apply(X = mat_far_ker, MARGIN = 1, FUN = function(r){
  pracma::trapz(x = (1:1500) / 1500, y = r)
}))
op_norm


