# Load data
data("data_far")

# Estimate risk function
## Using One bandwith
dt_autocov_risk <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
  lag = 3, bw_grid = NULL, use_same_bw = TRUE,
  center = TRUE, kernel_name = "epanechnikov")

dt_autocov_risk[, .(s,t, hs, ht, autocov_risk)]

## Plot the risk function
dt_dcast <- data.table::dcast(data = dt_autocov_risk,
                              formula = hs ~ s + t ,
                              value.var = "autocov_risk")

manipulateWidget::combineWidgets(
  list = list(
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.2, 0.25)" = `0.2_0.25`)],
                      main = "lag = 3 - (s, t) = (0.2, 0.25)",
                      xlab = "h",
                      ylab = "risk function"),
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.4, 0.5)" = `0.4_0.5`)],
                      main = "lag = 3 - (s, t) = (0.4, 0.5)",
                      xlab = "h",
                      ylab = "risk function"),
    dygraphs::dygraph(data = dt_dcast[, .(hs, "(s, t) = (0.8, 0.75)" = `0.8_0.75`)],
                      main = "lag = 3 - (s, t) = (0.8, 0.75)",
                      xlab = "h",
                      ylab = "risk function")
  ),
  nrow = 3
)

## Using Two bandwiths
dt_autocov_risk_2bw <- estimate_autocov_risk(
  data = data_far, idcol = "id_curve", tcol = "tobs", ycol = "X",
  s = c(1/5, 2/5, 4/5), t = c(1/4, 1/2, 3/4),
  lag = 3, bw_grid = NULL, use_same_bw = FALSE,
  center = TRUE, kernel_name = "epanechnikov")

dt_autocov_risk[, .(s,t, hs, ht, autocov_risk)]
