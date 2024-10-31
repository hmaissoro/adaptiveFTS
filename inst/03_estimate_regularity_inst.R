# Load data
data("data_far")

# Estimate local regularity at
t0 <- seq(0.2, 0.8, len = 8)

## If data is a data.table or a data. frame
dt_locreg <- estimate_locreg(data = data_far,
                             idcol = "id_curve",
                             tcol = "tobs",
                             ycol = "X",
                             t = t0,
                             Delta = NULL,
                             h = NULL,
                             kernel_name = "epanechnikov",
                             center = TRUE)
DT::datatable(dt_locreg)
