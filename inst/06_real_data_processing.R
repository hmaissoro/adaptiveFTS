library(data.table)
library(magrittr)
library(ggplot2)

# Import data ----
dt_raw <- fread("../data_electricity/household_power_consumption.txt", sep = ";")

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
  theme(legend.position = "top")
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
  theme(legend.position = "top")
ggsave(filename = file.path(figures_path, "real_data_selected_curves.png"), units = "px", dpi = 300)
