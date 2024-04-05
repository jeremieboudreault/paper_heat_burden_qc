# main.R


# Main script for the computation of the heat burden for multiple health outcomes.


# Project : paper_heat_burden_qc
# Author  : Jeremie Boudreault
# Email   : Jeremie.Boudreault@inrs.ca
# Depends : R [v4.3.0], rjutils [v0.1]
# Imports : See below.
# License : CC BY-NC-SA 4.0 DEED


# Packages ---------------------------------------------------------------------


library(data.table)
library(dlnm)
library(splines)
library(mixmeta)
library(rjutils) # To install: remotes::install_github("jeremieboudreault/rjutils")


# Functions --------------------------------------------------------------------


source("R/attrdl.R")
source("R/create_lagged_mat.R")


# Imports ----------------------------------------------------------------------


# Daily health data for each HO by RSS.
health_data <- fread("data/health_data_synthetic.csv")

# Daily lagged weather and air pollution daily data by RSS.
weather_data <- fread("data/weather_data.csv")

# Meta-predictors for each HO and RSS.
meta_data <- fread("data/meta_predictors.csv")

# Thresholds for extreme heat by RSS.
heat_thresh <- fread("data/heat_thresholds.csv")


# Parameters -------------------------------------------------------------------


# Maximum lag [8 values from day 0 (current) to day 7 lag]
max_lag <- 7 

# Number of degrees of freedom for date variable
df_date <- 4

# Knots for x_lag to put into arglag_tvar.
nk_xlag <- 2

# Position of knots for x main exposure in quantiles.
knots_xvar <- c(50, 90)/100

# Quantiles for the MMT search.
q_mmt <- seq(25L, 98L, by = 1L)/100

# Months for the analysis.
months <- 5:9

# Health outcome for the analysis.
hvar <- c("MOR", "HOS", "EDV", "AMB", "811")[1L]

# Temperature exposure variable.
tvar <- c("TMEAN", "TMAX", "TMIN", "HMDX", "TDEW")[1L]

# Control for air pollution.
polctrl  <- c("polno", "polo3", "polo3pm25")[1L]

# Number of simulations for AN/AF calculation.
nsim <- 1000L
cols_sim <- paste0("SIM", 1:nsim)

# Verbose.
verbose <- FALSE


# ------------------------------------------------------------------------------
# Step 0 : Data preparation ----------------------------------------------------
# ------------------------------------------------------------------------------


# Merge data based on selected "hvar".
merged_data <- merge(
    x     = health_data[HO == hvar, ], 
    y     = weather_data, 
    by    = c("RSS", "DATE"), 
    all.x = TRUE
)[order(RSS, DATE), ]

# Extract all RSS to be treated.
rss_sel <- unique(merged_data$RSS)

# Convert weekday and year variables to factor.
merged_data[, WEEKDAY_F := factor(WEEKDAY, levels = 1:7)]
merged_data[, YEAR_F := factor(YEAR)]

# Fix <DOY> for leap years.
merged_data[YEAR %% 4 == 0L, DOY := DOY - 1]

# Create lagged matrix of tvar and relh.
tvar_mat_l <- create_lagged_mat(tvar)
relh_mat_l <- create_lagged_mat("RELHMEAN")

# Create lagged matrix of o3 and pm25 (if needed).
if (grepl("o3", polctrl)) o3_mat_l <- create_lagged_mat("O3MEAN")
if (grepl("pm25", polctrl)) pm25_mat_l <- create_lagged_mat("PM25MEAN")

# Create final dataset using selected months only.
data_final <- merged_data[MONTH %in% months, ]
set(data_final, j = tvar, value = data_final[[paste0(tvar, "0")]])

# Select meta-predictor based on current HO.
meta_data <- meta_data[HO == hvar, -c("HO")]

# Check that matrices fits with results data.
sum(sapply(tvar_mat_l, nrow)) == nrow(data_final)
sum(sapply(relh_mat_l, nrow)) == nrow(data_final)


