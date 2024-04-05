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
source("R/fit_meta_stepwise_aic.R")


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
verbose <- TRUE


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

# Create lagged matrix of o3 and pm25 (only used if needed).
o3_mat_l <- create_lagged_mat("O3MEAN")
pm25_mat_l <- create_lagged_mat("PM25MEAN")

# Get RSS with air pollution data.
rss_pol_o3 <- rss_sel[sapply(o3_mat_l, function(x) sum(is.na(x)) != length(x))]
rss_pol_pm25 <- rss_sel[sapply(pm25_mat_l, function(x) sum(is.na(x)) != length(x))]

# Create final dataset using selected months only.
data_final <- merged_data[MONTH %in% months, ]
set(data_final, j = tvar, value = data_final[[paste0(tvar, "0")]])

# Select meta-predictor based on current HO.
meta_data <- meta_data[HO == hvar, -c("HO")]

# Check that matrices fits with results data.
sum(sapply(tvar_mat_l, nrow)) == nrow(data_final)
sum(sapply(relh_mat_l, nrow)) == nrow(data_final)


# ------------------------------------------------------------------------------
# Step 1.1 : DLNM model fitting by region (RSS) --------------------------------
# ------------------------------------------------------------------------------


# Model specification for temperature variable and control (relh, o3 and pm).
arglag_tvar <- list(fun = "ns", knots = logknots(max_lag, nk = nk_xlag))
argvar_ctrl <- list(fun = "lin")
arglag_ctrl <- arglag_tvar  

# Empty elements to store results.
coef_list  <- as.list(rep(NA, length.out = length(rss_sel)))
vcov_list  <- as.list(rep(NA, length.out = length(rss_sel)))
cr_list    <- as.list(rep(NA, length.out = length(rss_sel)))
simul_list <- as.list(rep(NA, length.out = length(rss_sel)))

# Fit DLNM by region (RSS).
for (rss in rss_sel){ 

    # Message.
    if (verbose) message("Fitting DLNM on RSS ", rss, ".")
    
    # Filter data based on <rss>.
    data_rss <- data_final[RSS == rss, ]
    nyears <- length(unique(data_rss$YEAR))
    
    # Extract indicator <rss_i>.
    rss_i <- which(rss_sel == rss)
    
    # Extract temperature and control variables.
    tvar_mat <- tvar_mat_l[[rss_i]]
    relh_mat <- relh_mat_l[[rss_i]]
    if (grepl("o3", polctrl) & rss %in% rss_pol_o3) o3_mat <- o3_mat_l[[rss_i]]
    if (grepl("pm25", polctrl) & rss %in% rss_pol_pm25) pm25_mat <- pm25_mat_l[[rss_i]]
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_rss[[tvar]], probs = knots_xvar), 
        Bound = range(data_rss[[tvar]])  
    ) 
    
    # Create cross-basis.
    cb_tvar <- crossbasis(tvar_mat, lag = max_lag, argvar = argvar_tvar, arglag = arglag_tvar)
    cb_relh <- crossbasis(relh_mat, lag = max_lag, argvar = argvar_ctrl, arglag = arglag_ctrl) 
    if (grepl("o3", polctrl) & rss %in% rss_pol_o3)
        cb_o3 <- crossbasis(o3_mat, lag = max_lag, argvar = argvar_ctrl, arglag = arglag_ctrl) 
    if (grepl("pm25", polctrl) & rss %in% rss_pol_pm25) 
        cb_pm25 <- crossbasis(pm25_mat, lag = max_lag, argvar = argvar_ctrl, arglag = arglag_ctrl) 
    
    # Create formula.
    form <- paste0("COUNT ~ cb_tvar + cb_relh + WEEKDAY_F + HOL + YEAR_F * ns(DOY, df = df_date)")
    if (grepl("o3", polctrl) & rss %in% rss_pol_o3) form <- paste0(form, " + cb_o3")
    if (grepl("pm25", polctrl) & rss %in% rss_pol_pm25) form <- paste0(form, " + cb_pm25")
    
    # Fitting DLNM with a quasiPoisson dsitribution.
    dlnm_rss <- glm(
        formula = as.formula(form), 
        family  = quasipoisson(), 
        data    = data_rss
    ) 
    
    # Reduce DLNM effect using mean values as MMT.
    cb_tvar_reduce <- crossreduce(
        basis = cb_tvar, 
        model = dlnm_rss, 
        cen   = mean(data_rss[[tvar]]),
        at    = quantile(data_rss[[tvar]], probs = q_mmt)
    )
    
    # Extract MMT and extreme temperature threshold.
    mmt <- quantile(data_rss[[tvar]], probs = q_mmt[which.min(cb_tvar_reduce$RRfit)])
    ext <- max(mmt, heat_thresh[RSS == rss, get(tvar)])
    
    # Set seed.
    set.seed(2912L)
    
    # Simulate AN for all heat.
    an_heat_all <- attrdl(
        x     = tvar_mat,
        basis = cb_tvar,
        cases = data_rss$COUNT,
        model = dlnm_rss,
        type  = "an",
        dir   = "back",
        range = c(mmt, max(tvar_mat[, 1L])),
        cen   = mmt,
        sim   = TRUE,
        nsim  = nsim
    )/nyears
    
    # Simulate AN for extreme heat.
    an_heat_ext <- attrdl(
        x     = tvar_mat,
        basis = cb_tvar,
        cases = data_rss$COUNT,
        model = dlnm_rss,
        type  = "an",
        dir   = "back",
        range = c(ext, max(tvar_mat[, 1L])),
        cen   = mmt,
        sim   = TRUE,
        nsim  = nsim
    )/nyears
    
    # Rearrange simulated AN into a data.table.
    simul <- data.table(
        RSS    = rss,
        TRANGE = c("heat_all", "heat_extreme")
    )
    simul <- cbind(simul, rbind(an_heat_all, an_heat_ext))
    colnames(simul)[3:(3+nsim-1)] <- cols_sim
    
    # Reduced function at MMT for plotting.
    cb_tvar_reduce <- crossreduce(
        basis = cb_tvar, 
        model = dlnm_rss, 
        cen   = mmt,
        at    = seq(min(data_rss[[tvar]]), max(data_rss[[tvar]]), by = 0.1)
    )
    
    # Save information for future steps.
    cr_list[[rss_i]]    <- cb_tvar_reduce
    coef_list[[rss_i]]  <- coef(cb_tvar_reduce) 
    vcov_list[[rss_i]]  <- vcov(cb_tvar_reduce)
    simul_list[[rss_i]] <- copy(simul)
    
}

# Note : 'cr_list' countains the reduced function for each RSS (i.e., regional DLNM).
# For example : 
# plot(cr_list[[1L]], xlab = tvar, y = "RR")

# Convert saved coefficient into a matrix.
coef_matrix <- matrix(
    data     = unlist(coef_list),
    nrow     = length(coef_list), 
    ncol     = length(coef_list[[1L]]),
    byrow    = TRUE, 
    dimnames = list(rss_sel)
)


# ------------------------------------------------------------------------------
# Step 1.2 : Burden quantification using DLNM ----------------------------------
# ------------------------------------------------------------------------------


# Combine all simulations.
data_sim <- do.call(rbind, simul_list)

# Compute grand total for Quebec by simul.
data_sim <- rbind(data_sim,
    data.table(RSS = 99L, data_sim[, lapply(.SD, sum), .SDcols = cols_sim, by = "TRANGE"])
)

# Compute results in terms of AN.
an_tbl <- do.call(rbind, lapply(c("heat_all", "heat_extreme"), function(trange) {
    
    # Compute mean, lower and higher value.
    do.call(rbind, lapply(c(rss_sel, 99L), function(rss) {
        
        # Extract values of interest.
        val <- ul(data_sim[RSS == rss & TRANGE == trange, ..cols_sim])
        
        # Extract mean values and quantiles.
        return(data.table(
            RSS     = rss,
            TRANGE  = trange,
            AN      = mean(val),
            AN_LOW  = quantile(val, probs = 0.025),
            AN_HIGH = quantile(val, probs = 0.975)
        ))
        
    }))
    
}))

# Count number per year.
count_peryear <- data_final[, .(COUNT = sum(COUNT)/length(unique(YEAR))), by = "RSS"]
count_peryear <- rbind(count_peryear, data.table(RSS = 99L, COUNT = sum(count_peryear$COUNT)))

# Compute AF results.
an_tbl <- merge(an_tbl, count_peryear, by = "RSS", all.x = TRUE)
an_tbl[, `:=`(AF = AN / COUNT, AF_LOW = AN_LOW / COUNT, AF_HIGH = AN_HIGH / COUNT)]
an_tbl[, COUNT := NULL]

# Note : 'an_tbl' contains AN/AF by RSS for total/extreme heat (regional DLNM).
# an_tbl


# ------------------------------------------------------------------------------
# Step 2.1 : Regional BLUPs from meta-regression -------------------------------
# ------------------------------------------------------------------------------


# Fit meta-regression with stepwise variable selection.
meta_reg <- fit_stepwise_meta_aic(
    coef       = coef_matrix, 
    vcov       = vcov_list,
    meta_preds = names(meta_data)[-1]
)
meta_reg_summ <- summary(meta_reg)

# I2 statistic.
meta_reg_summ$i2stat[1L]

# P-value for the Cochran Q-test.
meta_reg_summ$qstat$pvalue[1L]

# Formula of the meta-regression
meta_reg$formula

# Extract BLUPs.
blup_meta <- blup(meta_reg, vcov = TRUE) 

# Appending results to the following vectors.
qmin_rss <- rep(NA, length(rss_sel))
mmt_rss <- rep(NA, length(rss_sel))

# Extract MMT based on BLUP for each region.
for (rss in rss_sel){ 
    
    # Filter data based on <rss>.
    data_rss <- data_final[RSS == rss, ]
    rss_i <- which(rss_sel == rss)
    
    # Extract tvar values to extract MMT based on q_mmt.
    predvar <- quantile(data_rss[[tvar]], probs = q_mmt) 
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        x     = predvar,            
        fun   = "ns", 
        knots = quantile(data_rss[[tvar]], probs = knots_xvar), 
        Bound = range(data_rss[[tvar]]) 
    ) 
    
    # Apply one basis to the tvar.
    bvar <- do.call(onebasis, argvar_tvar) 
    
    # Extract quantile of MMT and MMT. 
    qmin_rss[rss_i] <- q_mmt[which.min((bvar %*% blup_meta[[rss_i]]$blup))]
    mmt_rss[rss_i] <- quantile(data_rss[[tvar]], qmin_rss[rss_i]) 
    
}

# Extract BLUPs of reduced functions at MMT for each RSS.
red_blups_mmt <- lapply(rss_sel, function(rss) {
    
    # Filter data based on <rss>.
    data_rss <- data_final[RSS == rss, ]
    rss_i <- which(rss_sel == rss)
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_rss[[tvar]], probs = knots_xvar), 
        Bound = range(data_rss[[tvar]])
    ) 
    
    # Extract coef and vcov from the BLUPs of meta-regression.
    coef <- blup_meta[[rss_i]]$blup
    vcov <- blup_meta[[rss_i]]$vcov
    
    # Extract MMT.
    mmt <- mmt_rss[rss_i]
    
    # Centered basis functions.
    bvar <- do.call(onebasis, c(list(x = data_rss[[tvar]]), argvar_tvar))
    cenvec <- do.call(onebasis, c(list(x = mmt), argvar_tvar))
    bvarcen <- scale(bvar, center = cenvec, scale = F)
    
    # Predict function values.
    red_blup <- crosspred(bvarcen, coef = coef, vcov = vcov, cen = mmt)
    
})

# Set names.
names(red_blups_mmt) <- rss_sel

# Note : 'red_blups_mmt' contains the BLUPs estimates for each RSS at MMT.
# For example :
# plot(red_blups_mmt[[1L]], xlab = tvar, ylab = "RR")


# ------------------------------------------------------------------------------
# Step 2.2 : Burden quantification with reduced regional BLUPs -----------------
# ------------------------------------------------------------------------------


# Get years.
years <- unique(data_final$YEAR)
nyears <- length(years)

# Compute simulated AN for each RSS based on BLUPs.
data_sim <- do.call(rbind, lapply(rss_sel, function(rss) {
    
    # Message.
    if (verbose) message("Computing AN for region ", rss, ".")
    
    # Filter data based on <rss>.
    data_rss <- data_final[RSS == rss, ]
    rss_i <- which(rss_sel == rss)
    
    # Argument for the tvar cross-basis.
    argvar_tvar <- list(
        fun   = "ns", 
        knots = quantile(data_rss[[tvar]], probs = knots_xvar), 
        Bound = range(data_rss[[tvar]]) 
    ) 
    
    # Extract MMT and threshold for extreme temperature
    mmt <- mmt_rss[rss_i]
    ext <- heat_thresh[RSS == rss, get(tvar)]
    
    # Extract coef and vcov from the meta-regression.
    coef <- blup_meta[[rss_i]]$blup
    vcov <- blup_meta[[rss_i]]$vcov
    
    # Centered basis functions.
    bvar <- do.call(onebasis, c(list(x = data_rss[[tvar]]), argvar_tvar))
    cenvec <- do.call(onebasis, c(list(x = mmt), argvar_tvar))
    bvarcen <- scale(bvar, center = cenvec, scale = F)
    
    # Indicators for heat and extreme heat days.
    data_rss[, heat_all     := get(tvar) > mmt]
    data_rss[, heat_extreme := get(tvar) > max(mmt, ext)]
    
    # Compute the contribution of daily values (average value).
    an_mean <- (1 - exp(-bvarcen %*% coef)) * data_rss[["COUNT"]]
    
    # Mean AN values for heat and extreme heat.
    an_heatall <- sum(an_mean[data_rss$heat_all])/nyears
    an_heatext <- sum(an_mean[data_rss$heat_extreme])/nyears
    
    # Sampling coefficient from BLUP assuming mvnorm.
    set.seed(2912L)
    coefsim <- MASS::mvrnorm(nsim, coef, vcov)
    
    # Compute simulations.
    an_sim <- sapply(1:nsim, function(s) {
        (1 - exp(-bvarcen %*% coefsim[s, ])) * data_rss[["COUNT"]]
    })
    
    # Export data.
    colnames(an_sim) <- cols_sim
    data_sim <- cbind(data_rss[, .(RSS, DATE, YEAR, heat_all, heat_extreme)], an_sim)
    
    # Return data_sim.
    return(data_sim)
    
}))

# Aggregate simulations by years.
x_year <- do.call(rbind, lapply(c("heat_all", "heat_extreme"), function(trange) {
    
    # Compute mean, lower and higher value.
    do.call(rbind, lapply(years, function(year) {
        
        # Extract values of interest.
        val <- ul(data_sim[YEAR == year & get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
        
        # Extract mean values and quantiles.
        return(data.table(
            RSS     = 99L,
            YEAR    = year,
            TRANGE  = trange,
            AN      = mean(val),
            AN_LOW  = quantile(val, probs = 0.025),
            AN_HIGH = quantile(val, probs = 0.975)
        ))
        
    }))
    
}))

# Aggregate simulations by region (RSS).
x_rss <- do.call(rbind, lapply(c("heat_all", "heat_extreme"), function(trange) {
    
    # Compute mean, lower and higher value.
    do.call(rbind, lapply(rss_sel, function(rss) {
        
        # Extract values of interest.
        val <- ul(data_sim[RSS == rss & get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
        
        # Extract mean values and quantiles.
        return(data.table(
            RSS     = rss,
            YEAR    = 9999L,
            TRANGE  = trange,
            AN      = mean(val)/nyears,
            AN_LOW  = quantile(val, probs = 0.025)/nyears,
            AN_HIGH = quantile(val, probs = 0.975)/nyears
        ))
        
    }))
    
}))

# Aggregate simulations for the whole province.
x_qc <- do.call(rbind, lapply(c("heat_all", "heat_extreme"), function(trange) {
    
    # Extract values of interest.
    val <- ul(data_sim[get(trange) == TRUE, lapply(.SD, sum), .SDcols = cols_sim])
    
    # Extract mean values and quantiles.
    return(data.table(
        RSS     = 99L,
        YEAR    = 9999L,
        TRANGE  = trange,
        AN      = mean(val)/nyears,
        AN_LOW  = quantile(val, probs = 0.025)/nyears,
        AN_HIGH = quantile(val, probs = 0.975)/nyears
    ))
    
}))

# Merge all tables.
an_tbl <- rbind(x_qc, x_rss, x_year)

# Count number per year and region (RSS).
count_year_rss <- rbind(
    data_final[, .(COUNT = sum(COUNT)/length(unique(YEAR)), YEAR = 9999L), by = c("RSS")],
    data_final[, .(COUNT = sum(COUNT), RSS = 99L), by = c("YEAR")],
    data_final[, .(COUNT = sum(COUNT)/length(unique(YEAR)), YEAR = 9999L, RSS = 99L)]
)

# Compute AF results.
an_tbl <- merge(an_tbl, count_year_rss, by = c("RSS", "YEAR"), all.x = TRUE)
an_tbl[, `:=`(AF = AN / COUNT, AF_LOW = AN_LOW / COUNT, AF_HIGH = AN_HIGH / COUNT)]
an_tbl[, COUNT := NULL]

# Note : 'an_tbl' contains AN/AF by RSS/Year based on regional BLUPS.
# an_tbl
