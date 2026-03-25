###############################################################
# CLEAN ENVIRONMENT
###############################################################
rm(list = ls())

###############################################################
# LOAD SIMULATION
###############################################################
sim <- readRDS("Rdata/simlations.rds")

df_long <- sim$df_long
true    <- sim$true

# truth components saved in the updated simulation
true_detection_raw      <- true$true_detection_raw
true_detection_contrast <- true$true_detection_contrast
true_spatial            <- true$true_spatial
true_year_coefficients  <- true$true_year_coefficients
true_year_totals        <- true$true_year_totals   # <-- includes year, Z_year_true, N_total_year

###############################################################
# 1) YEARLY POPULATION TREND (using saved truth)
###############################################################

years        <- true_year_totals$year
Z_year_true  <- true_year_totals$Z_year_true
N_total_year <- true_year_totals$N_total_year

# slope & intercept from saved truth
N_intercept_true <- true_year_coefficients$value[
  true_year_coefficients$parameter == "intercept_linear_true"
]

N_slope_true <- true_year_coefficients$value[
  true_year_coefficients$parameter == "slope_linear_true"
]

# expected linear trend
expected <- N_intercept_true + N_slope_true * Z_year_true

# plot truth vs expected
plot(years, N_total_year,
     pch = 19, col = "black",
     xlab = "Year", ylab = "Total population",
     main = "Yearly total population vs expected linear trend")

lines(years, expected, lty = 2, lwd = 2, col = "red")

###############################################################
# 2) DETECTION vs CONTINUOUS COVARIATES (x1, x2, x3)
###############################################################

p_true     <- df_long$p_true
beta0_true <- true_detection_raw$value[true_detection_raw$name == "beta0_raw"]

### X1 --------------------------------------------------------
x1 <- df_long$x1
beta_x1_true <- true_detection_raw$value[true_detection_raw$name == "x1"]

plot(x1, p_true,
     pch = 20, col = rgb(0,0,0,0.2),
     xlab = "x1 (scaled)", ylab = "p_true",
     main = "Detection probability vs x1")

x1_grid <- seq(min(x1), max(x1), length.out = 200)
p_hat_x1 <- plogis(beta0_true + beta_x1_true * x1_grid)
lines(x1_grid, p_hat_x1, col = "red", lwd = 2)

### X2 --------------------------------------------------------
x2 <- df_long$x2
beta_x2_true <- true_detection_raw$value[true_detection_raw$name == "x2"]

plot(x2, p_true,
     pch = 20, col = rgb(0,0,0,0.2),
     xlab = "x2 (scaled)", ylab = "p_true",
     main = "Detection probability vs x2")

x2_grid <- seq(min(x2), max(x2), length.out = 200)
p_hat_x2 <- plogis(beta0_true + beta_x2_true * x2_grid)
lines(x2_grid, p_hat_x2, col = "red", lwd = 2)

### X3 --------------------------------------------------------
x3 <- df_long$x3
beta_x3_true <- true_detection_raw$value[true_detection_raw$name == "x3"]

plot(x3, p_true,
     pch = 20, col = rgb(0,0,0,0.2),
     xlab = "x3 (scaled)", ylab = "p_true",
     main = "Detection probability vs x3 (true effect = 0)")

x3_grid <- seq(min(x3), max(x3), length.out = 200)
p_hat_x3 <- plogis(beta0_true + beta_x3_true * x3_grid)
lines(x3_grid, p_hat_x3, col = "red", lwd = 2)

###############################################################
# 3) MONTH EFFECT CHECK (saved raw seasonal effects)
###############################################################

month_levels <- sort(unique(df_long$month_f))

boxplot(p_true ~ df_long$month_f,
        las = 2,
        main = "Detection probability vs month",
        xlab = "Month", ylab = "p_true")

# theoretical seasonal raw coefficients
month_names <- paste0("month_", month_levels)
beta_month_raw <- true_detection_raw$value[
  match(month_names, true_detection_raw$name)
]

points(1:length(month_levels),
       plogis(beta0_true + beta_month_raw),
       col = "red", pch = 19)

###############################################################
# 4) FIXEDF2 EFFECT CHECK (saved raw values)
###############################################################

fixed_levels <- sort(unique(df_long$fixedF2))

boxplot(p_true ~ df_long$fixedF2,
        main = "Detection probability vs fixedF2",
        xlab = "fixedF2 level", ylab = "p_true")

fixed_names <- paste0("fixedF2_", fixed_levels)
beta_F2_raw <- true_detection_raw$value[
  match(fixed_names, true_detection_raw$name)
]

points(1:length(fixed_levels),
       plogis(beta0_true + beta_F2_raw),
       col = "red", pch = 19)

#
###############################################################
# 5) HEATMAP OF w (TRUE SPATIAL WEIGHTS)
###############################################################

# spatial truth table has rows for locations 1..K
w <- true_spatial$w
# determine grid dimensions automatically
coords_x <- sort(unique(true_spatial$x_coord))
coords_y <- sort(unique(true_spatial$y_coord))

nx <- length(coords_x)
ny <- length(coords_y)

# reshape into ny x nx matrix (row-wise fill, matching simulation)
w_mat <- matrix(w, nrow = ny, ncol = nx, byrow = TRUE)

image(coords_x, coords_y, w_mat,
      col = heat.colors(50),
      xlab = "x coordinate", ylab = "y coordinate",
      main = "Heatmap of true spatial weights w")
contour(coords_x, coords_y, w_mat, add = TRUE)

###############################################################
# 5) HEATMAP OF land
###############################################################

# reshape into ny x nx matrix (row-wise fill, matching simulation)
prop_land_mat <- matrix(true_spatial$prop_land, nrow = ny, ncol = nx, byrow = TRUE)

image(coords_x, coords_y, prop_land_mat,
      col = heat.colors(50),
      xlab = "x coordinate", ylab = "y coordinate",
      main = "Heatmap of true spatial weights w")
contour(coords_x, coords_y, prop_land_mat, add = TRUE)


###############################################################
# 6) HEATMAP OF n_latent  (FIXED)
###############################################################

# choose a year (first year)
chosen_year <- min(df_long$year)

# extract rows for that year
idx <- which(df_long$year == chosen_year)

# For that year, n_latent has length = 100 sites * 12 months = 1200
n_lat_year <- df_long$n_latent[idx]
loc_year   <- df_long$location[idx]

# Collapse to one value per site:
# Option A: mean across months
n_lat_site <- tapply(n_lat_year, loc_year, mean)

# reshape into matrix ny × nx
nlat_mat <- matrix(n_lat_site, nrow = ny, ncol = nx, byrow = TRUE)

# plot
image(coords_x, coords_y, nlat_mat,
      col = heat.colors(50),
      xlab = "x coordinate", ylab = "y coordinate",
      main = paste("Heatmap of n_latent (summed over months, year =", chosen_year, ")"))
contour(coords_x, coords_y, nlat_mat, add = TRUE)

###############################################################
# END OF DIAGNOSTICS
###############################################################
