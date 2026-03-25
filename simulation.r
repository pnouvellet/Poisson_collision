###############################################################
# FULL SIMULATION SCRIPT (LINEAR YEAR MODEL + CONTRAST DETECTION)
# WITH prop_land INTEGRATED
###############################################################
rm(list = ls())
library(dplyr)
library(ggplot2)
set.seed(156)
set.seed(14)


###############################################################
# 1) BASIC SETTINGS
###############################################################
grid_nx <- 10
grid_ny <- 10
K <- grid_nx * grid_ny
J <- 12
Y <- 5
month_levels <- tolower(month.abb)
levels_fixedF2 <- letters[1:5]

###############################################################
# 2) TRUE YEARLY POPULATION MODEL
###############################################################
year_numeric <- 1:Y
Z_year_true <- scale(year_numeric, center = TRUE, scale = TRUE)[,1]

N_intercept_true <- 4000
N_slope_true <- 300
N_total_year <- N_intercept_true + N_slope_true * Z_year_true

true_year_coefficients <- data.frame(
  parameter = c("intercept_linear_true","slope_linear_true"),
  value = c(N_intercept_true, N_slope_true)
)

true_year_totals <- data.frame(
  year = year_numeric,
  Z_year_true = Z_year_true,
  N_total_year = N_total_year
)

###############################################################
# 3) SPATIAL GAUSSIAN MIXTURE KERNEL
###############################################################
n_peaks <- 3
peak_heights <- c(2,3,2)
peak_widths <- rep(0.15, n_peaks)
peak_positions <- cbind(runif(n_peaks), runif(n_peaks))
kernel_floor <- 1e-6

xs <- seq(0, 1, length.out = grid_nx)
ys <- seq(0, 1, length.out = grid_ny)
x_coord <- rep(xs, times = grid_ny)
y_coord <- rep(ys, each = grid_nx)
coords <- cbind(x_coord, y_coord)
intensity <- rep(0, K)

for (p in 1:n_peaks) {
  dx <- coords[,1] - peak_positions[p,1]
  dy <- coords[,2] - peak_positions[p,2]
  r2 <- dx^2 + dy^2
  intensity <- intensity + peak_heights[p] * exp(-0.5 * r2 / peak_widths[p]^2)
}
intensity <- pmax(intensity, kernel_floor)

###############################################################
# 3b) ADD prop_land
###############################################################

prop_land <- rep(1, K)
idx <- sample(1:K, size = round(K/3), replace = FALSE)
prop_land[idx] <- runif(length(idx), min = 0.1, max = 0.5)

###############################################################
# Land-adjusted spatial weights
###############################################################
w_raw <- intensity / sum(intensity)
w_land <- w_raw * prop_land
w_land <- w_land / sum(w_land)
w <- w_land

# Latent per site-year
n_latent <- sapply(N_total_year, function(Ny) round(Ny * w))

###############################################################
# 4) DESIGN GRID
###############################################################
location <- rep(1:K, each = J*Y)
year <- rep(rep(1:Y, each = J), times = K)
month_id <- rep(rep(1:J, times = Y), times = K)
month_f <- factor(month_levels[month_id], levels = month_levels, ordered = TRUE)

###############################################################
# 5) COVARIATES
###############################################################
N_obs <- length(location)
x1 <- rnorm(N_obs)
x2 <- rnorm(N_obs)*0.5 + 0.2
x3 <- rnorm(N_obs)

# NEW: adjust x1 and x2 per land area, keeping their names
# land-area is site-level, so match via location index
x1 <- x1 / prop_land[location]   # or prop_land[location] if df_long not yet built
x2 <- x2 / prop_land[location]

x1s <- scale(x1)[,1]
x2s <- scale(x2)[,1]
x3s <- scale(x3)[,1]

###############################################################
# 6) FIXEDF2 FACTOR
###############################################################
fixedF2_by_loc <- sample(levels_fixedF2, K, replace = TRUE)
fixedF2 <- factor(fixedF2_by_loc[location], levels = levels_fixedF2)

###############################################################
# 7) BUILD df_long (INCLUDING prop_land)
###############################################################
df_long <- data.frame(
  location = location,
  year = year,
  month_f = month_f,
  x1 = x1s,
  x2 = x2s,
  x3 = x3s,
  N_total_year = N_total_year[year],
  n_latent = n_latent[cbind(location, year)],
  fixedF2 = fixedF2,
  x_coord = x_coord[location],
  y_coord = y_coord[location],
  w = w[location],
  prop_land = prop_land[location]
)

###############################################################
# 8) RAW DETECTION COEFFICIENTS
###############################################################
beta0_true <- -4
beta_x_raw <- c(x1 = 0.5, x2 = -0.3, x3 = 0.0)

season_amp <- 1.0
season_beta_raw <- season_amp * sin(2*pi*((1:J)-1)/J)
names(season_beta_raw) <- month_levels

beta_fixedF2_raw <- setNames(rep(0,length(levels_fixedF2)), levels_fixedF2)

###############################################################
# 9) DETECTION DESIGN MATRIX
###############################################################
detection_formula <- ~ x1 + x2 + x3 + month_f + fixedF2
mf_det <- model.frame(detection_formula, df_long)

is_num_det <- sapply(mf_det, is.numeric)
for (nm in names(mf_det)[is_num_det]) {
  sdv <- sd(mf_det[[nm]]); if (sdv <= 0) sdv <- 1
  mf_det[[nm]] <- (mf_det[[nm]] - mean(mf_det[[nm]])) / sdv
}

det_fac_names <- names(mf_det)[sapply(mf_det, is.factor)]
det_contr <- if (length(det_fac_names)) {
  setNames(
    lapply(det_fac_names, function(nm) contr.sum(nlevels(mf_det[[nm]]))),
    det_fac_names)
} else NULL

X_full <- model.matrix(attr(mf_det, "terms"), mf_det, contrasts.arg = det_contr)
X <- X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE]

###############################################################
# 10) TRUE CONTRAST-CODED COEFFICIENTS
###############################################################
beta_true_named <- rep(NA, ncol(X))
names(beta_true_named) <- colnames(X)

# fixedF2 contrasts
K_F2 <- length(levels_fixedF2)
C_F2 <- contr.sum(K_F2)
rownames(C_F2) <- levels_fixedF2
raw_F2 <- beta_fixedF2_raw - mean(beta_fixedF2_raw)
b_F2 <- qr.solve(C_F2, raw_F2)

# month_f contrasts
K_m <- length(month_levels)
C_m <- contr.sum(K_m)
rownames(C_m) <- month_levels
raw_month <- season_beta_raw - mean(season_beta_raw)
b_month <- qr.solve(C_m, raw_month)

for (nm in colnames(X)) {
  if (nm %in% c("x1","x2","x3")) {
    beta_true_named[nm] <- beta_x_raw[nm]; next
  }
  if (grepl("^month_f", nm)) {
    j <- as.integer(sub("^month_f","",nm))
    beta_true_named[nm] <- b_month[j]
    next
  }
  if (grepl("^fixedF2", nm)) {
    j <- as.integer(sub("^fixedF2","",nm))
    beta_true_named[nm] <- b_F2[j]
    next
  }
  stop("Unexpected design column: ", nm)
}

###############################################################
# 11) LINEAR PREDICTOR AND PROBABILITY
###############################################################
eta_true <- beta0_true + as.numeric(X %*% beta_true_named)
p_true <- plogis(eta_true)
mu_true <- df_long$n_latent * p_true

df_long$eta_true <- eta_true
df_long$p_true <- p_true
df_long$mu_true <- mu_true

###############################################################
# 12) OBSERVE COUNTS
###############################################################
df_long$y <- rbinom(N_obs, df_long$n_latent, df_long$p_true)

###############################################################
# 13) TRUTH TABLES
###############################################################
true_detection_contrast <- data.frame(
  name = c("beta0", names(beta_true_named)),
  value = c(beta0_true, beta_true_named)
)

true_detection_raw <- data.frame(
  name = c("beta0_raw",
           names(beta_x_raw),
           paste0("month_", month_levels),
           paste0("fixedF2_", levels_fixedF2)),
  value = c(beta0_true,
            beta_x_raw,
            season_beta_raw,
            beta_fixedF2_raw)
)

true_spatial <- data.frame(
  location = 1:K,
  w = w,
  prop_land = prop_land,
  x_coord = x_coord,
  y_coord = y_coord
)

###############################################################
# END OF SIMULATION
###############################################################
true <- list(
  true_detection_contrast = true_detection_contrast,
  true_detection_raw = true_detection_raw,
  true_spatial = true_spatial,
  true_year_coefficients = true_year_coefficients,
  true_year_totals = true_year_totals
)

simulation <- list(df_long = df_long,
                   true = true)

saveRDS(simulation, file = "Rdata/simlations.rds")
