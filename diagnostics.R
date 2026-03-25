###############################################################
# CLEAN ENVIRONMENT
###############################################################
rm(list=ls())

library(Hmisc)
library(posterior)
library(dplyr)
library(loo)

###############################################################
# LOAD SIMULATION AND FIT RESULTS
###############################################################
sim  <- readRDS("Rdata/simlations.rds")
fitR_full <- readRDS("Rdata/fitting_full.rds")
fitR_reduced <- readRDS("Rdata/fitting_reduced.rds")

# compare models
fitR_full$perf$loo_result
fitR_reduced$perf$loo_result

loo_compare(fitR_full$perf$loo_result,
            fitR_reduced$perf$loo_result)

show_full <- FALSE 

if(show_full){
  fitR <- fitR_full
}else{
  fitR <- fitR_reduced
}
# prepare for diagnostic of best model

df_long <- sim$df_long
true    <- sim$true
post_mat <- as_draws_matrix(fitR$fit$draws())
dat <- fitR$prep$dat

###############################################################
# 1. YEARLY POPULATION TREND (simple, unchanged)
###############################################################
Y <- dat$Y
nyears <- 1:Y
N_cols <- paste0("N_total_year[", nyears, "]")

post_N <- post_mat[, N_cols, drop=FALSE]
N_med  <- apply(post_N, 2, median)
N_low  <- apply(post_N, 2, quantile,0.025)
N_high <- apply(post_N, 2, quantile,0.975)

plot(nyears, N_med,
     ylim=range(c(N_low,N_high,true$true_year_totals$N_total_year)),
     pch=16, col="blue4",
     xlab="Year", ylab="Population size",
     main="Yearly population trend")
errbar(nyears, N_med, N_high, N_low, add=TRUE, col="blue4")
points(nyears,true$true_year_totals$N_total_year,col="red",pch=16)
lines(nyears,true$true_year_totals$N_total_year,col="red",lwd=2)

###############################################################
# 2. DETECTION: p_true vs p_hat
###############################################################
source("plot_detection_df_long.R")

# Extract posterior p_det draws
p_cols <- grep("^p_det\\[", colnames(post_mat))
post_p <- post_mat[, p_cols, drop=FALSE]
p_det_med <- apply(post_p,2,median)
p_det_low <- apply(post_p,2,quantile,0.025)
p_det_high<- apply(post_p,2,quantile,0.975)

# f <- which((df_long$year %in% c(1)) & (df_long$month_f %in% 'jul') & (df_long$fixedF2 %in% 'b'))
# f <- which((df_long$year %in% c(1)) & (df_long$month_f %in% 'jul') )
f <- which((df_long$year %in% c(2))  )
plot_detection_df_long(df_long[f,], p_det_med[f], p_det_low[f], p_det_high[f])

# d_temp <- cbind(df_long[f,], p_det_med[f], p_det_low[f], p_det_high[f])
###############################################################
# 3. CTN + FACTOR EFFECTS (month + fixedF2)
###############################################################
source("plot_factor_effect_stan.R")

plot_factor_effect_stan(df_long, "month_f", p_det_med, p_det_low, p_det_high)
plot_factor_effect_stan(df_long, "fixedF2", p_det_med, p_det_low, p_det_high)

source("plot_detection_curve.R")

# x1 diagnostics
plot_detection_curve("x1", df_long, true, post_mat, dat)
# x2 diagnostics
plot_detection_curve("x2", df_long, true, post_mat, dat)
if(show_full){
  # x3 diagnostics
  plot_detection_curve("x3", df_long, true, post_mat, dat)
}
###############################################################
# 4. N_LATENT DIAGNOSTICS USING STAN GQ OUTPUT
###############################################################
source("plot_n_latent_one_year.R")

# Extract n_latent from Stan
K <- dat$K
n_cols <- grep("^n_latent\\[", colnames(post_mat))
post_n <- post_mat[, n_cols, drop=FALSE]

# reshape Stan output into K × Y
n_lat_array <- array(NA, c(K,Y))
col_ids <- expand.grid(site=1:K, year=1:Y)

for (j in 1:ncol(post_n)) {
  s <- col_ids$site[j]
  y <- col_ids$year[j]
  n_lat_array[s,y] <- median(post_n[,j])
}

# But we want full summaries
n_lat_med  <- matrix(NA, K, Y)
n_lat_low  <- matrix(NA, K, Y)
n_lat_high <- matrix(NA, K, Y)

for (j in 1:ncol(post_n)) {
  s <- col_ids$site[j]
  y <- col_ids$year[j]
  n_lat_med[s,y]  <- median(post_n[,j])
  n_lat_low[s,y]  <- quantile(post_n[,j],0.025)
  n_lat_high[s,y] <- quantile(post_n[,j],0.975)
}

# Plot per year
for (y in 1){#:Y) {
  plot_n_latent_one_year(y, df_long, true, n_lat_med, n_lat_low, n_lat_high)
}

