###############################################################
# CLEANED INFERENCE SCRIPT
# - Uses uninformative priors on detection coefficients
# - No contrast-coded prior work
# - Minimal and consistent with simulation + make_stan_data
###############################################################

# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# cmdstanr::install_cmdstan()
library(cmdstanr)
library(posterior)
library(dplyr)
library(loo)


###############################################################
# CLEAN ENVIRONMENT
###############################################################
rm(list = ls())
source("make_stan_data.R") # load prep data function

###############################################################
# LOAD SIMULATION
###############################################################
sim <- readRDS("Rdata/simlations.rds")

df_long <- sim$df_long

###############################################################
# 1) SIMPLE BETA0 GUESS (MATCHING YOUR ORIGINAL METHOD)
###############################################################

Ey   <- mean(df_long$y)
E_Nw <- mean(df_long$N_total_year * df_long$w)
p_detect_guess <- (Ey / E_Nw)
beta0_guess <- log(p_detect_guess)

###############################################################
# 2) PREPARE STAN DATA (WITH UNINFORMATIVE PRIORS)
###############################################################

prep <- make_stan_data_formula_sumcoded_icar(
  df_long,
  detection_formula = "~ x1 + x2 + month_f ",
  year_formula      = "~ year",
  
  # year totals priors
  N_overall_mean = mean(df_long$N_total_year),
  N_overall_sd   = 150,
  
  # detection priors
  beta0_mean = beta0_guess,
  beta0_sd   = 1,
  
  # UNINFORMATIVE PRIORS ON ALL OTHER DETECTION BETAS:
  beta_mean = NULL,      # means = rep(0, P)
  beta_sd   = 2,         # weak/flat prior
  
  tau_prior = 1.0,
  verbose = TRUE
)

###############################################################
# 3) FIT STAN MODEL
###############################################################

mod <- cmdstan_model("Poisson_model.stan", force_recompile = FALSE)

fit <- mod$sample(
  data = prep$dat,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.85,
  max_treedepth = 11
)

###############################################################
# 4) SAVE POSTERIOR
###############################################################

post_mat <- as_draws_matrix(fit$draws())

###############################################################
# SUMMARY DATAFRAME FOR n_latent  (K × Y array)
###############################################################

# 1. Extract all n_latent draws from post_mat
nl_cols <- grep("^n_latent\\[", colnames(post_mat))
post_n_lat <- as.matrix(post_mat[, nl_cols, drop = FALSE])

# Number of sites and years
K <- length(unique(df_long$location))     # 100
Y <- length(unique(df_long$year))     # 5

# Sanity check: number of columns must be K*Y
stopifnot(ncol(post_n_lat) == K * Y)

# 2. Build index mapping: column j corresponds to (site, year)
# Column order in Stan is row-major: n_latent[1,1], n_latent[2,1], ..., n_latent[K,1], n_latent[1,2]...
col_map <- expand.grid(
  site = 1:K,
  year = 1:Y
)

# 3. Compute posterior summaries
median_vec <- apply(post_n_lat, 2, median)
low_vec    <- apply(post_n_lat, 2, quantile, 0.025)
high_vec   <- apply(post_n_lat, 2, quantile, 0.975)

# 4. Build final dataframe
n_latent_summary <- data.frame(
  site   = col_map$site,
  year   = col_map$year,
  median = median_vec,
  low95  = low_vec,
  high95 = high_vec
)

# Optional: sort by site, year
n_latent_summary <- n_latent_summary[order(n_latent_summary$year,
                                           n_latent_summary$site), ]

###############################################################
# SUMMARY DATAFRAME FOR p_detection (size N = rows of df_long)
###############################################################

# 1. Extract posterior draws for p_det
p_cols <- grep("^p_det\\[", colnames(post_mat))
post_pdet <- as.matrix(post_mat[, p_cols, drop = FALSE])

# N = number of observations = nrow(df_long)
N <- nrow(df_long)

# Sanity check
stopifnot(ncol(post_pdet) == N)

# 2. Posterior summaries per observation
p_med  <- apply(post_pdet, 2, median)
p_low  <- apply(post_pdet, 2, quantile, 0.025)
p_hig  <- apply(post_pdet, 2, quantile, 0.975)

# 3. Combine with site, year, month from df_long
p_detection_summary <- data.frame(
  row   = 1:N,
  site  = df_long$location,
  year  = df_long$year,
  month = df_long$month_f,
  median = p_med,
  low95  = p_low,
  high95 = p_hig
)

# Optionally reorder rows
p_detection_summary <- p_detection_summary[
  order(p_detection_summary$year,
        p_detection_summary$site,
        p_detection_summary$month),
]

# model performance

# post_mat <- as_draws_matrix(fitR$fit$draws())
loglik_cols <- grep("^log_lik\\[", colnames(post_mat))
log_lik <- post_mat[, loglik_cols]


loo_result  <- loo(log_lik)
waic_result <- waic(log_lik)

perf <- list(loo_result = loo_result,
             waic_result = waic_result)

# saving

results <- list(prep = prep,
                fit = fit,
                post_mat = post_mat,
                p_detection_summary = p_detection_summary,
                n_latent_summary = n_latent_summary,
                perf = perf)

saveRDS(object = results, file = 'Rdata/fitting_reduced.rds')
