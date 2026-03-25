# ===============================================================
# FUNCTION: plot_detection_curve
# Plots:
#   1) true p_true vs covariate x
#   2) posterior marginal prediction (mean + 95% CI)
#   3) true logistic curve from raw simulation parameters
#
# Arguments:
#   cov_name  = "x1", "x2", or "x3"
#   df_long   = simulated long dataframe
#   true      = truth list (from simulation_12.r)
#   post_mat  = posterior draws matrix (as_draws_matrix)
#   dat       = Stan data list (prep$dat)
#
# ===============================================================

plot_detection_curve <- function(cov_name, df_long, true, post_mat, dat) {
  
  # ----------------------------
  # Extract covariate and p_true
  # ----------------------------
  x <- df_long[[cov_name]]
  p_true <- df_long$p_true
  
  # ----------------------------
  # BASIC PLOT OF TRUE p_true
  # ----------------------------
  plot(x, p_true,
       pch = 20,
       col = rgb(0,0,0,0.2),
       xlab = paste0(cov_name, " (scaled)"),
       ylab = "Detection probability",
       main = paste("Detection vs", cov_name,
                    ": true p_true, posterior prediction, true logistic curve"))
  
  # ----------------------------
  # POSTERIOR PARAMETER DRAWS
  # ----------------------------
  beta0_draws <- post_mat[, "beta0"]
  
  beta_cols  <- grep("^beta\\[", colnames(post_mat))
  beta_draws <- post_mat[, beta_cols, drop = FALSE]
  
  beta_names <- colnames(dat$X)
  
  # which column corresponds to this covariate?
  cov_col <- which(beta_names == cov_name)
  if (length(cov_col) != 1)
    stop(paste("Cannot find", cov_name, "in design matrix column names."))
  
  # ----------------------------
  # Prediction grid
  # ----------------------------
  x_grid <- seq(min(x), max(x), length.out = 200)
  
  n_draws <- length(beta0_draws)
  n_grid  <- length(x_grid)
  
  eta_post <- matrix(NA, n_draws, n_grid)
  
  # SAFE, reliable computation: eta = beta0 + beta_cov * x
  for (j in 1:n_grid) {
    eta_post[, j] <- beta0_draws + beta_draws[, cov_col] * x_grid[j]
  }
  
  p_post <- plogis(eta_post)
  
  # posterior summaries
  p_post_mean <- apply(p_post, 2, mean)
  p_post_low  <- apply(p_post, 2, quantile, 0.025)
  p_post_high <- apply(p_post, 2, quantile, 0.975)
  
  # CI band
  polygon(
    x = c(x_grid, rev(x_grid)),
    y = c(p_post_low, rev(p_post_high)),
    col = rgb(0, 0, 1, 0.15),
    border = NA
  )
  
  # posterior mean
  lines(x_grid, p_post_mean, col = "blue4", lwd = 2)
  
  # ----------------------------
  # TRUE logistic curve
  # ----------------------------
  beta0_true <- true$true_detection_raw$value[
    true$true_detection_raw$name == "beta0_raw"
  ]
  
  beta_cov_true <- true$true_detection_raw$value[
    true$true_detection_raw$name == cov_name
  ]
  
  p_true_curve <- plogis(beta0_true + beta_cov_true * x_grid)
  
  lines(x_grid, p_true_curve, col = "red", lwd = 2, lty = 2)
  
  # ----------------------------
  # LEGEND
  # ----------------------------
  legend("topleft",
         legend = c("True p_true", "Posterior mean", "Posterior 95% CI", "True logistic curve"),
         col = c(rgb(0,0,0,0.4), "blue4", rgb(0,0,1,0.15), "red"),
         pch = c(20, NA, 15, NA),
         lwd = c(NA, 2, NA, 2),
         lty = c(0, 1, 0, 2),
         pt.cex = c(1, NA, 2, NA),
         bty = "n")
}
