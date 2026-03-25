plot_factor_effect_stan <- function(df_long, factor_name, p_det_med, p_det_low, p_det_high) {
  library(Hmisc)
  
  fac <- df_long[[factor_name]]
  fac_levels <- levels(fac)
  
  # TRUE p_true by factor
  p_true <- df_long$p_true
  boxplot(p_true ~ fac,
          las=2,
          xlab = factor_name,
          ylab = "TRUE detection probability",
          main = paste0("Effect of ", factor_name, ": TRUE p_true + posterior CIs"))
  
  # posterior p_det by factor
  n_levels <- length(fac_levels)
  xpos <- 1:n_levels
  
  post_m <- tapply(p_det_med, fac, median)
  post_l <- tapply(p_det_low, fac, median)
  post_h <- tapply(p_det_high, fac, median)
  
  errbar(xpos+0.1, post_m, post_h, post_l,
         add=TRUE,
         pch=15,
         col="blue4",
         errbar.col="blue4")
  
  # true effect from simulation (raw)
  true_raw <- true$true_detection_raw
  beta0_true <- true_raw$value[true_raw$name == "beta0_raw"]
  
  if (factor_name == "month_f") {
    true_names <- paste0("month_", fac_levels)
  } else {
    true_names <- paste0(factor_name, "_", fac_levels)
  }
  
  fac_raw <- true_raw$value[match(true_names, true_raw$name)]
  p_true_curve <- plogis(beta0_true + fac_raw)
  
  points(xpos - 0.1, p_true_curve, pch=16, col="red", cex=1.2)
  
  legend("topright",
         legend=c("Posterior median + 95% CI", "TRUE factor effect"),
         col=c("blue4","red"),
         pch=c(15,16),
         bty="n")
}