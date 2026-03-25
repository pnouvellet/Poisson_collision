plot_detection_df_long <- function(df_long, p_det_med, p_det_low, p_det_high) {
  library(Hmisc)
  
  p_true <- df_long$p_true
  
  plot(p_true, p_det_med,
       pch = 16,
       col = rgb(0,0,0,0.4),
       xlab = "TRUE p_true",
       ylab = "Posterior p_hat (median)",
       main = "TRUE vs ESTIMATED detection probability")
  
  errbar(p_true, p_det_med, p_det_high, p_det_low,
         add=TRUE,
         pch=16,
         col=rgb(0,0,1,0.3),
         errbar.col=rgb(0,0,1,0.3))
  
  abline(0, 1, col="red", lwd=2, lty=2)
  
  legend("topleft",
         legend=c("Posterior median + 95% CI", "1:1 line"),
         col=c("blue","red"),
         pch=c(16, NA),
         lty=c(0, 2),
         bty="n")
}