plot_n_latent_one_year <- function(year_idx, df_long, true, n_lat_med, n_lat_low, n_lat_high) {
  library(Hmisc)
  
  true_sp <- true$true_spatial
  w_true  <- true_sp$w
  
  coords_x <- sort(unique(true_sp$x_coord))
  coords_y <- sort(unique(true_sp$y_coord))
  nx <- length(coords_x)
  ny <- length(coords_y)
  
  K <- length(w_true)
  
  # TRUE n_latent
  true_N_y <- true$true_year_totals$N_total_year[year_idx]
  true_n <- true_N_y * w_true
  true_mat <- matrix(true_n, ny, nx, byrow=TRUE)
  
  image(coords_x, coords_y, true_mat,
        col = heat.colors(50),
        main = paste("TRUE n_latent (year =", year_idx, ")"),
        xlab="x", ylab="y")
  contour(coords_x, coords_y, true_mat, add=TRUE)
  
  # Estimated
  est_n <- n_lat_med[, year_idx]
  est_mat <- matrix(est_n, ny, nx, byrow=TRUE)
  
  image(coords_x, coords_y, est_mat,
        col = heat.colors(50),
        main = paste("Estimated n_latent (median, year =", year_idx, ")"),
        xlab="x", ylab="y")
  contour(coords_x, coords_y, est_mat, add=TRUE)
  
  # Scatter
  plot(true_n, est_n,
       pch=16,
       col=rgb(0,0,0,0.5),
       xlab="TRUE n_latent",
       ylab="Posterior median n_latent",
       main=paste("TRUE vs ESTIMATED n_latent (year =", year_idx, ")"))
  
  errbar(true_n, est_n, n_lat_high[,year_idx], n_lat_low[,year_idx],
         add=TRUE,
         col=rgb(0,0,1,0.4),
         errbar.col=rgb(0,0,1,0.4))
  
  abline(0,1,col="red",lwd=2,lty=2)
}