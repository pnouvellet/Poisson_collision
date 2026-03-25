make_stan_data_formula_sumcoded_icar <- function(
    df_long,
    detection_formula,
    year_formula,
    
    # Priors
    beta0_mean = 0, beta0_sd = 1,
    beta_mean = NULL, beta_sd = 0.5,
    N_overall_mean, N_overall_sd,
    gamma_mean = NULL, gamma_sd = 200,
    sigma_N_prior = 50,
    
    # ICAR prior scale
    tau_prior = 1.0,
    verbose = TRUE
){
  if (verbose) cat("\n=== PREPARING DATA (ICAR spatial smoothing) ===\n")
  
  ###############################################################
  # BASIC INDEXING
  ###############################################################
  df_long$location <- as.integer(as.factor(df_long$location))
  df_long$year     <- as.integer(as.factor(df_long$year))
  
  K <- length(unique(df_long$location))
  Y <- length(unique(df_long$year))
  N <- nrow(df_long)
  
  years_unique <- sort(unique(df_long$year))
  df_year <- data.frame(year = years_unique)
  
  ###############################################################
  # DETECTION DESIGN MATRIX
  ###############################################################
  as_rhs <- function(f) {
    f <- if (is.character(f)) as.formula(f) else f
    if (length(f)==3) as.formula(paste("~", deparse(f[[3]]))) else f
  }
  drop_intercept <- function(M) {
    if ("(Intercept)" %in% colnames(M))
      M[, setdiff(colnames(M), "(Intercept)"), drop=FALSE]
    else
      M
  }
  
  f_det <- as_rhs(detection_formula)
  mf_det <- model.frame(f_det, df_long, na.action = na.fail)
  
  is_fac_det <- vapply(mf_det, is.factor,  logical(1))
  is_num_det <- vapply(mf_det, is.numeric, logical(1))
  
  # Standardize numeric predictors
  for (nm in names(mf_det)[is_num_det]) {
    sdv <- sd(mf_det[[nm]])
    if (sdv <= 0) sdv <- 1
    mf_det[[nm]] <- (mf_det[[nm]] - mean(mf_det[[nm]])) / sdv
  }
  
  det_fac_names <- names(mf_det)[is_fac_det]
  det_contr <- if (length(det_fac_names)) {
    setNames(
      lapply(det_fac_names, function(nm)
        contr.sum(nlevels(mf_det[[nm]]))
      ),
      det_fac_names
    )
  } else NULL
  
  X_full <- model.matrix(attr(mf_det, "terms"), mf_det,
                         contrasts.arg = det_contr)
  X <- drop_intercept(X_full)
  P <- ncol(X)
  if (P == 0) X <- matrix(0, nrow = N, ncol = 0)
  
  ###############################################################
  # YEAR TOTALS DESIGN MATRIX
  ###############################################################
  f_year <- as_rhs(year_formula)
  mf_year <- model.frame(f_year, df_year, na.action = na.fail)
  
  is_fac_year <- vapply(mf_year, is.factor, logical(1))
  is_num_year <- vapply(mf_year, is.numeric, logical(1))
  
  for (nm in names(mf_year)[is_num_year]) {
    sdv <- sd(mf_year[[nm]])
    if (sdv <= 0) sdv <- 1
    mf_year[[nm]] <- (mf_year[[nm]] - mean(mf_year[[nm]])) / sdv
  }
  
  year_fac_names <- names(mf_year)[is_fac_year]
  year_contr <- if (length(year_fac_names)) {
    setNames(
      lapply(year_fac_names, function(nm)
        contr.sum(nlevels(mf_year[[nm]]))
      ),
      year_fac_names
    )
  } else NULL
  
  Z_full <- model.matrix(attr(mf_year, "terms"), mf_year,
                         contrasts.arg = year_contr)
  Z_year <- drop_intercept(Z_full)
  Q <- ncol(Z_year)
  if (Q == 0) Z_year <- matrix(0, Y, 0)
  
  ###############################################################
  # PRIORS
  ###############################################################
  if (is.null(beta_mean))
    beta_mean <- rep(0, P)  # uninformative mean
  if (length(beta_sd) == 1)
    beta_sd <- rep(beta_sd, P)
  
  if (is.null(gamma_mean))
    gamma_mean <- rep(0, Q)
  if (length(gamma_sd) == 1)
    gamma_sd <- rep(gamma_sd, Q)
  
  alpha_N_mean <- N_overall_mean
  alpha_N_sd   <- N_overall_sd
  
  ###############################################################
  # ICAR ADJACENCY (ROOK GRID)
  ###############################################################
  if (!("x_coord" %in% names(df_long)) ||
      !("y_coord" %in% names(df_long)))
    stop("Need x_coord and y_coord in df_long for ICAR adjacency.")
  
  coords_by_loc <- df_long %>%
    group_by(location) %>%
    summarise(x = x_coord[1], y = y_coord[1], .groups = "drop") %>%
    arrange(location)
  
  xs <- sort(unique(coords_by_loc$x))
  ys <- sort(unique(coords_by_loc$y))
  nx <- length(xs)
  ny <- length(ys)
  
  edges <- list()
  N_neighbors <- integer(K)
  loc_mat <- matrix(1:K, nrow = ny, ncol = nx, byrow = TRUE)
  
  for (iy in 1:ny) {
    for (ix in 1:nx) {
      k <- loc_mat[iy, ix]
      neigh <- c()
      if (ix > 1)  neigh <- c(neigh, loc_mat[iy, ix-1])
      if (ix < nx) neigh <- c(neigh, loc_mat[iy, ix+1])
      if (iy > 1)  neigh <- c(neigh, loc_mat[iy-1, ix])
      if (iy < ny) neigh <- c(neigh, loc_mat[iy+1, ix])
      N_neighbors[k] <- length(neigh)
      for (j in neigh)
        if (k < j) edges[[length(edges)+1]] <- c(k, j)
    }
  }
  
  edges_mat <- do.call(rbind, edges)
  N_edges <- nrow(edges_mat)
  
  ###############################################################
  # BUILD FINAL DATA LIST
  ###############################################################
  dat <- list(
    K = K, N = N,
    y = as.integer(df_long$y),
    location = as.integer(df_long$location),
    
    prop_land = tapply(df_long$prop_land, df_long$location, mean),
    
    Y = Y,
    year = as.integer(df_long$year),
    
    P = P, X = X,
    Q = Q, Z_year = Z_year,
    
    alpha_N_mean = alpha_N_mean,
    alpha_N_sd   = alpha_N_sd,
    gamma_mean   = as.vector(gamma_mean),
    gamma_sd     = as.vector(gamma_sd),
    sigma_N_prior = sigma_N_prior,
    
    # ICAR adjacency info
    N_edges = N_edges,
    edges = edges_mat,
    N_neighbors = N_neighbors,
    tau_prior = tau_prior,
    
    # Detection priors
    beta0_mean = beta0_mean,
    beta0_sd   = beta0_sd,
    beta_mean  = as.vector(beta_mean),
    beta_sd    = as.vector(beta_sd)
  )
  
  list(dat = dat)
}
