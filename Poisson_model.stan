// =============================================================
// Poisson detection with spatial ICAR smoothing on w
// w = softmax(f), f ~ ICAR(τ)
// YEAR TOTALS: LINEAR ON ORIGINAL SCALE (NO LOG)
// =============================================================
data {
  int<lower=1> K;                  // number of sites
  int<lower=1> N;                  // observations
  array[N] int<lower=0> y;         // counts
  array[N] int<lower=1> location;  // obs -> site
  int<lower=1> Y;
  array[N] int<lower=1,upper=Y> year;

  // Detection design (no intercept in X)
  int<lower=0> P;
  matrix[N, P] X;

  // Year totals design (centered & scaled year)
  int<lower=0> Q;
  matrix[Y, Q] Z_year;
  
  // Proportion of land available
  array[K] real<lower=0, upper=1> prop_land;

  // Priors for YEAR TOTALS on ORIGINAL scale
  real alpha_N_mean;                 // mean(N_total_year) prior (original scale)
  real<lower=0> alpha_N_sd;          // sd(N_total_year) prior (original scale)
  vector[Q] gamma_mean;              // slope(s) prior mean (original scale)
  vector<lower=0>[Q] gamma_sd;       // slope(s) prior sd   (original scale)
  real<lower=0> sigma_N_prior;       // prior scale for residual sd on original scale

  // ---- ICAR structure ----
  int<lower=1> N_edges;              // number of adjacency edges
  array[N_edges, 2] int<lower=1> edges; // (i, j) neighbors (i<j)
  array[K] int<lower=1> N_neighbors; // provided for completeness

  // ICAR prior scale
  real<lower=0> tau_prior;

  // Priors for detection
  real beta0_mean;
  real<lower=0> beta0_sd;
  vector[P] beta_mean;
  vector<lower=0>[P] beta_sd;
}
parameters {
  // ICAR latent field (non-centered)
  vector[K] f_raw;              // unscaled ICAR
  real<lower=0> tau;            // marginal SD scale

  // Detection parameters
  real beta0;
  vector[P] beta;

  // Year totals (ORIGINAL SCALE)
  real alpha_N;                 // intercept on original scale
  vector[Q] gamma;              // slope(s) on original scale (per SD of year)
  real<lower=0> sigma_N;        // residual SD on original scale
  vector<lower=0>[Y] N_total_year;  // yearly totals (constrained positive)
}
transformed parameters {
  vector[K] f;        // spatial field
  vector[K] w_raw;
  simplex[K] w;       // site proportions

  f = tau * f_raw;    // scaled ICAR
  
  for (k in 1:K) {
    w_raw[k] = exp(f[k]) * prop_land[k];
  }

  w = w_raw / sum(w_raw);     // maps to simplex (sum w = 1)
}
model {
  // ---- ICAR prior (Besag) ----
  for (e in 1:N_edges) {
    int i = edges[e, 1];
    int j = edges[e, 2];
    target += -0.5 * square(f_raw[i] - f_raw[j]); // Normal(0,1) on differences
  }
  f_raw ~ normal(0, 0.1);       // soft sum-to-zero anchoring
  tau ~ normal(0, tau_prior);   // half-normal via <lower=0>

  // ---- Detection priors ----
  beta0 ~ normal(beta0_mean, beta0_sd);
  if (P > 0) beta ~ normal(beta_mean, beta_sd);

  // ---- Year totals (ORIGINAL SCALE, LINEAR in Z_year) ----
  alpha_N ~ normal(alpha_N_mean, alpha_N_sd);
  if (Q > 0) gamma ~ normal(gamma_mean, gamma_sd);
  sigma_N ~ normal(0, sigma_N_prior); // half-normal via <lower=0>

  for (yid in 1:Y) {
    real mu_lin = alpha_N + (Q > 0 ? Z_year[yid] * gamma : 0);
    N_total_year[yid] ~ normal(mu_lin, sigma_N);
  }

  // ---- Likelihood ----
  for (i in 1:N) {
    int k = location[i];
    real eta = beta0 + (P > 0 ? dot_product(X[i], beta) : 0);
    real q = exp(eta);  // detection multiplier
    real lambda = (N_total_year[year[i]] * w[k]) * q;
    y[i] ~ poisson(lambda);
  }
}

generated quantities {
    // 1. Posterior n_latent for each site-year
    //    n_latent[k,y] = N_total_year[y] * w[k]
    array[K, Y] real n_latent;

    for (k in 1:K) {
        for (yid in 1:Y) {
            n_latent[k, yid] = N_total_year[yid] * w[k];
        }
    }

    // 2. Posterior detection probability for each observation
    //    p_det[i] = inv_logit(beta0 + X[i] * beta)
    vector[N] p_det;

    for (i in 1:N) {
        real eta = beta0;
        if (P > 0) {
            eta += dot_product(X[i], beta);
        }
        p_det[i] = inv_logit(eta);
    }
    
    // likelihood
    
    vector[N] log_lik;
    for (i in 1:N) {
      int k = location[i];
      real eta = beta0 + dot_product(X[i], beta);
      real lambda = N_total_year[year[i]] * w[k] * exp(eta);
      log_lik[i] = poisson_lpmf(y[i] | lambda);
    }

}
