data {
  // --- Sizes ---
  int<lower=1> L;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1> A;
  int<lower=1> K;
  int<lower=1> p;

  // All rows (observed + missing)
  int<lower=1> N_all;

  // Observed subset
  int<lower=1> Nobs;
  array[Nobs] int<lower=1, upper=N_all> obs_idx;
  array[Nobs] int<lower=0> D_obs;

  // Missing subset (convenience)
  int<lower=0> Nmis;
  array[Nmis] int<lower=1, upper=N_all> mis_idx;

  // Indices/covariates for ALL rows
  array[N_all] int<lower=1, upper=L> ii_l_all;
  array[N_all] int<lower=1, upper=S> ii_s_all;
  array[N_all] int<lower=1, upper=T> ii_t_all;
  array[N_all] int<lower=1, upper=A> ii_a_all;
  vector[N_all] log_P_all;
  matrix[N_all, p] Xmat_all;   // linear covariates only
  vector[N_all] Wrow_all;      // linear W only

  // Age basis
  matrix[A, K] B_age;

  // Intercept prior for deaths
  real mu_mean;
  real<lower=0> mu_sd;

  // Beta(a_alpha,b_alpha) prior for p0 via meta-analysis (on p0)
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;

  // BYM2 per country (component 1)
  int<lower=1> L1;
  int<lower=0> E1;
  array[E1] int<lower=1, upper=L1> node1_c1;
  array[E1] int<lower=1, upper=L1> node2_c1;
  array[L1] int<lower=1, upper=L> loc_ids_c1;
  real<lower=0> scale_icar_c1;

  // BYM2 per country (component 2)
  int<lower=1> L2;
  int<lower=0> E2;
  array[E2] int<lower=1, upper=L2> node1_c2;
  array[E2] int<lower=1, upper=L2> node2_c2;
  array[L2] int<lower=1, upper=L> loc_ids_c2;
  real<lower=0> scale_icar_c2;
}

parameters {
  // Global / fixed effects
  real mu;
  real alpha;                 // logit for p0
  vector[p] beta;             // linear X covariates
  real gamma;                 // linear W effect

  // SDs
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;  // hierarchical scale for gamma

  // Sex/time random effects + SDs
  vector[S] delta;
  real<lower=0> sigma_s;

  vector[T] zeta;
  real<lower=0> sigma_t;

  // Spatial BYM2
  real<lower=0> sigma_c;
  real<lower=0, upper=1> rho;

  vector[L1] theta_c1;
  vector[L1] phi_c1;

  vector[L2] theta_c2;
  vector[L2] phi_c2;

  // Age spline coefficients (RW2)
  vector[K] lambda;
  real<lower=0> sigma_rw2;
}

transformed parameters {
  // Average reporting rate when covariates are at reference values
  real<lower=0, upper=1> p0 = inv_logit(alpha);

  vector[L1] bym2_c1 = sqrt(1 - rho) * theta_c1
                     + sqrt(rho / scale_icar_c1) * phi_c1;

  vector[L2] bym2_c2 = sqrt(1 - rho) * theta_c2
                     + sqrt(rho / scale_icar_c2) * phi_c2;

  vector[L] phi_global = rep_vector(0, L);

  for (i in 1:L1)
    phi_global[loc_ids_c1[i]] = sigma_c * bym2_c1[i];

  for (i in 1:L2)
    phi_global[loc_ids_c2[i]] = sigma_c * bym2_c2[i];
}

model {
  // ---- Priors ----
  mu ~ normal(mu_mean, mu_sd);

  // Beta(a_alpha, b_alpha) prior on p0 = inv_logit(alpha) with Jacobian
  p0 ~ beta(a_alpha, b_alpha);
  target += alpha - 2 * log1p_exp(alpha);  // Jacobian term

  sigma_s     ~ student_t(5, 0, 1);
  sigma_t     ~ student_t(5, 0, 1);
  sigma_c     ~ student_t(5, 0, 1);
  sigma_beta  ~ student_t(5, 0, 1);
  sigma_rw2   ~ student_t(5, 0, 1);
  sigma_gamma ~ student_t(5, 0, 1);

  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);

  delta ~ normal(0, sigma_s);
  zeta  ~ normal(0, sigma_t);

  // Soft sum-to-zero for fixed effects
  sum(delta) ~ normal(0, 0.001 * S);
  sum(zeta)  ~ normal(0, 0.001 * T);

  // BYM2 hyper prior and components
  rho ~ beta(0.5, 0.5);

  theta_c1 ~ normal(0, 1);
  theta_c2 ~ normal(0, 1);

  // ICAR priors (up to proportionality) + soft centering
  if (E1 > 0)
    target += -0.5 * dot_self(phi_c1[node1_c1] - phi_c1[node2_c1]);
  sum(phi_c1) ~ normal(0, 0.001 * L1);

  if (E2 > 0)
    target += -0.5 * dot_self(phi_c2[node1_c2] - phi_c2[node2_c2]);
  sum(phi_c2) ~ normal(0, 0.001 * L2);

  // RW2(lambda)
  if (K >= 1) lambda[1] ~ normal(0, 0.5);
  if (K >= 2) lambda[2] ~ normal(0, 0.5);
  if (K >= 3) {
    for (k in 3:K)
      target += normal_lpdf(lambda[k] - 2 * lambda[k - 1] + lambda[k - 2] | 0, sigma_rw2);
  }

  // ---- Likelihood on observed rows only ----
  {
    // Reporting model: alpha + gamma * W
    vector[N_all] eta_pi_all  = alpha + gamma * Wrow_all;
    vector[N_all] log_pi_all  = -log1p_exp(-eta_pi_all);

    vector[A] f_age = B_age * lambda;

    vector[N_all] eta_true_all = log_P_all
                               + mu
                               + Xmat_all * beta
                               + delta[ii_s_all]
                               + zeta[ii_t_all]
                               + phi_global[ii_l_all]
                               + f_age[ii_a_all];

    target += poisson_log_lpmf(D_obs |
                eta_true_all[obs_idx] + log_pi_all[obs_idx]);
  }
}

generated quantities {
  // Posterior means (already had these)
  vector[N_all] pi_row_all;       // reporting probability per row
  vector[N_all] mu_true_all;      // expected true deaths per row
  vector[N_all] mu_reported_all;  // expected reported deaths per row

  // NEW: posterior predictive integer counts for ALL rows
  int<lower=0> y_true_all[N_all]; // predictive "true" deaths
  int<lower=0> y_rep_all[N_all];  // predictive reported deaths

  // Age effect reused here
  vector[A] f_age = B_age * lambda;

  {
    vector[N_all] eta_pi_all  = alpha + gamma * Wrow_all;
    vector[N_all] log_pi_all  = -log1p_exp(-eta_pi_all);

    vector[N_all] eta_true_all = log_P_all
                               + mu
                               + Xmat_all * beta
                               + delta[ii_s_all]
                               + zeta[ii_t_all]
                               + phi_global[ii_l_all]
                               + f_age[ii_a_all];

    for (n in 1:N_all) {
      // Means
      pi_row_all[n]      = inv_logit(eta_pi_all[n]);
      mu_true_all[n]     = exp(eta_true_all[n]);
      mu_reported_all[n] = exp(eta_true_all[n] + log_pi_all[n]);

      // Posterior predictive integer counts (for ALL rows)
      y_true_all[n] = poisson_rng(mu_true_all[n]);
      y_rep_all[n]  = poisson_rng(mu_reported_all[n]);
    }
  }
}

