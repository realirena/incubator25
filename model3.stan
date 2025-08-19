// homicide_reporting_simple_AR1_tau.stan
data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> T;
  int<lower=1> P_alpha;
  int<lower=1> P_beta;

  int<lower=0> Z[N];
  vector[N] log_pop;
  matrix[N, P_alpha] X_alpha;
  matrix[N, P_beta]  X_beta;
  int<lower=1,upper=I> area_id[N];
  int<lower=1,upper=T> year_id[N];
}

parameters {
  // ---- fixed effects ----
  real alpha0;
  vector[P_alpha] alpha;
  real beta0;
  vector[P_beta] beta;

  // ---- intensity random effects (non-centered) ----
  vector[I] eta_raw;                 // standard normal
  real<lower=0> eps;                 // sd(eta)

  // ---- reporting AR(1): non-centered with marginal SD ----
  matrix[I, T] z;                    // i.i.d. N(0,1) innovations
  real<lower=0> tau;                 // marginal SD of theta
  real<lower=0, upper=1> phi_raw;    // map to phi in (-1,1)
}

transformed parameters {
  real phi = 2 * phi_raw - 1;        // AR(1) persistence in (-1,1)
  vector[I] eta = eta_raw * eps;
  vector[I] eta_c = eta - mean(eta); // center to separate from beta0

  // Build AR(1) process with marginal SD tau:
  matrix[I, T] theta;
  for (i in 1:I) {
    theta[i, 1] = tau * z[i, 1];
    for (t in 2:T) {
      theta[i, t] = phi * theta[i, t - 1] + tau * sqrt(1 - square(phi)) * z[i, t];
    }
  }

  vector[N] logit_pi;
  vector[N] log_lambda;
  vector[N] log_mu;
  for (n in 1:N) {
    logit_pi[n]   = alpha0 + X_alpha[n] * alpha + theta[ area_id[n], year_id[n] ];
    log_lambda[n] = log_pop[n] + beta0 + X_beta[n] * beta + eta_c[ area_id[n] ];
    log_mu[n]     = log_inv_logit(logit_pi[n]) + log_lambda[n]; // log(pi * lambda)
  }
}

model {
  // ---- priors ----
  alpha0 ~ normal(0.8, 0.08);

  alpha  ~ normal(0, 1);
  beta0  ~ normal(-2, 1);
  beta   ~ normal(0, 0.5);

  eps     ~ normal(0, 0.5);
  eta_raw ~ normal(0, 1);

  phi_raw      ~ beta(2, 2);
  tau          ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);

  // ---- likelihood ----
  Z ~ poisson_log(log_mu);
}

generated quantities {
  vector[N] pi_it;
  vector[N] lambda_it;
  int Z_rep[N];                      // posterior predictive for reported counts

  for (n in 1:N) {
    pi_it[n]     = inv_logit(logit_pi[n]);
    lambda_it[n] = exp(log_lambda[n]);
    Z_rep[n]     = poisson_log_rng(log_mu[n]);   // draw Z~Poisson(pi*lambda*pop)
  }
}
