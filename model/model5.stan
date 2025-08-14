data {
  int<lower=1> n;                  // number of observations
  int<lower=1> R;                  // number of areas
  int<lower=1, upper=R> area[n];   // area index
  int<lower=0> z[n];               // observed counts
  vector[n] Youth;                 // observation-level
  vector[R] Poverty;               // area-level mean Poverty
  vector[n] Time;                  // observation-level
  vector[n] Missing;               // observation-level
}
parameters {
  // Fixed effects
  real beta0;       // global intercept
  real beta1;       // Youth
  real beta2;       // Poverty (area-level)
  real beta3;       // Time
  //random effects
  real<lower=0> tau;
  real u[R]; 
  // Reporting probability
  real alpha;
  real kappa2;
}
transformed parameters {
  vector[n] log_lambda;
  vector[n] lambda;
  vector[n] pi;
  for (i in 1:n) {
    log_lambda[i] = beta0 + u[area[i]] + beta1 * Youth[i] + beta2 * Poverty[area[i]] + beta3 * Time[i];
    lambda[i] = exp(log_lambda[i]);
    pi[i] = inv_logit(alpha + kappa2 * Missing[i]);
  }
}
model {
  // Priors
  beta0 ~ normal(0, 1);
  u ~ normal(0, tau); 
  tau ~ cauchy(0, 2.5);
  beta1 ~ normal(0, 0.5);
  beta2 ~ normal(0, 0.5);
  beta3 ~ normal(0, 0.5);
  alpha  ~ normal(0.85, 0.5);
  kappa2 ~ normal(0, 1);
  // Likelihood
  z ~ poisson(pi .* lambda);
}






