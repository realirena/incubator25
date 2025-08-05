data {
  int<lower=1> A;                     // Number of areas
  int<lower=1> K;                     // Number of clusters
  int<lower=0> y[A];                  // Observed counts
  vector<lower=0>[A] E;               // Expected counts
  matrix[A, P] X;                     // Covariate matrix (P covariates)
  matrix[K, A] H;                     // Split-coded matrix H[k, i] = h_ki
}

parameters {
  real beta_0; //intercept 
  vector[P] beta; //regression coefficients 
  vector<lower=0>[K] alpha; //for probability of reporting
  vector<lower=0>[K] nu; //for probability of reporting
  vector<lower=0, upper=1>[K] Z; //if death is reported or not
}

transformed parameters {
  vector<lower=0, upper=1>[K] gamma; //underreporting probability
  vector<lower=0, upper=1>[A] eps; //reporting probability
  vector<lower=0>[A] lambda; //mean rate in each area 
  vector[A] eta; //relative risk(?)

  gamma[1] = Z[1];
  for (k in 2:K) {
    real prod = 1.0;
    for (j in 1:(k - 1)) {
      prod *= (1.0 - Z[j]);
    }
    gamma[k] = Z[k] * prod;
  }

  for (i in 1:A) {
    eps[i] = 1.0 - dot_product(H[, i], gamma);
    eta[i] = beta_0 + dot_product(X[i], beta);
    lambda[i] = E[i] * exp(eta[i]) * eps[i];
  }
}

model {
  // Priors
  beta_0 ~ normal(0, 5);
  beta ~ normal(0, 2);

  alpha ~ gamma(2, 0.5);
  nu ~ gamma(2, 0.5);

  for (k in 1:K)
    Z[k] ~ beta(alpha[k], nu[k]);

  // Likelihood
  for (i in 1:A)
    y[i] ~ poisson(lambda[i]);
}
