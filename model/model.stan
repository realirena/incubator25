data {
  int<lower=1> A;                     // Number of areas
  int<lower=1> K;                     // Number of clusters
  int<lower=1> J;                     // timepoints
  int<lower=1> b_s; // number of spline regression coeffs
  int<lower=1> ids[J];              // regions over time 
  int y[J];                  // Observed counts
  vector<lower=0>[J] E;               // Expected counts 
  int<lower=1> P;                     // Number of covariates
  matrix[A, P] X;                     // Covariate matrix (P covariates) - can be time varying or not
  matrix[J, b_s] time;               // time variable (as b-splines)
  matrix[K, A] H;                     // Split-coded matrix H[k, i] = h_ki
}

parameters {
  real beta_0; //intercept 
  vector[P] beta; //regression coefficients 
  vector[b_s] b_k; // spline coefficients
  vector<lower=0>[K] alpha; //for probability of reporting
  vector<lower=0>[K] nu; //for probability of reporting
  vector<lower=0, upper=1>[K] Z; //prob death is reported or not
}

transformed parameters {
  vector<lower=0, upper=1>[K] gamma; //underreporting probability
  vector<lower=0, upper=1>[A] eps; //reporting probability
  vector<lower=0>[J] lambda; //mean rate in each area 
  vector[J] eta; //relative risk(?)


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

  }
  for(j in 1:J){
    eta[j] = beta_0 + dot_product(X[ids[j]], beta) + dot_product(time[j],b_k); 
    lambda[j] = E[j] * exp(eta[j]) * eps[ids[j]]; 
  }

}

model {
  // Priors
  beta_0 ~ normal(0, 1);
  beta ~ normal(0, 2);
  b_k ~ normal(0,1); 
  alpha ~ gamma(2, 0.5); 
  nu ~ gamma(2, 0.5); 

  for (k in 1:K)
    Z[k] ~ beta(alpha[k], nu[k]);  

  // Likelihood
  for (j in 1:J)
    y[j] ~ poisson(lambda[j]);
}
