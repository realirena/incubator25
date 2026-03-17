data{
  int<lower=1> N; // observations
  int<lower=1> A; // age
  int<lower=1> T; // year
  int<lower=1> L; // location
  array[N] int<lower=1,upper=A> an; // age index
  array[N] int<lower=1,upper=T> tn; // year index
  array[N] int<lower=1,upper=L> ln; // location
  
  int<lower=1> Nobs; // observed subset
  array[Nobs] int<lower=1,upper=N> obs; // observation index
  array[Nobs] int<lower=0>y; // observed homicide counts
  
  vector<lower=0>[N] E;      // Exposure 
  matrix[N,2] X; // covariates
  int<lower=0> K;
  matrix[A,K] B; // spline basis
  
  // neighbour for region effect
  int<lower=0> N_edges;
  array[N_edges] int<lower=1,upper=L> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1,upper=L> node2; // and node1[i] < node2[i]
  real<lower=0> scaling_factor; // scale the variances of spatial effect
}
transformed data{
  vector[N] log_E = log(E); //log of exposure
}
parameters{
  real alpha0;  // intercept for mortality
  vector[K] lambda; // RW2 spline coefficient
  vector[T] kappa; // year effect
  real drift; // drift for RW1 kappa
  
  real<lower=0> sigma_lambda;  // variances
  real<lower=0> sigma_kappa;
  real<lower=0> sigma_eps;
  real<lower=0> sigma_BYM2;
  
  real<lower=0, upper=1> rho; // proportion of structured vs unstructure variance
  vector[L] v;  // non-spatial effect
  vector[L] u;  // spatial effect
  
  vector[N] eps;     // overdispersion 
  
  real beta0; // intercept for  
  vector[2] beta; // beta coefficient for the covariates
}
transformed parameters{
  vector[L] convolved_re;  // total BYM2 effect
  vector[L] phi;    // total regional effect
  // variance of each component should be ~1
  convolved_re = sqrt(1-rho)*v+sqrt(rho/scaling_factor)*u; // Riebler et. al (2017)
  phi = convolved_re*sigma_BYM2; // noncentred parameterisation
  
  vector[A] alpha = B*lambda; // age effect
}
model{
  sigma_lambda ~ student_t(5,0,1);
  sigma_kappa ~ student_t(5,0,1);
  sigma_eps ~ student_t(5,0,1);
  sigma_BYM2 ~ student_t(5,0,1);
  
  drift ~ normal(0,2);
  
  alpha0 ~ normal(-7,1); // mild informative prior (male)
  //alpha0 ~ normal(-8,1); // mild informative prior (female)
  eps ~ normal(0,1);
  rho ~ beta(0.5,0.5);
  
  target += beta_lpdf(inv_logit(beta0)|29,2); // informative beta prior (male)
  //target += beta_lpdf(inv_logit(beta0)|18,1); // informative beta prior (female)
  target += log_inv_logit(beta0)+log1m_inv_logit(beta0); // Jacobian adjustment 
  beta ~ normal(0,1);
  
  // RW2 spline coefficient
  lambda[1:2] ~ normal(0, sigma_lambda);
  lambda[3:K] ~ normal(2*lambda[2:(K-1)]-lambda[1:(K-2)], sigma_lambda);
  
  // RW1 temporal trend
  kappa[1] ~ normal(drift, sigma_kappa);
  kappa[2:T] ~ normal(drift+kappa[1:(T-1)], sigma_kappa);
  sum(kappa) ~ normal(0, 0.001*T);
  
  // location effect (BYM2 model)
  target += -0.5 * dot_self(u[node1]-u[node2]);
  sum(u) ~ normal(0, 0.001*L);
  v ~ normal(0,1);
  
  vector[N] mu = alpha0 + alpha[an] + kappa[tn] + phi[ln];
  vector[N] epsilon = log_inv_logit(beta0 + X*beta);
  y ~ poisson_log(log_E[obs] + mu[obs] + sigma_eps*eps[obs] + epsilon[obs]);
}
generated quantities{
  array[N] int y_true;
  array[N] int y_rep;
  vector[A] alpha_gq = B*lambda;
  vector[L] convolved_re_gq = sqrt(1-rho)*v+sqrt(rho/scaling_factor)*u;
  vector[L] phi_gq = sigma_BYM2*convolved_re_gq;
  vector[N] mu_gq = alpha0+alpha_gq[an]+kappa[tn]+phi_gq[ln];
  vector[N] rep_prob = inv_logit(beta0+X*beta);
  vector[N] lambda_true = E .* exp(mu_gq+sigma_eps*eps);
  vector[N] lambda_obs = lambda_true .* rep_prob;
  for(n in 1:N){
    y_true[n] = poisson_rng(lambda_true[n]);
    y_rep[n] = poisson_rng(lambda_obs[n]);
  }
}

