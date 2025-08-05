library(rstan)
library(dplyr)

set.seed(123)

# Dimensions
A <- 50       
K <- 5       
P <- 3       

# Simulate covariates
X <- matrix(rnorm(A * P), nrow = A, ncol = P)

# regression parameters
beta_0_true <- 0.5
beta_true <- c(0.3, -0.2, 0.1)

# Expected counts
E <- runif(A, 5, 15)

# Stick-breaking parameters
alpha_true <- rep(2, K)
nu_true <- rep(2, K)

Z_true <- rbeta(K, alpha_true, nu_true)

# gamma
gamma_true <- numeric(K)
gamma_true[1] <- Z_true[1]
for(k in 2:K) {
  gamma_true[k] <- Z_true[k] * prod(1 - Z_true[1:(k-1)])
}

# Split coding scheme
cluster_ids <- rep(1:K, length.out = A)
H <- matrix(0, nrow = K, ncol = A)
for(i in 1:A) {
  H[, i] <- as.numeric(1:K <= cluster_ids[i])
}

# reporting probabilities
eps_true <- numeric(A)
for(i in 1:A) {
  eps_true[i] <- 1 - sum(H[, i] * gamma_true)
}

# relative risk
eta_true <- beta_0_true + X %*% beta_true
theta_true <- exp(eta_true)

# lambda
lambda_true <- E * theta_true * eps_true

# Simulate observed counts 
y <- rpois(A, lambda_true)
stan_data <- list(
  A = A,
  K = K,
  P = P,
  y = y,
  E = E,
  X = X,
  H = H
)

# fit
fit <- stan(
  file = "model.stan",
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95)
)
