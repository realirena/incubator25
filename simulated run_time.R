rm(list=ls())
library(rstan)
library(dplyr)
library(splines)
set.seed(123)

# Dimensions
A <- 50       
K <- 5       
P <- 3       

# Simulate covariates
X <- matrix(rnorm(A * P), nrow = A, ncol = P)

## generate time variable for mortality 

## simulate time
Ti = sample(3:12, A, replace=T) ##no visits per subject
J = sum(Ti)
## placeholder for the simulated time to FMPs
time = unlist(lapply(seq_along(Ti), function(i) seq(1:Ti[i])))
ids <- unlist(lapply(seq_along(Ti), function(i) rep(i, Ti[i])))

time_bs = bs(time, knots =c(5, 10), degree=1, intercept = FALSE)
b_k = c(-0.5, 0.6, 1.2)
b_s = ncol(time_bs)
mu_time <-time_bs%*%b_k

time_df <- data.frame(ids, time_bs, mu_time)
plot(time_df[time_df$ids==1,]$mu_time)


# regression parameters
beta_0_true <- 0.5
beta_true <- c(0.3, -0.2, 0.1)

# Expected counts
E <- runif(J, 5, 15)

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
eta_true <- drop(t(sapply(seq_along(ids), function(i){beta_0_true + X[ids[i],]%*% beta_true + mu_time[i]})))

theta_true <- exp(eta_true)

# lambda
lambda_true <-  drop(t(sapply(seq_along(ids), function(i){E[i] * theta_true[i] * eps_true[ids[i]]})))

# Simulate observed counts 
y <- rpois(J, lambda_true)


stan_data <- list(
  A = A,
  J = J,
  b_s = b_s,
  K = K,
  P = P,
  y = y,
  E = E,
  X = X,
  ids = ids,
  time = time_bs,
  H = H
)

# fit
fit <- stan(
  file = "U:/Documents/repos/incubator25/model/model.stan",
  data = stan_data,
  iter = 200,
  warmup = 100,
  chains = 1,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

traceplot(fit, pars=c("alpha", "nu"))
traceplot(fit, pars=c("beta_0", "beta", "b_k"))

plot(fit, pars=c("beta_0", "beta", "b_k"))
plot(fit, pars=c("alpha", "nu"))

summary(fit, pars=c("alpha", "nu"))
summary(fit, pars=c("beta_0", "beta"))$summary


