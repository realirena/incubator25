# ---- 1) Extract posterior samples ----
posterior_df <- as.data.frame(fit_time)  # each row = draw

n_draws <- nrow(posterior_df)
n <- stan_data$n

# ---- 2) Storage for posterior predictive draws ----
y_rep <- matrix(NA, nrow = n_draws, ncol = n)

# ---- 3) Loop over posterior draws ----
for (s in 1:n_draws) {
  
  # Compute log_lambda for all observations
  log_lambda <- posterior_df$beta0[s] +
    posterior_df$beta1[s] * stan_data$Youth +
    posterior_df$beta2[s] * stan_data$Poverty[stan_data$area] +
    posterior_df$beta3[s] * stan_data$Time
  
  lambda <- exp(log_lambda)
  
  # Reporting probability
  pi <- 1 / (1 + exp(-(posterior_df$alpha[s] + posterior_df$kappa2[s] * stan_data$Missing)))
  
  # Posterior predictive Poisson draw
  y_rep[s, ] <- rpois(n, lambda * pi)
}

# ---- 4) Convert to data frame and export ----
y_rep_df <- as.data.frame(y_rep)
colnames(y_rep_df) <- paste0("obs_", 1:n)  

write.csv(y_rep_df, "posterior_predictive_y.csv", row.names = FALSE)

