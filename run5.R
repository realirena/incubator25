library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores(logical= FALSE))
library(readxl)
# ---- 0) Prepare data ----
df <- read_xlsx("U:/Documents/repos/incubator25/data/modeldata.xlsx")

df0 <- df %>%
  select(state_name, year, deaths, yth, pov, propmis) 

stopifnot(nrow(df0) > 0)

# ---- 1) Numeric area index ----
df1 <- df0 %>%
  mutate(area = as.integer(factor(state_name)))
R <- length(unique(df1$area))

# ---- 2) Area-mean Poverty (center/scale) ----
poverty_area <- df1 %>%
  group_by(area) %>%
  summarise(Poverty_area = mean(pov, na.rm = TRUE), .groups = "drop") %>%
  arrange(area) %>%
  mutate(Poverty_c = scale(Poverty_area, center = TRUE, scale = TRUE)[,1])

# ---- 3) Center/scale covariates ----
df1 <- df1 %>%
  mutate(
    Youth_c   = scale(yth,  center = TRUE, scale = TRUE)[,1],
    Time_c    = scale(year, center = TRUE, scale = TRUE)[,1],
    Missing_c = scale(propmis, center = TRUE, scale = TRUE)[,1]
  )

# ---- 4) Stan data list ----
stan_data <- list(
  n       = nrow(df1),
  R       = R,
  area    = df1$area,
  z       = as.integer(df1$deaths),
  Youth   = df1$Youth_c,
  Poverty = poverty_area$Poverty_c,
  Time    = df1$Time_c,
  Missing = df1$Missing_c
)

# ---- 5) Fit the model ----
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit_time <- stan(
  file   = "U:/Documents/repos/incubator25/model/model5.stan",
  data   = stan_data,
  chains = 1,
  iter   = 1000,
  warmup = 500,
  cores  = 4,
  seed   = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  init_r = 0.1
)

# ---- 6) Inspect results ----
print(fit_time, pars = c("beta0","beta1","beta2","beta3","alpha","kappa2", "tau"), probs = c(0.05,0.5,0.95))

# ---- 7) Posterior traceplots ----
library(bayesplot)
posterior <- as.array(fit_time)
trace_params <- c("beta0","beta1","beta2","beta3","alpha","kappa2", "tau")
mcmc_trace(posterior, pars = trace_params, facet_args = list(ncol = 2))
mcmc_pairs(posterior, pars = trace_params, facet_args = list(ncol = 2))

rstan::traceplot(fit_time, pars=c("tau"))

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
