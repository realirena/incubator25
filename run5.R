library(dplyr)
library(rstan)

# ---- 0) Prepare data ----
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
  file   = "model5.stan",
  data   = stan_data,
  chains = 4,
  iter   = 5000,
  warmup = 1000,
  cores  = 4,
  seed   = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  init_r = 0.1
)

# ---- 6) Inspect results ----
print(fit_time, pars = c("beta1","beta2","beta3","alpha","kappa2"), probs = c(0.05,0.5,0.95))

# ---- 7) Posterior traceplots ----
library(bayesplot)
posterior <- as.array(fit_time)
trace_params <- c("beta1","beta2","beta3","alpha","kappa2")
mcmc_trace(posterior, pars = trace_params, facet_args = list(ncol = 2))
mcmc_pairs(posterior, pars = trace_params, facet_args = list(ncol = 2))
