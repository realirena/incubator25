# ==================== Packages ====================
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(stringi)
library(rstan)
library(bayesplot)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ==================== Helpers ====================
fix_mojibake <- function(x) {
  x <- as.character(x)
  x2 <- suppressWarnings(iconv(x, "UTF-8","latin1"))
  x3 <- suppressWarnings(iconv(x2, "latin1", "UTF-8"))
  x3[is.na(x3)] <- x[is.na(x3)]
  x3
}
norm_dept <- function(x) {
  x %>%
    fix_mojibake() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    tolower() %>%
    str_squish() %>%
    str_replace_all(c(
      "^n\\.?\\s*de\\s*santander$" = "norte de santander",
      "^valle$"                    = "valle del cauca",
      "^bogota d\\.?c\\.?$"        = "bogota",
      "^san andres y providencia$" = "san andres"
    ))
}
scale_within <- function(x, grp) {
  # z-score within grouping factor 'grp'
  unlist(by(x, grp, function(v) as.numeric(scale(v))), use.names = FALSE)
}

# ==================== Build joint dataset (Colombia + Mexico) ====================
countries_use <- c("colombia", "mexico")

admin0 <- df %>%
  filter(tolower(country) %in% countries_use) %>%
  rename(
    Z = deaths,
    mis = propmis,
    unem = unemployment,
    yth = yth,
    pop = pop,
    area_name = state_name,
    area_code = state
  ) %>%
  mutate(
    country_key   = tolower(country),
    area_name_key = norm_dept(area_name),
    area_code     = as.integer(area_code)
  ) %>%
  # keep rows with essential fields
  filter(!is.na(Z), !is.na(pop), !is.na(mis))

# ---- area-level means per (country, area) ----
area_means <- admin0 %>%
  group_by(country_key, area_code, area_name_key) %>%
  summarize(
    hdi_i   = mean(hdi,  na.rm = TRUE),
    pov_i   = mean(pov,  na.rm = TRUE),
    unem_i  = mean(unem, na.rm = TRUE),
    indig_i = mean(indig, na.rm = TRUE),
    .groups = "drop"
  )

# Fill NAs with country-wise means, then z-score within country
country_fill <- area_means %>%
  group_by(country_key) %>%
  summarize(
    hdi_fill   = mean(hdi_i,  na.rm = TRUE),
    pov_fill   = mean(pov_i,  na.rm = TRUE),
    unem_fill  = mean(unem_i, na.rm = TRUE),
    indig_fill = mean(indig_i,na.rm = TRUE),
    .groups = "drop"
  )

area_means2 <- area_means %>%
  left_join(country_fill, by = "country_key") %>%
  mutate(
    hdi_i   = ifelse(is.na(hdi_i),   hdi_fill,   hdi_i),
    pov_i   = ifelse(is.na(pov_i),   pov_fill,   pov_i),
    unem_i  = ifelse(is.na(unem_i),  unem_fill,  unem_i),
    indig_i = ifelse(is.na(indig_i), indig_fill, indig_i)
  ) %>%
  group_by(country_key) %>%
  mutate(
    hdi_i_z   = as.numeric(scale(hdi_i)),
    pov_i_z   = as.numeric(scale(pov_i)),
    unem_i_z  = as.numeric(scale(unem_i)),
    indig_i_z = as.numeric(scale(indig_i))
  ) %>%
  ungroup() %>%
  select(country_key, area_code, area_name_key, hdi_i_z, pov_i_z, unem_i_z, indig_i_z)

# ---- join back & standardize time-varying covariates within country ----
adminX <- admin0 %>%
  left_join(area_means2, by = c("country_key","area_code","area_name_key")) %>%
  group_by(country_key) %>%
  mutate(
    mis_z   = as.numeric(scale(mis)),     # within-country
    yth_z   = as.numeric(scale(yth)),     # within-country
    log_pop = log(pop/1e2)
  ) %>%
  ungroup()

# ---- indices: areas are unique by (country, area_code); years are global union ----
areas_tbl <- adminX %>%
  distinct(country_key, area_code, area_name_key) %>%
  arrange(country_key, area_code) %>%
  mutate(area_uid = paste(country_key, area_code, sep = "::"))

years_tbl <- adminX %>%
  distinct(year) %>%
  arrange(year)

adminX <- adminX %>%
  mutate(
    area_uid = paste(country_key, area_code, sep = "::"),
    area_id  = match(area_uid, areas_tbl$area_uid),
    year_id  = match(year, years_tbl$year)
  ) %>%
  arrange(area_id, year_id)

# ---- design matrices ----
X_alpha <- as.matrix(transmute(adminX, mis_z, hdi_i_z, indig_i_z))
X_beta  <- as.matrix(transmute(adminX, yth_z, pov_i_z, unem_i_z))

stan_data <- list(
  N = nrow(adminX),
  I = nrow(areas_tbl),
  T = nrow(years_tbl),
  P_alpha = ncol(X_alpha),
  P_beta  = ncol(X_beta),
  Z = as.integer(adminX$Z),
  log_pop = adminX$log_pop,
  X_alpha = X_alpha,
  X_beta  = X_beta,
  area_id = as.integer(adminX$area_id),
  year_id = as.integer(adminX$year_id)
)

cat(sprintf("Joint fit: %d rows, %d areas, %d years\n", stan_data$N, stan_data$I, stan_data$T))

# ==================== Compile & sample ====================
stan_file <- "model3.stan"  
sm <- rstan::stan_model(stan_file)

fit <- rstan::sampling(
  sm, data = stan_data,
  chains = 2, iter = 4000, warmup = 2000, cores = 2, seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit, pars = c("alpha0","beta0","eps","tau","phi_raw"))

# ==================== Diagnostics ====================
x <- as.array(fit)
pn <- dimnames(x)[[3]]

bayesplot::mcmc_trace(
  x,
  pars       = c("alpha0","beta0","eps","tau","phi_raw"),
  regex_pars = c("^alpha\\[", "^beta\\[")
)
ggsave("traceplot.pdf", height = 18, width = 20, units = "cm", dpi = 300)

alphas <- grep("^alpha\\[", pn, value = TRUE)
betas  <- grep("^beta\\[",  pn, value = TRUE)
keep   <- c("alpha0","beta0","eps","tau","phi_raw", head(alphas, 2), head(betas, 2))
bayesplot::mcmc_pairs(x, pars = keep)

# ==================== Posterior prediction ====================

library(dplyr)
library(tidyr)
library(tibble)

# 0) Map row n -> (area_id, year_id)
index_map <- adminX %>%
  transmute(n = row_number(), area_id, year_id)

# 1) Extract full draws arrays: [iterations x chains x N]
Zarr <- rstan::extract(fit, pars = "Z_rep",      permuted = FALSE)
Larr <- rstan::extract(fit, pars = "lambda_it",  permuted = FALSE)
Parr <- rstan::extract(fit, pars = "pi_it",      permuted = FALSE)

# Sanity: dimensions must match
stopifnot(identical(dim(Zarr), dim(Larr)), identical(dim(Zarr), dim(Parr)))

it <- dim(Zarr)[1]   # iterations per chain (post-warmup)
ch <- dim(Zarr)[2]   # number of chains
N  <- dim(Zarr)[3]   # rows in adminX

# 2) Build base index over all draws without relying on column names
post_long <- tibble(
  .iter  = rep(seq_len(it), times = ch * N),
  .chain = rep(rep(seq_len(ch), each = it), times = N),
  n      = rep(seq_len(N), each = it * ch)
) %>%
  mutate(
    .draw   = (.chain - 1L) * it + .iter,
    Z_rep   = as.vector(Zarr),       # vectorizes in the same (iter, chain, n) order
    lambda  = as.vector(Larr),
    pi      = as.vector(Parr)
  ) %>%
  left_join(index_map, by = "n") %>%
  # Area labels
  left_join(
    areas_tbl %>%
      mutate(area_id = row_number()) %>%
      left_join(
        admin0 %>% distinct(country_key, area_code, area_name),
        by = c("country_key","area_code")
      ) %>%
      mutate(state = coalesce(area_name, area_name_key)) %>%
      select(area_id, country_key, area_code, state),
    by = "area_id"
  ) %>%
  # Year labels
  left_join(
    years_tbl %>% mutate(year_id = row_number()) %>% select(year_id, year),
    by = "year_id"
  ) %>%
  mutate(Z_obs = adminX$Z[n]) %>%
  arrange(state, year, .draw) %>%
  relocate(state, year, .draw, .chain, .iter, country_key, area_code,
           Z_obs, Z_rep, lambda, pi)

# 3) (Optional) keep the same random 50 draws for everyone
# set.seed(123)
# keep_draws <- sample(unique(post_long$.draw), size = min(50, dplyr::n_distinct(post_long$.draw)))
# post_long_50 <- post_long %>% filter(.draw %in% keep_draws)

# 4) Save
write.csv(post_long, "posterior_estimates.csv", row.names = FALSE)
# If you made the 50-draw subset:
# write.csv(post_long_50, "posterior_estimates_50draws_with_lambda_pi.csv", row.names = FALSE)

# select 50 draws at random
set.seed(123)  # for reproducibility

# 1) pick the same random draws for everyone
all_draws <- sort(unique(postZ_long$.draw))
n_keep    <- min(50, length(all_draws))
keep_draws <- sample(all_draws, n_keep, replace = FALSE)

# 2) filter your long data to those draws
postZ_50_same <- postZ_long %>%
  dplyr::filter(.draw %in% keep_draws) %>%
  dplyr::arrange(country_key, state, year, .draw)
write.csv(postZ_50_same, "Posterior_draws_50.csv")
