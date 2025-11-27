## ---- packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  sf, dplyr, tidyr, stringi, rnaturalearth, rnaturalearthdata,
  spdep, Matrix, INLA, splines, rstan, igraph, purrr, tibble, metafor
)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## ---- inputs present? ----
stopifnot(exists("df"))
stopifnot(exists("rp"))
stopifnot(all(c("country","state","year","age","sex","deaths","pop",
                "propyth","pov","propmis_bothsex","propcomp_bothsex") %in% names(df)))
stopifnot(all(c("replicate","state","year","period3y","sample_num","prop_completeness") %in% names(rp)))

## ---- helpers ----
norm_txt <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    tolower() |>
    gsub("[^a-z0-9]+", " ", x = _) |>
    trimws()
}

manual_map <- c(
  # Colombia
  "valle"               = "valle del cauca",
  "san andres"          = "san andres",
  "bogota"              = "bogota",
  "narino"              = "narino",
  "choco"               = "choco",
  "cordoba"             = "cordoba",
  "guainia"             = "guainia",
  "vaupes"              = "vaupes",
  "norte de santander"  = "norte de santander",
  "la guajira"          = "la guajira",
  # Mexico
  "mexico city"         = "distrito federal",
  "mexico state"        = "mexico",
  "nuevo leon"          = "nuevo leon",
  "san luis potosi"     = "san luis potosi",
  "michoacan"           = "michoacan",
  "queretaro"           = "queretaro",
  "veracruz"            = "veracruz",
  "quintana roo"        = "quintana roo"
)

build_A_for_country <- function(country_name, user_states_vec) {
  adm <- rnaturalearth::ne_states(country = country_name, returnclass = "sf") |>
    dplyr::select(name, name_en, gn_name, postal, type, geometry) |>
    mutate(
      name_std    = norm_txt(coalesce(name, "")),
      name_en_std = norm_txt(coalesce(name_en, "")),
      gn_std      = norm_txt(coalesce(gn_name, "")),
      any_std     = paste(name_std, name_en_std, gn_std)
    )
  user_std <- norm_txt(user_states_vec)
  map <- manual_map
  
  match_idx <- integer(length(user_std))
  for (i in seq_along(user_std)) {
    u_std <- user_std[i]
    key <- if (!is.na(map[u_std])) map[u_std] else u_std
    cand <- which(
      vapply(adm$name_std,     \(nm) grepl(key, nm, fixed = TRUE), TRUE) |
        vapply(adm$name_en_std,  \(nm) grepl(key, nm, fixed = TRUE), TRUE) |
        vapply(adm$gn_std,       \(nm) grepl(key, nm, fixed = TRUE), TRUE) |
        vapply(adm$any_std,      \(nm) grepl(key, nm, fixed = TRUE), TRUE)
    )
    if (length(cand) == 0)
      stop(sprintf("No admin-1 match for '%s' in %s", user_states_vec[i], country_name))
    exact <- cand[adm$name_std[cand] == key]
    if (length(exact) >= 1) cand <- exact
    match_idx[i] <- cand[1]
  }
  
  shp <- adm[match_idx, ]
  shp$user_state <- user_states_vec
  row.names(shp) <- shp$user_state
  
  nb_q <- spdep::poly2nb(shp, row.names = shp$user_state, queen = TRUE, snap = 1e-8)
  A <- spdep::nb2mat(nb_q, style = "B", zero.policy = TRUE)
  deg <- rowSums(A)
  
  if (any(deg == 0)) {
    cent   <- sf::st_point_on_surface(shp$geometry)
    coords <- sf::st_coordinates(cent); rownames(coords) <- shp$user_state
    nb_knn1 <- spdep::knn2nb(spdep::knearneigh(coords, k = 1, longlat = FALSE))
    Ak <- spdep::nb2mat(nb_knn1, style = "B", zero.policy = TRUE)
    Ak <- (Ak + t(Ak) > 0) * 1L; diag(Ak) <- 0L
    iso <- which(deg == 0)
    A[iso, ] <- Ak[iso, ]; A[, iso] <- Ak[, iso]
    deg <- rowSums(A)
  }
  stopifnot(all(deg > 0))
  rownames(A) <- colnames(A) <- shp$user_state
  list(A = A, shp = shp)
}

CARData4Stan <- function(NeighborhoodMatrix) {
  N <- nrow(NeighborhoodMatrix)
  A <- (NeighborhoodMatrix > 0) * 1L
  diag(A) <- 0L
  A_ut <- A * (row(A) < col(A))
  edge_idx <- which(A_ut == 1L, arr.ind = TRUE)
  node1 <- as.integer(edge_idx[, 1])
  node2 <- as.integer(edge_idx[, 2])
  list(N = N, N_edges = length(node1), node1 = node1, node2 = node2)
}

ScalingFacBYM2 <- function(Nodes, AdjMat) {
  Q <- Diagonal(Nodes$N, rowSums(AdjMat)) - AdjMat
  Q_pert <- Q + Diagonal(Nodes$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  Q_inv <- INLA::inla.qinv(Q_pert, constr = list(A = matrix(1, 1, Nodes$N), e = 0))
  exp(mean(log(diag(Q_inv))))
}

## ---- clean df and build keys ----
df <- df %>%
  mutate(country = tolower(country),
         state   = tolower(state),
         sex     = tolower(sex)) %>%
  arrange(country, state, year, age)

stopifnot(all(c("colombia","mexico") %in% unique(df$country)))

loc_key <- df %>% distinct(country, state) %>%
  arrange(country, state) %>%
  mutate(l_id = row_number(),
         c_id = as.integer(factor(country, levels = c("colombia","mexico"))))

sex_key  <- df %>% distinct(sex)  %>% arrange(sex)  %>% mutate(s_id = row_number())
year_key <- df %>% distinct(year) %>% arrange(year) %>% mutate(t_id = row_number())
age_key  <- df %>% distinct(age)  %>% arrange(age)  %>% mutate(a_id = row_number())

dfi <- df %>%
  inner_join(loc_key, by = c("country","state")) %>%
  inner_join(sex_key,  by = "sex") %>%
  inner_join(year_key, by = "year") %>%
  inner_join(age_key,  by = "age") %>%
  arrange(l_id, t_id, s_id, a_id)

## ---- panel ids (not used by Stan now but kept) ----
panel_key <- dfi %>%
  distinct(l_id, t_id) %>%
  arrange(l_id, t_id) %>%
  mutate(panel_id = row_number())

dfi <- dfi %>% left_join(panel_key, by = c("l_id","t_id"))

## ---- BYM2 inputs (per country) ----
locs_col <- loc_key %>% filter(country == "colombia") %>% arrange(state)
locs_mx  <- loc_key %>% filter(country == "mexico")   %>% arrange(state)

user_states_col <- locs_col %>% pull(state)
user_states_mx  <- locs_mx  %>% pull(state)

co_out <- build_A_for_country("Colombia", user_states_col)
mx_out <- build_A_for_country("Mexico",   user_states_mx)

A_col <- co_out$A; A_mx <- mx_out$A
nodes_col <- CARData4Stan(A_col)
nodes_mx  <- CARData4Stan(A_mx)
sf_col <- ScalingFacBYM2(nodes_col, A_col)
sf_mx  <- ScalingFacBYM2(nodes_mx,  A_mx)

L1 <- nrow(locs_col); L2 <- nrow(locs_mx)
loc_ids_c1 <- as.integer(locs_col$l_id)
loc_ids_c2 <- as.integer(locs_mx$l_id)

## ---- Age cubic spline with 3 internal knots (columns centered) ----
ages <- sort(unique(dfi$age))
ik <- as.numeric(quantile(ages, probs = c(0.25, 0.50, 0.75), type = 1))
ik <- ik[ik > min(ages) & ik < max(ages)]
stopifnot(length(ik) >= 3)
B_age <- splines::bs(x = ages, degree = 3,
                     knots = ik,
                     Boundary.knots = range(ages),
                     intercept = FALSE)
B_age <- scale(B_age, center = TRUE, scale = FALSE)
B_age <- as.matrix(B_age)
A_len <- length(ages); K <- ncol(B_age)

## ---- Observed vs missing masks ----
is_obs <- !is.na(dfi$deaths)
obs_idx <- which(is_obs)
mis_idx <- which(!is_obs)

N_all <- nrow(dfi)
Nobs  <- length(obs_idx)
Nmis  <- length(mis_idx)
stopifnot(Nobs > 0)

## ---- Design & scaling: use scale() with stats from observed rows ----

# propyth
x1_mu <- mean(dfi$propyth[is_obs], na.rm = TRUE)
x1_sd <- sd(dfi$propyth[is_obs],   na.rm = TRUE)
if (!is.finite(x1_sd) || x1_sd < .Machine$double.eps) x1_sd <- 1
x1 <- as.numeric(scale(dfi$propyth, center = x1_mu, scale = x1_sd))

# pov
x2_mu <- mean(dfi$pov[is_obs], na.rm = TRUE)
x2_sd <- sd(dfi$pov[is_obs],   na.rm = TRUE)
if (!is.finite(x2_sd) || x2_sd < .Machine$double.eps) x2_sd <- 1
x2 <- as.numeric(scale(dfi$pov, center = x2_mu, scale = x2_sd))

# ---- UPDATED: W uses completeness (propcomp_bothsex), not missingness ----
w_mu <- mean(dfi$propcomp_bothsex[is_obs], na.rm = TRUE)
w_sd <- sd(dfi$propcomp_bothsex[is_obs],   na.rm = TRUE)
if (!is.finite(w_sd) || w_sd < .Machine$double.eps) w_sd <- 1
w1 <- as.numeric(scale(dfi$propcomp_bothsex, center = w_mu, scale = w_sd))

# X matrix for ALL rows: linear terms only
Xmat_all <- cbind(
  propyth = x1,
  pov     = x2
)
Xmat_all <- as.matrix(Xmat_all)
p <- ncol(Xmat_all)   # p = 2

# W rows for reporting model (linear only)
Wrow_all <- w1

## ---- Indices and offsets for ALL rows ----
ii_l_all <- as.integer(dfi$l_id)
ii_s_all <- as.integer(dfi$s_id)
ii_t_all <- as.integer(dfi$t_id)
ii_a_all <- as.integer(dfi$a_id)

D_obs     <- as.integer(dfi$deaths[obs_idx])
log_P_all <- log(pmax(dfi$pop, 1e-12))

## ---- Intercept prior for deaths ----
mu_mean <- -9
mu_sd   <- 0.5

## ===================================================================
## ====  CR → Beta prior for p0 via RE meta-analysis on logit scale ===
## ===================================================================

logit  <- function(p) log(p / (1 - p))
ilogit <- function(x) 1 / (1 + exp(-x))
eps <- 1e-8

rp2 <- rp %>%
  mutate(
    p      = pmin(pmax(prop_completeness, eps), 1 - eps),
    logitp = logit(p)
  )

cr_units <- rp2 %>%
  group_by(state, year) %>%
  summarise(
    yi      = mean(logitp),
    vi      = var(logitp),
    n_draws = dplyr::n(),
    .groups = "drop"
  )

stopifnot(nrow(cr_units) > 0)

fit_meta <- metafor::rma.uni(yi = yi, vi = vi, method = "REML", data = cr_units)

mu_alpha <- as.numeric(fit_meta$b)
var_mu_a <- as.numeric(vcov(fit_meta))

set.seed(123)
B <- 20000
alpha_draws <- rnorm(B, mean = mu_alpha, sd = sqrt(var_mu_a))
p_draws     <- ilogit(alpha_draws)

mu_prob <- mean(p_draws)
v_prob  <- var(p_draws)
v_prob  <- max(min(v_prob, mu_prob * (1 - mu_prob) - 1e-12), 1e-12)

kappa_raw <- mu_prob * (1 - mu_prob) / v_prob - 1
kappa <- min(max(kappa_raw, 2), 1000)  # temper ESS

a_alpha <- as.numeric(mu_prob * kappa)
b_alpha <- as.numeric((1 - mu_prob) * kappa)

cat("Beta prior for p0 (from CR, RE meta-analysis):\n")
cat(sprintf("  pooled mean (logit) = %.4f\n", mu_alpha))
cat(sprintf("  pooled mean (prob)  = %.4f\n", mu_prob))
cat(sprintf("  ESS raw = %.1f, tempered = %.1f\n", kappa_raw, kappa))
cat(sprintf("  a = %.3f, b = %.3f\n", a_alpha, b_alpha))
cat(sprintf("  95%% prior interval ≈ [%.3f, %.3f]\n",
            qbeta(0.025, a_alpha, b_alpha), qbeta(0.975, a_alpha, b_alpha)))

## ---- Assemble stan_data ----
stan_data <- list(
  # sizes
  L  = nrow(loc_key),
  S  = nrow(sex_key),
  T  = nrow(year_key),
  A  = A_len,
  K  = K,
  p  = p,
  
  # counts
  N_all  = N_all,
  Nobs   = Nobs,
  obs_idx = as.integer(obs_idx),
  D_obs   = as.integer(D_obs),
  Nmis    = Nmis,
  mis_idx = as.integer(mis_idx),
  
  # ALL rows
  ii_l_all  = as.integer(ii_l_all),
  ii_s_all  = as.integer(ii_s_all),
  ii_t_all  = as.integer(ii_t_all),
  ii_a_all  = as.integer(ii_a_all),
  log_P_all = as.numeric(log_P_all),
  Xmat_all  = Xmat_all,
  Wrow_all  = as.numeric(Wrow_all),
  
  # age basis
  B_age = B_age,
  
  # intercept prior
  mu_mean = mu_mean,
  mu_sd   = mu_sd,
  
  # Beta prior hyperparameters for p0
  a_alpha = as.numeric(a_alpha),
  b_alpha = as.numeric(b_alpha),
  
  # BYM2 inputs (per country)
  L1 = as.integer(L1),
  E1 = as.integer(nodes_col$N_edges),
  node1_c1   = as.integer(nodes_col$node1),
  node2_c1   = as.integer(nodes_col$node2),
  loc_ids_c1 = as.integer(loc_ids_c1),
  scale_icar_c1 = as.numeric(sf_col),
  
  L2 = as.integer(L2),
  E2 = as.integer(nodes_mx$N_edges),
  node1_c2   = as.integer(nodes_mx$node1),
  node2_c2   = as.integer(nodes_mx$node2),
  loc_ids_c2 = as.integer(loc_ids_c2),
  scale_icar_c2 = as.numeric(sf_mx)
)

## ---- Compile & sample ----
stan_file <- "model_BYM2.stan"

sm <- rstan::stan_model(file = stan_file)

fit <- rstan::sampling(
  sm, data = stan_data,
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 12),
  seed = 123
)

print(fit, pars = c("mu","alpha","p0","beta","gamma","rho",
                    "sigma_s","sigma_t","sigma_gamma",
                    "sigma_c", "sigma_beta","sigma_rw2"),
      probs = c(.1,.5,.9))

#------------------------------------------------------------------------------#
# Diagnostic Plots
#------------------------------------------------------------------------------#
library(bayesplot)
# ---- Parameters to check ----
pars_main <- c(
  "mu", "alpha",
  "beta[1]", "beta[2]", "gamma", 
  "sigma_s", "sigma_t",
  "sigma_c", "sigma_rw2",
  "rho"
)

# ---- Traceplot ----
posterior <- rstan::extract(fit, pars = pars_main)
mcmc_trace(as.array(fit), pars = pars_main)

# ---- Pairs plot ----
mcmc_pairs(as.array(fit), pars = pars_main)

#------------------------------------------------------------------------------#
# Construct Posterior Data Frame (only 50 draws)
#------------------------------------------------------------------------------#

library(dplyr)
library(tidyr)
library(tibble)

stopifnot(exists("fit"))
stopifnot(exists("dfi"))

# ---- Extract full posterior arrays first ----
post_full <- rstan::extract(
  fit,
  pars = c("pi_row_all",
           "mu_true_all",
           "mu_reported_all",
           "y_true_all",
           "y_rep_all")
)

n_draws_full <- dim(post_full$pi_row_all)[1]
N_all        <- nrow(dfi)

# ---- Sample 50 draws ----
set.seed(123)
keep_draws <- sample(seq_len(n_draws_full), size = 50, replace = FALSE)

# ---- Subset arrays to 50 draws ----
post <- list(
  pi_row_all      = post_full$pi_row_all[keep_draws, ],
  mu_true_all     = post_full$mu_true_all[keep_draws, ],
  mu_reported_all = post_full$mu_reported_all[keep_draws, ],
  y_true_all      = post_full$y_true_all[keep_draws, ],
  y_rep_all       = post_full$y_rep_all[keep_draws, ]
)

n_draws <- length(keep_draws)  # now = 50

# ---- Vectorize across draws × rows ----
pi_vec      <- as.vector(t(post$pi_row_all))
mu_true_vec <- as.vector(t(post$mu_true_all))
mu_rep_vec  <- as.vector(t(post$mu_reported_all))
y_true_vec  <- as.vector(t(post$y_true_all))
y_rep_vec   <- as.vector(t(post$y_rep_all))

draw_index <- rep(seq_len(n_draws), each = N_all)
row_index  <- rep(seq_len(N_all),   times = n_draws)

# ---- Observed indicator + observed deaths ----
is_obs     <- !is.na(dfi$deaths)
deaths_obs <- dfi$deaths

# ---- Build draw-level dataframe ----
draws_df <- tibble(
  draw    = draw_index,
  row     = row_index,
  
  country = rep(dfi$country, times = n_draws),
  state   = rep(dfi$state,   times = n_draws),
  sex     = rep(dfi$sex,     times = n_draws),
  year    = rep(dfi$year,    times = n_draws),
  age     = rep(dfi$age,     times = n_draws),
  
  is_obs     = rep(is_obs,     times = n_draws),
  deaths_obs = rep(deaths_obs, times = n_draws),
  
  pi      = pi_vec,
  mu_true = mu_true_vec,
  mu_rep  = mu_rep_vec,
  y_true  = y_true_vec,
  y_rep   = y_rep_vec
)

# ---- Sanity checks ----
print(nrow(draws_df))         # should be 50 * N_all
print(n_draws * N_all)
draws_df %>% count(draw)      # each draw should have N_all rows

# ---- Save ----
save(draws_df, file = "posteriorBYM2_50draws.RData")

#------------------------------------------------------------------------------#
# Fin.
#------------------------------------------------------------------------------#

