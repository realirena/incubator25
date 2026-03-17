pacman::p_load(tidyverse, sf, spdep, geodata, purrr, splines,
               INLA, cmdstanr, Matrix, posterior, bayesplot)

load("data/df.Rdata")

# Indices
df1 <- df %>% 
  filter(sex == "male") %>%
  mutate(age = as.numeric(sub(" .*","", age))) %>% 
  arrange(state, year, age) %>% 
  mutate(
    an = match(age, sort(unique(age))),
    tn = match(year, sort(unique(year))),
    ln = match(state, sort(unique(state)))
  )

# BYM2 model 
# Transform Neighborhood Matrix into List for Stan
CARData4Stan <- function(NeighborhoodMatrix) {
  
  # Get upper triangular part (excluding diagonal)
  upper_tri <- as.matrix(NeighborhoodMatrix) * 
    upper.tri(NeighborhoodMatrix, diag = FALSE) 
  
  
  N <- nrow(NeighborhoodMatrix) #Amount of Regions
  N_edges <- sum(upper_tri) #Amount of Edges
  
  # Get row and column indices of 1s
  edges <- which(upper_tri == 1, 
                 arr.ind = TRUE)
  
  #save results in vector (orderd)
  node1 <- edges[order(edges[,1]),1]
  node2 <- edges[order(edges[,2]),2]
  
  node1Names <- rownames(NeighborhoodMatrix)[node1] 
  node2Names <- rownames(NeighborhoodMatrix)[node2] 
  
  return(list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2,
              "node1Names"=node1Names, node2Names=node2Names))
}
#Compute Scaling Factor for BYM2 Model
#Code taken from Morris (2019) 
#https://mc-stan.org/users/documentation/case-studies/icar_stan.html
ScalingFacBYM2 <- function(Nodes, AdjMat){
  #Add a small jitter to the diagonal for numerical stability 
  Q <- Diagonal(Nodes$N, rowSums(AdjMat)) - AdjMat
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_pert <-  Q + Diagonal(Nodes$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the 
  # constraint that the entries of the ICAR sum to zero.
  #See the inla.qinv function help for further details.
  Q_inv <-  INLA::inla.qinv(Q_pert, constr=list(A = matrix(1,1,Nodes$N),e=0))
  
  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
  scaling_factor <-  exp(mean(log(diag(Q_inv))))
  
  return(scaling_factor)
}

# polygons (states)
mx <- st_as_sf(geodata::gadm("MEX", level = 1, path = tempdir()))
mx <- st_make_valid(mx)
# make names match the data
mx$state <- tolower(iconv(mx$NAME_1, from="UTF-8", to="ASCII//TRANSLIT"))
mx$state[mx$state == "distrito federal"] <- "mexico city"
mx$state[mx$state == "mexico"]           <- "mexico state"
# reorder polygons to match the data
ord <- unique(df1$state)
mx  <- mx[match(ord, mx$state), ]
stopifnot(all(mx$state == ord))
# adjacency
nb  <- poly2nb(mx, queen = TRUE)
AdjMatMexico <- as(nb2mat(nb, style="B", zero.policy=TRUE), "Matrix")
#Get Nodes for Stan
Nodes <- CARData4Stan(NeighborhoodMatrix = AdjMatMexico) 
#scaling Factor BYM2 Model
scalingFac <- ScalingFacBYM2(Nodes = Nodes, AdjMat = AdjMatMexico) 

# age spline
ages <- sort(unique(df1$age))
internal_knots = seq(15,75, by = 15)
boundary_knots = c(0,90)
B <- bs(
  x = ages,
  knots = internal_knots,
  Boundary.knots = boundary_knots,
  degree = 3,
  intercept = FALSE
)
B <- scale(B, center = TRUE, scale = TRUE) # zero-mean spline
B <- as.matrix(B)
K <- ncol(B)
dim(B)

# observation index
obs <- which(!is.na(df1$deaths))
Nobs <- length(obs)

# covariates
disap <- scale(df1$disap, center = TRUE, scale = TRUE)
pov <- scale(df1$pov, center = TRUE, scale = TRUE)

X <- cbind(disap, pov)
X <- as.matrix(X)


# stan data
stan_data <- list(
  N = nrow(df1),
  A = max(df1$an),
  T = max(df1$tn),
  L = max(df1$ln),
  an = df1$an,
  tn = df1$tn,
  ln = df1$ln,
  obs = obs,
  Nobs = Nobs, 
  
  y = as.integer(df1$deaths[obs]),
  E = pmax(df1$pop, 1e-6),
  
  X = X,
  K = K,
  B = B,
  
  N_edges = Nodes$N_edges,
  node1 = Nodes$node1,
  node2 = Nodes$node2,
  scaling_factor = scalingFac
)

# compile model
mod <- cmdstan_model("Other/model_BYM2.stan")

# fit model
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.90,
  max_treedepth = 10
)

# diagnostic
fit$cmdstan_diagnose()

# posterior data
vars <- c("rep_prob","lambda_true","lambda_obs","y_true","y_rep")
draws  <- as_draws_matrix(fit$draws(vars))
ndraw <- nrow(draws)
N     <- stan_data$N 
# helper: pull var[1:N] columns in correct numeric order
pull_var <- function(var, draws, N){
  cols <- grep(paste0("^", var, "\\["), colnames(draws), value = TRUE)
  idx  <- as.integer(sub(".*\\[(\\d+)\\].*", "\\1", cols))
  cols <- cols[order(idx)]
  stopifnot(length(cols) == N)
  draws[, cols, drop = FALSE]  # ndraw x N
}
# Make long posterior df: (draw x n) with 5 columns
rp <- pull_var("rep_prob",    draws, N)
lt <- pull_var("lambda_true", draws, N)
lo <- pull_var("lambda_obs",  draws, N)
yt <- pull_var("y_true",      draws, N)
yr <- pull_var("y_rep",       draws, N)
df_post <- tibble(
  .draw       = rep(seq_len(ndraw), each = N),
  n           = rep(seq_len(N), times = ndraw),
  rep_prob    = as.vector(t(rp)),
  lambda_true = as.vector(t(lt)),
  lambda_obs  = as.vector(t(lo)),
  y_true      = as.integer(round(as.vector(t(yt)))),
  y_rep       = as.integer(round(as.vector(t(yr))))
)
# Join to df1
df1 <- df1 %>% mutate(n = row_number())
post_df <- df_post %>% left_join(df1, by = "n")

# select 100 random draws
set.seed(1234)
draw_ids <- sample(unique(post_df$.draw), 100)
post_df_100 <- post_df %>%
  filter(.draw %in% draw_ids)
# append posteriors for males and females
post_df_m <- post_df_100
# post_df_f <- post_df_100
# post_draws <- rbind.data.frame(post_df_m, post_df_f)
# save(post_draws, file = "results/post_draws.RData")

# PPC
y_rep <- as_draws_matrix(fit$draws("y_rep"))
y_obs <- df1$deaths[obs]

p1 <- ppc_dens_overlay(y_obs, y_rep[1:200, obs])
p2 <- ppc_ecdf_overlay(y_obs, y_rep[1:200, obs])
p3 <- ppc_stat(y_obs, y_rep[, obs], stat = "mean")
p4 <- ppc_stat(y_obs, y_rep[, obs], stat = "sd")
p5 <- ppc_stat(y_obs, y_rep[, obs], stat = function(x) mean(x == 0))

ggsave("figure/A1_dens_overlay_m.png", p1, height = 16, width = 18, dpi = 300)
ggsave("figure/A1_dens_overlay_m.pdf", p1, height = 16, width = 18, dpi = 300)
ggsave("figure/A2_ecdf_overlay_m.png", p2, height = 16, width = 18, dpi = 300)
ggsave("figure/A2_ecdf_overlay_m.pdf", p2, height = 16, width = 18, dpi = 300)
ggsave("figure/A3_mean_m.png", p3, height = 16, width = 18, dpi = 300)
ggsave("figure/A3_mean_m.pdf", p3, height = 16, width = 18, dpi = 300)
ggsave("figure/A4_sd_m.png", p4, height = 16, width = 18, dpi = 300)
ggsave("figure/A4_sd_m.pdf", p4, height = 16, width = 18, dpi = 300)
ggsave("figure/A5_prop_zero_m.png", p5, height = 16, width = 18, dpi = 300)
ggsave("figure/A5_prop_zero_m.pdf", p5, height = 16, width = 18, dpi = 300)
