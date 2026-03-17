pacman::p_load(tidyverse, data.table, future.apply)
# cause of death data
load("data/oth_mx.RData")
colnames(deaths)[5] <- "other"
deaths$age <- as.numeric(sub(" .*", "", deaths$age))
# posteriors
load("results/post_draws.Rdata")
colnames(post_draws)[1] <- "draw"

df1 <- merge(post_draws, deaths, 
              by = c("state", "sex", "year", "age"),
              all.x = T)
df <- df1 %>% 
  mutate(mx_h = lambda_true/pop,
         mx_o = other/pop) %>% 
  select(state, sex, year, age, draw, mx_h, mx_o) %>% 
  arrange(state, sex, year, age, draw)

source("01_functions.R")

setDT(df)
causes <- c("mx_h", "mx_o")
ages <- sort(unique(df$age))
years <- sort(unique(df$year))

groups <- split(df, by = c("state", "sex", "draw"), keep.by = TRUE)
# future::plan(multisession, workers = parallelly::availableCores()-1)
future::plan(multisession, workers = 32)
decomp_result <- rbindlist(
  future_lapply(
    groups,
    decomp_grp,
    N = 50,
    causes = causes,
    future.seed = TRUE
  ),
  use.names = TRUE,
  fill = TRUE
)
future::plan(sequential)
save(decomp_result, file = "results/decomp_result.RData")

# contribution to change in e(0) between 2007-2023
decomp_sum <- decomp_result %>% 
  group_by(state, sex, draw, age, cause) %>% 
  summarise(contribution = sum(contribution, na.rm = TRUE), 
            .groups = "drop") %>% 
  arrange(state, sex, draw, cause)

