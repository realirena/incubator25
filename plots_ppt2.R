
rm(list=ls())
setwd("N:/Incubator2025_Conflict/Figures/final_plots")

library(tidyverse)
library(readxl)

# Load


# load statenames 
load("mx_statesnames.rdata")

lookup <- mx_statenames$state
names(lookup) <- mx_statenames$statenum

# load data 
data <- read_csv("posterior_estimates.csv")

# load age proportions 
#aged <- read_xlsx("agepropmexico.xlsx")
load("ageprop.rdata")
aged <- ageprop

# load observed homicides
load("N:/Incubator2025_Conflict/Figures/hom_mx.rdata")

# Convert to state names
aged$state <- lookup[as.numeric(aged$statenum)]
aged <- aged%>% select(-statenum)
#data$death <- data$Z_rep
data$death <- data$lambda
#data$state <- data$state_name
data$draw <- data$.draw


# Load population
load("m_agestructure.RData")

pop <- 
  long_mx_pop %>% 
  select(-pob_total) %>% 
  pivot_longer(cols = -c("state", "sex", "year")) %>% 
  group_by(state, year, name) %>% 
  summarise(
    value = sum(value, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  mutate(
    name = gsub("pob_|_", "", name)
    , name = gsub("([[:digit:]]{2})(\\d+)", "\\1-\\2", name)
    , name = gsub("(?<=-)0+(\\d+)", "\\1", name, perl = TRUE)
    , name = ifelse(name == "00-4", "0-4", name)
    , name = ifelse(name == "05-9", "5-9", name) 
    , name = ifelse(name == "85mm", "85", name) 
  ) %>% 
  rename(pop = value, age = name)
  
  
# 2. Create data structure ------------

df_all <- 
  # Create all combinations 
  crossing(
  state = unique(data$state),
  year = unique(data$year),
  draw = unique(data$draw),
  age = unique(aged$age)
) %>% 
  # Add deaths
  left_join(data, by = c("draw", "state", "year")) %>%
  # Add proportions
  left_join(aged, by = c("age", "state", "year")) %>% 
  # Add the population counts
  left_join(pop, by = c("age", "state", "year")) %>%
  mutate(
    # estimate age-specific deaths
    death2 = death * age_prop
    # , death2 = round(death)
    # Get rates
    , rate = death2 / pop
    ) # %>% 
  # filter(!is.na(age_prop)) %>% 
  # select(-redistributed, -age_prop)
  # select(-age_prop)

df_all$death2 <- ifelse(is.na(df_all$death2), 0, df_all$death2)
df_all$rate <- ifelse(is.na(df_all$rate), 0, df_all$rate)

  

# Get percentiles

df_perc <- 
  df_all %>% 
  group_by(state, age, year) %>% 
  summarise(
    median = quantile(rate, .5, na.rm = T)
    , low = quantile(rate, .01, na.rm = T)
    , high = quantile(rate, .09, na.rm = T)
  ) %>% 
  ungroup()

# Plot ----

ages_lev <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
          "40-44", "45-49",  "50-54", "55-59", "60-64", "65-69", 
          "70-74", "75-79")

df_perc %>% 
  filter(state == "chihuahua") %>% 
  mutate(age = factor(age, levels = ages_lev)) %>% 
  ggplot(aes(x = year, y = age, fill = median)) + 
  geom_tile() + 
  scale_fill_viridis_c(direction = -1) +
  theme_bw()
  
# geohamp
# 15-29: observed and estiamted by geomape (focus on median)

# overall heatmap
library(dplyr)
library(ggplot2)
library(geofacet)

df_perc2 <- df_perc %>%
  left_join(mx_statenames, by = "state") %>%
  select(-geom) %>%
  mutate(
    NAME_1 = recode(NAME_1,
                      "Distrito Federal" = "Ciudad de México"),
    age = factor(age, levels = ages_lev)
  ) %>%
  filter(!is.na(NAME_1))


# 5) heatmap Plot
df_perc2 %>%
  # filter(age %in% c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
  #                   "40-44", "45-49")) %>%
  filter(year>2007) %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = year, y = age, fill = median)) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() +
  facet_geo(~ NAME_1, grid = "mx_state_grid3")

df_perc2 %>% 
  filter(NAME_1 == "Chihuahua") %>% 
  mutate(age = factor(age, levels = ages_lev)) %>% 
  ggplot(aes(x = year, y = age, fill = median)) + 
  geom_tile() + 
  scale_fill_viridis_c(direction = -1) +
  theme_bw()


###############
# Recompute percentiles using death COUNTS (death2), not rates
df_perc_deaths <- df_all %>%
  group_by(state, age, year) %>%
  summarise(
    median_deaths = quantile(death2, 0.5, na.rm = TRUE),
    low  = quantile(death2, 0.10, na.rm = TRUE),  
    high = quantile(death2, 0.90, na.rm = TRUE),  
    .groups = "drop"
  )

# Filter to Chihuahua & age 25-29 (case-robust on state)
plot_dat <- df_perc_deaths %>%
  filter(tolower(state) == "chihuahua", age == "25-29") %>%
  arrange(year)

# Line plot: Y = median death count, X = year
ggplot(plot_dat, aes(x = year, y = median_deaths)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(size = 1) +
  geom_point() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Median Estimated Deaths — Chihuahua, Age 25–29",
    x = "Year",
    y = "Median deaths (count)",
    caption = "Ribbon shows 10th–90th percentile across posterior draws"
  ) +
  theme_bw()



###
library(dplyr)
library(ggplot2)

# 1) Filter to Chihuahua & age 25–29; keep draw-level death COUNTS
df_chih_25_29 <- df_all %>%
  filter(tolower(state) == "chihuahua", age == "25-29") %>%
  select(year, draw, death2)

# 2) Sample up to 50 distinct draws (handle cases with < 50 draws)
n_draws <- df_chih_25_29 %>% distinct(draw) %>% nrow()

set.seed(123) 
sampled_draws <- df_chih_25_29 %>%
  distinct(draw) %>%
  slice_sample(n = min(50, n_draws))

# 3) Keep only sampled draws
df_sampled <- df_chih_25_29 %>%
  semi_join(sampled_draws, by = "draw") %>%
  arrange(draw, year)

# 4) Plot: transparent lines for each draw + LOESS trend
p_loess <- ggplot(df_sampled, aes(x = year, y = death2, group = draw)) +
  geom_line(alpha = 0.2, color = "blue") +  # <- light blue lines
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 1.2) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Posterior Draws (≤50) with LOESS — Chihuahua, Age 25–29",
    x = "Year",
    y = "Deaths (count)"
  ) +
  theme_minimal()

print(p_loess)


#####
library(dplyr)
library(ggplot2)
library(geofacet)
library(stringr)

# Build a clean state-name map from your lookup (mx_statenames was already loaded earlier)
state_map <- mx_statenames %>%
  distinct(state, NAME_1) %>%
  mutate(
    # Fix the grid mismatch
    NAME_1 = dplyr::recode(NAME_1,
                             "Distrito Federal" = "Ciudad de México")
  )

# Keep draw-level death COUNTS and the grid facet name
df_age <- df_all %>%
  filter(age == "25-29") %>%
  left_join(state_map, by = "state") %>%
  filter(!is.na(NAME_1)) %>%     # drop anything that can't facet
  select(state, NAME_1, year, draw, death2)

# ----- 2) Sample up to 50 draws PER STATE to keep things readable -----

# Distinct draws per state, then sample within each state
set.seed(123)  # reproducibility
# Sample up to 50 draws PER STATE safely
sampled_draws <- df_age %>%
  dplyr::distinct(NAME_1, draw) %>%
  dplyr::group_by(NAME_1) %>%
  dplyr::reframe(draw = sample(draw, size = min(50, dplyr::n()), replace = FALSE)) %>%
  dplyr::ungroup()

# Keep only sampled draws and plot
df_sampled <- df_age %>%
  dplyr::semi_join(sampled_draws, by = c("NAME_1", "draw")) %>%
  dplyr::arrange(NAME_1, draw, year) %>%
  filter(NAME_1 != "México") %>%
  filter(NAME_1 != "Ciudad de México") %>%
  filter(NAME_1 != "Jalisco") 

ggplot(df_sampled, aes(x = year, y = death2, group = draw)) +
  geom_line(alpha = 0.2, color = "lightblue") +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 0.8, span = 0.6) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Posterior Draws (≤50 per state) with LOESS — Age 25–29",
    subtitle = "Each panel is a state; lines are sampled posterior draws of deaths (counts)",
    x = "Year", y = "Deaths (count)"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8)) +
  geofacet::facet_geo(~ NAME_1, grid = "mx_state_grid3")
table(df_sampled$NAME_1)



## add observed homicides 

## filter observed homicides to "25 to 29"
hom_mx_25 <- hom_mx %>% filter(agecat == "25 to 29") %>%
  group_by(statenum, year) %>%
  summarise(observed = sum(Homicides))

hom_mx_25 <- hom_mx_25 %>% left_join( mx_statenames, by = "statenum")


###
library(dplyr)
library(ggplot2)
library(geofacet)
library(stringr)

# ----- 1) Prep: filter to age 25–29 and attach geofacet-friendly names -----
state_map <- mx_statenames %>%
  distinct(state, NAME_1) %>%
  mutate(
    NAME_1 = dplyr::recode(NAME_1, "Distrito Federal" = "Ciudad de México")
  )

df_age <- df_all %>%
  filter(age == "25-29") %>%
  left_join(state_map, by = "state") %>%
  filter(!is.na(NAME_1)) %>%
  select(state, NAME_1, year, draw, death) %>%
  filter(NAME_1 %in% c("Chihuahua"))

# ----- 2) Sample up to 50 draws PER STATE -----
set.seed(123)
sampled_draws <- df_age %>%
  distinct(NAME_1, draw) %>%
  group_by(NAME_1) %>%
  reframe(draw = sample(draw, size = min(50, n()), replace = FALSE)) %>%
  ungroup()

df_sampled <- df_age %>%
  semi_join(sampled_draws, by = c("NAME_1", "draw")) %>%
  arrange(NAME_1, draw, year)

# ----- 3) Observed count for age 25–29 -----
hom_mx_25 <- hom_mx %>%
  filter(agecat == "25 to 29") %>%
  filter(year > 2006 ) %>%
  group_by(statenum, year) %>%
  summarise(observed = sum(Homicides), .groups = "drop") %>%
  left_join(mx_statenames, by = "statenum") %>%
  mutate(
    NAME_1 = dplyr::recode(NAME_1, "Distrito Federal" = "Ciudad de México")
  ) %>%
  filter(NAME_1 %in% c("Chihuahua"))

# ----- 4) Plot -----
ggplot(df_sampled, aes(x = year, y = death, group = draw)) +
  geom_line(alpha = 0.2, color = "blue") +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "black", size = 0.8, span = 0.6) +
  geom_line(data = hom_mx_25, aes(x = year, y = observed, group = NAME_1),
            inherit.aes = FALSE, color = "red", size = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Chihuahua: Age 25–29",
    subtitle = "Blue: posterior draws, Black: LOESS trend, Red: observed counts",
    x = "Year", y = "Deaths (count)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )



## try geofacet

# sampled_draws2 <- df_all %>%
#   left_join(state_map, by = "state") %>%
#   distinct(NAME_1, draw) %>%
#   group_by(NAME_1) %>%
#   reframe(draw = sample(draw, size = min(50, n()), replace = FALSE)) %>%
#   ungroup()
# 
# df_sampled2 <- df_all %>%
#   left_join(state_map, by = "state") %>%
#   semi_join(sampled_draws2, by = c("NAME_1", "draw")) %>%
#   arrange(NAME_1, draw, year)
# 
# # ----- 3) Observed count for age 25–29 -----
# hom_mx_25 <- hom_mx %>%
#   filter(agecat == "25 to 29") %>%
#   filter(year > 2006 ) %>%
#   group_by(statenum, year) %>%
#   summarise(observed = sum(Homicides), .groups = "drop") %>%
#   left_join(mx_statenames, by = "statenum") %>%
#   mutate(
#     NAME_1 = dplyr::recode(NAME_1, "Distrito Federal" = "Ciudad de México")
#   ) %>%
#   filter(NAME_1 %in% c("Chihuahua"))
# 
# # ----- 4) Plot -----
# ggplot(df_sampled2, aes(x = year, y = death, group = draw)) +
#   geom_line(alpha = 0.2, color = "blue") +
#   geom_smooth(aes(group = 1), method = "loess", se = FALSE,
#               color = "black", size = 0.8, span = 0.6) +
#   #geom_line(data = hom_mx, aes(x = year, y = Homicides, group = NAME_1),
#   #          inherit.aes = FALSE, color = "red", size = 0.8) +
#   scale_y_continuous(labels = scales::comma) +
#   facet_geo(~ NAME_1, grid = "mx_state_grid3") +
#   labs(
#     title = "Chihuahua: Age 25–29",
#     subtitle = "Blue: posterior draws, Black: LOESS trend, Red: observed counts",
#     x = "Year", y = "Deaths (count)"
#   ) +
#   theme_minimal(base_size = 10) +
#   theme(
#     panel.grid.minor = element_blank(),
#     strip.text = element_text(size = 10)
#   ) 


# ----- 1) Prep: filter to age 25–29 and attach geofacet-friendly names -----
state_map <- mx_statenames %>%
  distinct(state, NAME_1) %>%
  mutate(
    NAME_1 = dplyr::recode(NAME_1, "Distrito Federal" = "Ciudad de México")
  )

df_plot <- df_all %>%
  #filter(age == "25-29") %>%
  left_join(state_map, by = "state") %>%
  filter(!is.na(NAME_1)) %>%
  select(state, NAME_1, year, draw, death) %>%
  filter(NAME_1 %in% c("Oaxaca"))

# ----- 2) Sample up to 50 draws PER STATE -----
set.seed(123)
sampled_draws3 <- df_plot %>%
  distinct(NAME_1, draw) %>%
  group_by(NAME_1) %>%
  reframe(draw = sample(draw, size = min(50, n()), replace = FALSE)) %>%
  ungroup()

df_sampled3 <- df_plot %>%
  semi_join(sampled_draws3, by = c("NAME_1", "draw")) %>%
  arrange(NAME_1, draw, year)

# ----- 3) Observed count for age 25–29 -----
hom_mx_all <- hom_mx %>%
  #filter(agecat == "25 to 29") %>%
  filter(year > 2006 ) %>%
  group_by(statenum, year) %>%
  summarise(observed = sum(Homicides), .groups = "drop") %>%
  left_join(mx_statenames, by = "statenum") %>%
  mutate(
    NAME_1 = dplyr::recode(NAME_1, "Distrito Federal" = "Ciudad de México")
  ) %>%
  filter(NAME_1 %in% c("Oaxaca"))

# ----- 4) Plot -----
ggplot(df_sampled3, aes(x = year, y = death, group = draw)) +
  geom_line(alpha = 0.2, color = "blue") +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "black", size = 0.8, span = 0.6) +
  geom_line(data = hom_mx_all, aes(x = year, y = observed, group = NAME_1),
            inherit.aes = FALSE, color = "red", size = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Oaxaca: Age: all",
    subtitle = "Blue: posterior draws, Black: LOESS trend, Red: observed counts",
    x = "Year", y = "Deaths (count)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )
