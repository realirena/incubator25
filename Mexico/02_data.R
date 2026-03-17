Sys.setenv(LANGUAGE = "en")
pacman::p_load(dplyr, tidyr, stringr)

# get data
# deaths 
load("data/hom_mx.rdata")
df <- hom_mx %>% 
  select(state, sex, year, age, deaths) %>% 
  mutate(sex=ifelse(sex==0,"male", "female")) %>% 
  arrange(state, sex, year, age)

# Joint redistribution within year to observed (state, sex, age) cells
# then drop NAs
df_clean <- df %>%
  mutate(deaths = replace_na(deaths, 0)) %>%
  group_by(year) %>%
  mutate(miss = is.na(state) | is.na(sex) | is.na(age),
         miss_deaths = sum(deaths[miss])) %>%
  ungroup() %>%
  filter(!miss) %>%
  group_by(year, state, sex, age) %>%
  summarise(deaths = sum(deaths), miss_deaths = first(miss_deaths), .groups = "drop") %>%
  group_by(year) %>%
  mutate(deaths = deaths + deaths / sum(deaths) * first(miss_deaths)) %>%
  ungroup() %>%
  select(state, sex, year, age, deaths)

# optional integer deaths
df_clean$deaths <- round(df_clean$deaths)
# complete the state-sex-year grid
deaths <- df_clean %>%
  complete(state, sex, year, age)

# population
load("data/mx_pops.rdata")
df <- mx_pops %>% 
  select(state, sex, year, age, pop) %>% 
  mutate(sex=ifelse(sex==0, "male", "female")) %>% 
  mutate(age_lower = suppressWarnings(as.numeric(sub(" to.*", "", age))),
         age_lower = ifelse(grepl("and above", age), 80, age_lower)) %>%
  arrange(state, sex, year, age_lower) %>%
  select(-age_lower)

pop <- df %>%
  group_by(state, sex, year) %>%
  mutate(
    prop_youth = sum(pop[age %in% c("15 to 19","20 to 24","25 to 29",
                                    "30 to 34","35 to 39","40 to 44",
                                    "45 to 49","50 to 54")]) /
      sum(pop)
  ) %>%
  ungroup()

# disappearances
load("data/m_disappearances.rdata")
df <- mx_disap %>% 
  select(state, sex, year, age, starts_with("total_")) %>%
  mutate(state = ifelse(state=="unkown", NA, state),
         sex = ifelse(sex==0, "male", "female"),
         age = str_replace(age, "-", " to "),
         age = if_else(age == "15 to19", "15 to 19", age)) %>% 
  arrange(state,sex,year,age)

df_clean <- df %>%
  pivot_longer(starts_with("total_"), names_to = "metric", values_to = "value") %>%
  mutate(value = replace_na(value, 0)) %>%
  group_by(sex, year, metric) %>%
  mutate(miss = is.na(state) | is.na(age),
         miss_value = sum(value[miss])) %>%
  ungroup() %>%
  filter(!miss) %>%
  group_by(sex, year, metric, state, age) %>%
  summarise(value = sum(value), miss_value = first(miss_value), .groups = "drop") %>%
  group_by(sex, year, metric) %>%
  mutate(value = value + value / sum(value) * first(miss_value)) %>%
  ungroup() %>%
  select(-miss_value) %>%
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(across(starts_with("total_"), round)) %>% 
  arrange(state, sex, year, age)
disap <- df_clean %>% 
  complete(state, sex, year, age) %>%
  mutate(disap = ifelse(total_dnll>0, total_nf/total_dnll, 0),
         age_lower = suppressWarnings(as.numeric(sub(" to.*", "", age))),
         age_lower = ifelse(grepl("and above", age), 80, age_lower)) %>%
  arrange(state, sex, year, age_lower) %>%
  select(-c(age_lower,starts_with("total_")))

# poverty
load("data/m_povertyl.rdata")
pov <- m_povertyl %>%
  group_by(state) %>%
  summarise(pov = mean(pov, na.rm = TRUE)) %>%
  ungroup()

# merge the data
df <- deaths %>% 
  filter(year>=2007) %>% 
  left_join(pop, by = c("state", "sex", "year", "age")) %>%
  left_join(disap, by = c("state", "sex", "year", "age")) %>%
  left_join(pov, by = "state") %>% 
  arrange(state, sex, year, age)

rm(list = setdiff(ls(), "df"))  
save(df, file = "data/df.rdata")
