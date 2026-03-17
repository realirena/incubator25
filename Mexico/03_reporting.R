#------------------------------------------------------------------------------#
# Paper: Under-reported homicide count in Mexico
# Title: Program from constructing informative beta prior 
# Data:
# Author:
# Date:
#------------------------------------------------------------------------------#
Sys.setenv(LANGUAGE = "en")
pacman::p_load(tidyverse, readxl)
# GBD
gbd <- read.csv("other/gbd.csv")
gbd <- gbd[c(8,14:18)]
gbd <- gbd %>% 
  rename_with(~ sub("_name$", "", .x)) %>% 
  filter(metric == "Rate") %>% 
  select(sex, year, val) %>% 
  rename("gbd" = "val") %>% 
  arrange(sex, year)

# WHO Mortality Database
whomd <- read.csv(text = sub(",$", "", readLines("other/whomd.csv")), skip = 6)
whomd <- whomd[c(5:6, 12)]
colnames(whomd) <- c("year", "sex", "whomd")
whomd <- whomd %>% 
  select(sex, year, whomd) %>% 
  filter(sex!="All") %>% 
  arrange(sex, year)
# GHE
ghe <- read_xlsx("other/ghe.xlsx")
colnames(ghe) <- ghe[2,]
ghe <- ghe[c(10,13,24,27,29)]
colnames(ghe) <- c("year", "sex", "ghe", "lower", "upper")
ghe <- ghe %>% 
  filter(sex %in% c("Male", "Female")) %>% 
  select(sex, year, ghe) %>% 
  mutate(across(c(year, ghe), as.numeric)) %>% 
  arrange(sex, year)
# UNDOC
unodc <- read_xlsx("other/unodc.xlsx")
colnames(unodc) <- unodc[2,]
unodc <- unodc %>% filter(Country =="Mexico" &
                            Indicator =="Victims of intentional homicide" &
                            Dimension == "Total" &
                            Category == "Total" &
                            Age == "Total" &
                            `Unit of measurement`=="Rate per 100,000 population")
unodc <- unodc[c(8,10,12)]
colnames(unodc) <- c("sex","year","unodc")
unodc <- unodc %>% 
  filter(sex!="Total") %>% 
  mutate(across(c(year, unodc), as.numeric)) %>% 
  arrange(sex, year)
# INEGI
inegi <- read.csv("other/inegi.csv")
inegi <- inegi %>% filter(STATE_ID==0)
inegi <- inegi[c(1,11,12)]
names(inegi) <- tolower(names(inegi))
inegi <- inegi %>% 
  pivot_longer(cols = c(male_rate, female_rate), 
               names_to = "sex", names_pattern = "(.*)_rate",
               values_to = "inegi") %>% 
  mutate(sex = str_to_sentence(sex)) %>% 
  arrange(sex, year)

# merge
data_list <- list(inegi, gbd, ghe, whomd, unodc)
df1 <- reduce(data_list, full_join, by = c("sex", "year"))
df1 <- df1 %>% 
  arrange(sex, year)
save(df1, file = "other/df1.RData")

load("other/df1.Rdata")
dat <- df %>% 
  group_by(year,sex) %>% 
  summarise(deaths = sum(deaths, na.rm = TRUE),
            pop = sum(pop, na.rm = TRUE), 
            .groups = "drop") %>% 
  mutate(rate = deaths/pop*100000,
         sex = str_to_sentence(sex)) %>% 
  select(sex, year, everything()) %>% 
  arrange(sex,year) %>% 
  left_join(df1, by = c("sex", "year")) %>% 
  pivot_longer(
    cols = c(inegi, gbd, ghe, whomd, unodc),
    names_to = "source", values_to = "rate1"
  ) %>% 
  mutate(rp=rate/rate1)

# moment-matching beta prior specification (Jersakova et al. 2022)
beta_mom <- function(p, eps = 1e-6) {
  
  p <- pmin(pmax(p, eps), 1 - eps)
  
  m <- mean(p, na.rm = TRUE)
  v <- var(p, na.rm = TRUE)
  
  # sanity checks
  if (v <= 0) stop("Variance is zero or negative; can't fit Beta by moments.")
  if (v >= m * (1 - m)) stop("Variance too large for a Beta distribution.")
  
  a <- (m^2 * (1 - m)) / v - m
  b <- a * (1 - m) / m
  c(a = a, b = b)
}
rp <- dat$rp[dat$sex=="Female"]
rp <- rp[is.finite(rp) & rp> 0 & rp < 1]
ab <- beta_mom(rp)
ab

# beta_elicit function (Chen et al. 2022)
rp <- dat$rp[dat$sex=="Male"]
rp <- rp[is.finite(rp) & rp> 0 & rp < 1]

# mode
d <- density(rp, from = 0, to = 1)
mode_hat <- d$x[which.max(d$y)]
mode_hat

# cutoff
c <- as.numeric(quantile(rp, 0.95))
c
alpha <- 0.001

# beta prior
beta_elicit <- function(mode, c, alpha){
  b <- 1 # initial value
  a <- (1+mode*(b-2))/(1-mode)
  f <- function(b){
    a <- (1+mode*(b-2))/(1-mode)
    qbeta(1-alpha,a,b)-c
  }
  x <- uniroot(f, lower = 1, upper = 100)
  b <- x$root
  a <- (1+mode*(b-2))/(1-mode)
  list(a,b)
}
ab <- beta_elicit(mode = mode_hat, c = c, alpha = alpha)
ab
