#------------------------------------------------------------------------------#
# Title:
# Paper:
# Data:
# Author:
# Date:
#------------------------------------------------------------------------------#

# load the required libraries
library(dplyr)
library(srtingr)
library(stringi)

# load the death data
load("df.Rdata")

#load the reporting estimates from Garguilo et al. (2025)
rp <- read.csv(gzfile("data/CO-underreporting.csv.gz"))
rp <- rp %>% rename("state" = "name") %>%
  dplyr::select(replicate, state, period3y, sample_num, prop_completeness) %>%
  mutate(
    year = str_extract_all(period3y, "\\d{4}") %>%
      lapply(as.numeric) %>%
      sapply(function(x) mean(x)),
    year = round(year)
  ) %>%
  relocate(year, .after = state) %>%
  mutate(
    state = tolower(state),
    state = stringi::stri_trans_general(state, "Latin-ASCII")
  )

rp <- rp %>%
  mutate(state = case_when(
    state == "atla!ntico" ~ "atlantico",
    state == "bogota!" ~ "bogota",
    state == "bola-var" ~ "bolivar",
    state == "boyaca!" ~ "boyaca",
    state == "caqueta!" ~ "caqueta",
    state == "ca³rdoba" ~ "cordoba",
    state == "choca³" ~ "choco",
    state == "naria+/-o" ~ "narino",
    state == "n. de santander" ~ "norte de santander",
    state == "quinda-o" ~ "quindio",
    state == "san andra(C)s" ~ "san andres",
    state == "guaina-a" ~ "guainia",
    state == "vaupa(C)s" ~ "vaupes",
    TRUE ~ state
  ))
rp$state[rp$state == "valle del cauca"] <- "valle"
rp$state <- gsub("ca³rdoba", "cordoba", rp$state)
rp$state <- gsub("choca³", "choco", rp$state)

setdiff(unique(rp$state), unique(df$state[df$country=="colombia"]))

rm(list = setdiff(ls(), c("df", "rp")))