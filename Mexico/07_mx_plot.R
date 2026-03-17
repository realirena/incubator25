# mortality plot
pacman::p_load(tidyverse, ggplot2, RColorBrewer, viridis, geofacet)
load("results/post_draws2.RData")

# plot for overall death rates by year
df1 <- post_draws %>%
  filter(sex == "male") %>%
  group_by(state, year, .draw) %>%
  summarise(
    dx = sum(deaths, na.rm = TRUE),
    dx_true = round(sum(lambda_true)),
    pop = sum(pop),
    # mx = log(sum(deaths, na.rm = TRUE) / sum(pop)),
    # mx_true = log(sum(lambda_true) / sum(pop)),
    .groups = "drop"
  ) %>%
  mutate(
    mx = log(dx/pop),
    mx_true = log(dx_true/pop),
  ) %>% 
  mutate(state = str_to_title(state)) %>%
  mutate(state = recode(state,
                        "Mexico City" = "Ciudad de México",
                        "Mexico State" = "México",
                        "Michoacan" = "Michoacán",
                        "Nuevo Leon" = "Nuevo León",
                        "Queretaro" = "Querétaro",
                        "San Luis Potosi" = "San Luis Potosí",
                        "Yucatan" = "Yucatán"
  ))
display.brewer.all(colorblindFriendly = TRUE)
col1 <- brewer.pal(9, "Reds")[5]
col2 <- brewer.pal(9, "Reds")[9]

ggplot(df1, aes(x = year)) +
  geom_line(
    aes(y = mx_true, group = .draw, colour = "mu"),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_line(
    aes(y = mx, colour = "muhat"),
    linewidth = 0.3
  ) +
  facet_geo(~state, grid = mx_state_grid2,
            labeller=as_labeller(function(x) str_wrap(x,width=14))) +
  scale_x_continuous(breaks = c(2010, 2015, 2020),labels = c("2010", "'15", "'20"))+
  scale_y_continuous(n.breaks = 4)+
  scale_colour_manual(
    values = c(mu = col1, muhat = col2),
    labels = c(
      mu = expression("True Rate (" * mu * ")"),
      muhat = expression("Observed Rate (" * hat(mu) * ")")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
  labs(x = "Year", y = "Log Death Rate", colour = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = margin(.2,.2,.2,.2,"cm")
  )
ggsave("figure/1_mx_m.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/1_mx_m.pdf", height = 22, width = 26, units = "cm", dpi = 300)

# plot deaths
ggplot(df1, aes(x = year)) +
  geom_line(
    aes(y = dx_true, group = .draw, colour = "mu"),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_line(
    aes(y = dx, colour = "muhat"),
    linewidth = 0.3
  ) +
  facet_geo(~state, grid = mx_state_grid2,
            labeller=as_labeller(function(x) str_wrap(x,width=14))) +
  scale_x_continuous(breaks = c(2010, 2015, 2020),labels = c("2010", "'15", "'20"))+
  scale_y_continuous(n.breaks = 4)+
  scale_colour_manual(
    values = c(mu = col1, muhat = col2),
    labels = c(
      mu = expression("True counts (" * D * ")"),
      muhat = expression("Observed counts (y)")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
  labs(x = "Year", y = "Homicide counts", colour = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = margin(.2,.2,.2,.2,"cm")
  )
ggsave("figure/1_dx_m.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/1_dx_m.pdf", height = 22, width = 26, units = "cm", dpi = 300)

# plot by broad age-groups
df2 <- post_draws %>% 
  filter(sex=="female") %>% 
  mutate(age_group = case_when(
    age < 15 ~ "<15",
    age >= 15 & age < 50 ~ "15-49",
    age >= 50 ~ ">=50"
  )) %>% 
  group_by(state, sex, year, .draw, age_group) %>% 
  summarise(
    mx = log(sum(deaths, na.rm = TRUE) / sum(pop)),
    mx_true = log(sum(lambda_true) / sum(pop)),
    .groups = "drop"
  ) %>% 
  mutate(state = str_to_title(state)) %>%
  mutate(state = recode(state,
                        "Mexico City" = "Ciudad de México",
                        "Mexico State" = "México",
                        "Michoacan" = "Michoacán",
                        "Nuevo Leon" = "Nuevo León",
                        "Queretaro" = "Querétaro",
                        "San Luis Potosi" = "San Luis Potosí",
                        "Yucatan" = "Yucatán")) %>% 
  select(state, sex, year, age_group, .draw, everything()) %>% 
  arrange(state, sex, year, age_group, .draw)

ggplot(subset(df2, age_group == ">=50"), aes(x = year)) +
  geom_line(
    aes(y = mx_true, group = .draw, colour = "mu"),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_line(
    aes(y = mx, colour = "muhat"),
    linewidth = 0.3
  ) +
  facet_geo(~state, grid = mx_state_grid2,
            labeller=as_labeller(function(x) str_wrap(x,width=14))) +
  scale_x_continuous(breaks = c(2010, 2015, 2020),labels = c("2010", "'15", "'20"))+
  scale_y_continuous(n.breaks = 4)+
  scale_colour_manual(
    values = c(mu = col1, muhat = col2),
    labels = c(
      mu = expression("True Rate (" * mu * ")"),
      muhat = expression("Observed Rate (" * hat(mu) * ")")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
  labs(x = "Year", y = "Log Death Rate", colour = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = margin(.2,.2,.2,.2,"cm")
  )
ggsave("figure/1_mx80_f.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/1_mx80_f.pdf", height = 22, width = 26, units = "cm", dpi = 300)

# plot age pattern over periods
# period: 2007-11, 2011-14, 2015-2018, 2018-2023

df3 <- post_draws %>%
  filter(sex=="female") %>% 
  mutate(period = case_when(
    year >= 2007 & year <= 2010 ~ "2007-10",
    year >= 2011 & year <= 2014 ~ "2011-14",
    year >= 2015 & year <= 2018 ~ "2015-18",
    year >= 2019 & year <= 2023 ~ "2019-23"
  )) %>% 
  group_by(state, age, .draw, period) %>% 
  summarise(
    mx = log(sum(deaths, na.rm = TRUE) / sum(pop)),
    mx_true = log(sum(lambda_true) / sum(pop)),
    .groups = "drop"
  ) %>% 
  mutate(state = str_to_title(state)) %>%
  mutate(state = recode(state,
                        "Mexico City" = "Ciudad de México",
                        "Mexico State" = "México",
                        "Michoacan" = "Michoacán",
                        "Nuevo Leon" = "Nuevo León",
                        "Queretaro" = "Querétaro",
                        "San Luis Potosi" = "San Luis Potosí",
                        "Yucatan" = "Yucatán")) %>% 
  select(state, period, age, .draw, everything()) %>% 
  arrange(state, period, age, .draw)

display.brewer.all(colorblindFriendly = TRUE)
col1 <- brewer.pal(9, "Reds")[5]
col2 <- brewer.pal(9, "Reds")[9]

ggplot(subset(df3, period=="2015-18"), aes(x = age)) +
  geom_line(
    aes(y = mx_true, group = .draw, colour = "mu"),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_line(
    aes(y = mx, colour = "muhat"),
    linewidth = 0.3
  ) +
  facet_geo(~state, grid = mx_state_grid2,
            labeller=as_labeller(function(x) str_wrap(x,width=14))) +
  scale_x_continuous(breaks = seq(0,80,20))+
  scale_y_continuous()+
  scale_colour_manual(
    values = c(mu = col1, muhat = col2),
    labels = c(
      mu = expression("True Rate (" * mu * ")"),
      muhat = expression("Observed Rate (" * hat(mu) * ")")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
  labs(x = "Age", y = "Log Death Rate", colour = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = margin(.2,.2,.2,.2,"cm")
  )
ggsave("figure/1_mx2015_18_f.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/1_mx2015_18_f.pdf", height = 22, width = 26, units = "cm", dpi = 300)

# age pattern
df3 <- post_draws %>% 
  filter(sex=="female") %>% 
  mutate(mx = log(deaths/pop),
         mx_true = log(lambda_true/pop)) %>% 
  mutate(state = str_to_title(state)) %>%
  mutate(state = recode(state,
                        "Mexico City" = "Ciudad de México",
                        "Mexico State" = "México",
                        "Michoacan" = "Michoacán",
                        "Nuevo Leon" = "Nuevo León",
                        "Queretaro" = "Querétaro",
                        "San Luis Potosi" = "San Luis Potosí",
                        "Yucatan" = "Yucatán"))
ggplot(subset(df3, year==2018), aes(x = age)) +
  geom_line(
    aes(y = mx_true, group = .draw, colour = "mu"),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_line(
    aes(y = mx, colour = "muhat"),
    linewidth = 0.3
  ) +
  facet_geo(~state, grid = mx_state_grid2,
            labeller=as_labeller(function(x) str_wrap(x,width=14))) +
  scale_x_continuous(breaks = seq(0,80,20))+
  scale_y_continuous()+
  scale_colour_manual(
    values = c(mu = col1, muhat = col2),
    labels = c(
      mu = expression("True Rate (" * mu * ")"),
      muhat = expression("Observed Rate (" * hat(mu) * ")")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
  labs(x = "Age", y = "Log Death Rate", colour = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.margin = margin(.2,.2,.2,.2,"cm")
  )
ggsave("figure/1_mx2018.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/1_mx2018.pdf", height = 22, width = 26, units = "cm", dpi = 300)
