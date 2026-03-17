pacman::p_load(tidyverse, ggplot2, viridis, geofacet)
decomp_summary <- decomp_sum %>%
  filter(sex == "male", cause == "mx_h") %>%
  group_by(state, age) %>%
  summarise(
    mean = mean(contribution),
    lower = quantile(contribution, 0.025),
    upper = quantile(contribution, 0.975),
    .groups = "drop"
  ) %>%
  # filter(age<=70) %>%
  mutate(age = factor(age, levels = sort(unique(age))),
         state=str_to_title(state),
         state=recode(state,
                      "Mexico City"="Ciudad de México",
                      "Mexico State"="México",
                      "Michoacan"="Michoacán",
                      "Nuevo Leon"="Nuevo León",
                      "Queretaro"="Querétaro",
                      "San Luis Potosi"="San Luis Potosí",
                      "Yucatan"="Yucatán"))

ggplot(decomp_summary,aes(x = age, y = mean, fill = age, colour = age)) +
  geom_col(width = 0.7, alpha = 0.5, linewidth = 0.2) +
  geom_linerange(aes(ymin = lower, ymax = upper),linewidth = 0.2) +
  scale_x_discrete(breaks = seq(0,75,15))+
  scale_y_continuous(breaks = c(0.00, -0.05, -0.09)#n.breaks = 4
    )+
  scale_fill_viridis_d(option = "C") +
  scale_colour_viridis_d(option = "C", end = 0.8) +
  geom_hline(yintercept = 0, lwd = 0.3) +
  facet_geo(~state,grid="mx_state_grid2",
            labeller=as_labeller(\(x) stringr::str_wrap(x,width=14)))+
  labs(x = "Age", y = "Contribution") +
  theme_minimal() +
  theme(
    axis.title.x =  element_text(size = 12, face = "bold", 
                                 margin = margin(b = 2)),
    axis.title.y = element_text(size = 12, face = "bold", 
                                margin = margin(t=2)),
    axis.text.x = element_text(size = 10, angle = 90, 
                               hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
    legend.position = "none")
ggsave("figure/2_decomp_m.png", height = 22, width = 26, units = "cm", dpi = 300)
ggsave("figure/2_decomp_m.pdf", height = 22, width = 26, units = "cm", dpi = 300)
