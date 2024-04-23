########################################################################################
## Moving window Spearman cross-correlations between Rt and mobility, Individual COVID-19 waves
## Figure 5: SARS-CoV-2 moving window cross-correlations
## Significant cross-correlations between SARS-CoV-2 Rt and mobility during each wave
########################################################################################

library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
library(gghighlight)
library(cdcfluview)

########################################################################################
## Import 5mo block bootstrap results
########################################################################################
all_path_ccf <- read_rds("4_Block_Bootstrap_Mobility_vs_Rt/3_make_block_bootstrap_figures/all_5mo_block_bootstrap_results.rds")

all_endemic_results_avg <- read_csv("2_Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")
all_endemic_results_avg$organism <- as.factor(all_endemic_results_avg$organism)
levels(all_endemic_results_avg$organism)

levels(all_endemic_results_avg$organism) <- c(
  "Adenovirus",
  "Enterovirus",
  "Human metapneumovirus",
  "Human parainfluenza 1 + 2",
  "Human parainfluenza 3 + 4",
  "Influenza A/H1N1",
  "Influenza A/H3N2",
  "Influenza B",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Rhinovirus",
  "SARS-CoV-2",
  "Seasonal CoV 229E + OC43",
  "Seasonal CoV HKU1 + NL63"
)

levels(all_endemic_results_avg$organism)

unique(all_endemic_results_avg$organism)
all_endemic_results_avg$organism <- factor(all_endemic_results_avg$organism,
  levels = c(
    "Influenza A/H1N1",
    "Influenza A/H3N2",
    "Influenza B",
    "Respiratory syncytial virus (RSV) A",
    "Respiratory syncytial virus (RSV) B",
    "Human metapneumovirus",
    "Seasonal CoV 229E + OC43",
    "Seasonal CoV HKU1 + NL63",
    "Human parainfluenza 1 + 2",
    "Human parainfluenza 3 + 4",
    "Rhinovirus",
    "Enterovirus",
    "Adenovirus",
    "SARS-CoV-2"
  )
)
levels(all_endemic_results_avg$organism)

all_endemic_results_avg <- bind_cols(all_endemic_results_avg, mmwr_week(all_endemic_results_avg$date)[, 1:2])
all_endemic_results_avg$epi_date <- mmwr_week_to_date(all_endemic_results_avg$mmwr_year, all_endemic_results_avg$mmwr_week)
head(all_endemic_results_avg)
unique(all_endemic_results_avg$organism)

rt_weekly <- all_endemic_results_avg %>%
  group_by(epi_date, tag, level, organism) %>%
  summarize_at(c("lower", "upper", "median"), ~ mean(.x, na.rm = T)) %>%
  arrange(organism, epi_date) %>%
  filter(organism != "Influenza A/H3N2") %>%
  droplevels()

unique(rt_weekly$organism)
levels(all_path_ccf$pathogen)
levels(rt_weekly$organism)
levels(all_path_ccf$pathogen)
levels(rt_weekly$organism) <- levels(all_path_ccf$pathogen)

########################################################################################
######## Figure 5: Cross-correlations between Rt and mobility across SARS-CoV-2 waves
########################################################################################
all_mob_plot2 <- ggplot(
  all_path_ccf %>%
    filter(obs_max_ccf_lag < 1) %>%
    filter(start_week >= as.Date("2020-02-16") & start_week < as.Date("2022-02-13") & pathogen == "SARS-CoV-2"),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  annotate("rect",
    xmin = as.Date("2020-03-23"),
    xmax = as.Date("2020-06-05"),
    ymin = -Inf,
    ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  annotate("rect",
    xmin = as.Date("2020-09-25"),
    xmax = as.Date("2021-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-03-01"),
    xmax = as.Date("2021-05-31"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-07-01"),
    xmax = as.Date("2021-11-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-12-01"),
    xmax = as.Date("2022-03-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_point(pch = 21, size = 3, stroke = 0.2, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 2, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 7) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )

all_mob_plot2

color_vec_lim <- c("#ABB065")

cont_rt <- ggplot() +
  annotate("rect",
    xmin = as.Date("2020-03-23"),
    xmax = as.Date("2020-06-05"),
    ymin = -Inf,
    ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  annotate("rect",
    xmin = as.Date("2020-09-25"),
    xmax = as.Date("2021-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-03-01"),
    xmax = as.Date("2021-05-31"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-07-01"),
    xmax = as.Date("2021-11-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-12-01"),
    xmax = as.Date("2022-03-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 0.5, alpha = 0.5) +
  geom_line(
    data = rt_weekly %>%
      filter(organism == "SARS-CoV-2") %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date <= as.Date("2022-02-13")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 0.8, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(organism == "SARS-CoV-2") %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date <= as.Date("2022-02-13")) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.5
  ) +
  theme_bw(base_size = 7) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02), limits = c(as.Date("2020-02-16"), as.Date("2022-03-01"))) +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "none",
    strip.text.y.right = element_text(angle = 0, size = 6)
  ) +
  scale_color_manual(values = color_vec_lim) +
  scale_fill_manual(values = color_vec_lim) +
  ylab("Rt") +
  xlab("Epidemic Week")
cont_rt

rt_top <- plot_grid(NULL, cont_rt, NULL, rel_widths = c(0.3, 9, 1.8), nrow = 1)
rt_mob_cont <- plot_grid(rt_top, all_mob_plot2, nrow = 2, rel_heights = c(2.5, 12), labels = "AUTO", label_size = 7)
rt_mob_cont

unique(all_path_ccf$mobility_metric)

all_path_ccf %>%
  filter(start_week >= as.Date("2020-02-16") & start_week < as.Date("2022-02-13") &
    pathogen == "SARS-CoV-2" &
    mobility_metric %in% c(
      "Within-neighborhood\nmovement",
      # "Between-neighborhood\nmovement",
      "Infux visitors\nother WA counties",
      "Influx out-of-state\nvisitors",
      "% Devices\nleaving home",
      "Religious\norganizations",
      "Restaurants",
      "Transit",
      "Groceries and\npharmacies",
      "% not wearing masks",
      "Oxford Stringency\nIndex"
    )) %>%
  pull(start_week) %>%
  range()

covid_mob_red <- ggplot(
  all_path_ccf %>% filter(start_week >= as.Date("2020-02-16") & start_week <= as.Date("2022-02-13") & pathogen == "SARS-CoV-2" &
    mobility_metric %in% c(
      "Within-neighborhood\nmovement",
      "Infux visitors\nother WA counties",
      "Influx out-of-state\nvisitors",
      "% Devices\nleaving home",
      "Religious\norganizations",
      "Restaurants",
      "Transit",
      "Groceries and\npharmacies",
      "% not wearing masks",
      "Oxford Stringency\nIndex"
    )),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  annotate("rect",
    xmin = as.Date("2020-03-23"),
    xmax = as.Date("2020-06-05"),
    ymin = -Inf,
    ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  annotate("rect",
    xmin = as.Date("2020-09-25"),
    xmax = as.Date("2021-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-03-01"),
    xmax = as.Date("2021-05-31"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-07-01"),
    xmax = as.Date("2021-11-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  annotate("rect",
    xmin = as.Date("2021-12-01"),
    xmax = as.Date("2022-03-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_point(pch = 21, size = 3, stroke = 0.2, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 2, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed", size = 0.5, alpha = 0.5) +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 7) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0, size = 6),
    axis.text.y = element_text(size = 5)
  )

## add wave labels
cont_rt2 <- cont_rt +
  annotate("text", x = as.Date("2020-12-01"), y = 3.5, size = 2, label = "Winter\n2020-2021", fontface = "bold") +
  annotate("text", x = as.Date("2021-04-15"), y = 3.5, size = 2, label = "Alpha\nvariant", fontface = "bold") +
  annotate("text", x = as.Date("2021-09-01"), y = 3.5, size = 2, label = "Delta\nvariant", fontface = "bold") +
  annotate("text", x = as.Date("2022-01-15"), y = 3.5, size = 2, label = "Omicron BA.1\nvariant", fontface = "bold")

## combine rt and mobility plots
rt_top <- plot_grid(NULL, cont_rt2, NULL, rel_widths = c(0.15, 9, 1.3), nrow = 1)
rt_mob_cont <- plot_grid(rt_top, covid_mob_red, nrow = 2, rel_heights = c(2, 10), labels = NULL)
rt_mob_cont

save_plot(rt_mob_cont,
  file = "figures/fig_5_block_bootstrap_covid_5mo_moving_window_select_mobility_indicators_neg_lags.pdf",
  units = "mm", base_width = 180, base_height = 160, dpi = 300
)
