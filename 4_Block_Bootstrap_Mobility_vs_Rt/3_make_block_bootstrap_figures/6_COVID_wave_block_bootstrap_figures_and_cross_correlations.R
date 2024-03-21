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
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient") +
  # scale_fill_viridis_c() +
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
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism == "SARS-CoV-2") %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date < as.Date("2022-02-13")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.2, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(organism == "SARS-CoV-2") %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date <= as.Date("2022-02-15")) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.9
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec_lim) +
  scale_fill_manual(values = color_vec_lim) +
  ylab("Rt") +
  xlab("Epidemic Week")
cont_rt

rt_top <- plot_grid(NULL, cont_rt, NULL, rel_widths = c(0.3, 9, 1.8), nrow = 1)
rt_mob_cont <- plot_grid(rt_top, all_mob_plot2, nrow = 2, rel_heights = c(2.5, 12), labels = "AUTO")
rt_mob_cont
# save_plot(rt_mob_cont, file = "figures/fig_5_block_bootstrap_covid_5mo_moving_window_spearman_neg_lags.png", base_width = 14, base_height = 16)

unique(all_path_ccf$mobility_metric)

covid_mob_red <- ggplot(
  all_path_ccf %>% filter(start_week >= as.Date("2020-02-16") & start_week < as.Date("2022-02-13") & pathogen == "SARS-CoV-2" &
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
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient") +
  # scale_fill_viridis_c() +
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

# rt_top <- plot_grid(NULL, cont_rt, NULL, rel_widths = c(0.25, 10, 2), nrow = 1)
cont_rt2 <- cont_rt +
  annotate("text", x = as.Date("2020-12-01"), y = 3.5, label = "Winter\n2020-2021", fontface = "bold") +
  annotate("text", x = as.Date("2021-04-15"), y = 3.5, label = "Alpha\nvariant", fontface = "bold") +
  annotate("text", x = as.Date("2021-09-01"), y = 3.5, label = "Delta\nvariant", fontface = "bold") +
  annotate("text", x = as.Date("2022-01-15"), y = 3.5, label = "Omicron BA.1\nvariant", fontface = "bold")
rt_top <- plot_grid(NULL, cont_rt2, NULL, rel_widths = c(0.3, 9, 1.6), nrow = 1)
rt_mob_cont <- plot_grid(rt_top, covid_mob_red, nrow = 2, rel_heights = c(2, 7), labels = NULL)
rt_mob_cont
save_plot(rt_mob_cont, file = "figures/fig_5_block_bootstrap_covid_5mo_moving_window_select_mobility_indicators_spearman_neg_lags.png", base_width = 14, base_height = 14)

########################################################################################
## SARS-CoV-2 waves 2020 - 2022
## statistically significant cross-correlations between Rt and mobility
########################################################################################

########################################################################################
## SAH period
########################################################################################
sc2_rt <- all_endemic_results_avg %>% filter(organism == "SARS-CoV-2")

sc2_rt %>%
  filter(date > as.Date("2020-02-15") & date < as.Date("2020-07-30")) %>%
  filter(median < 1)

ggplot(sc2_rt %>% filter(date > as.Date("2020-02-15") & date < as.Date("2020-07-30"))) +
  geom_line(aes(x = date, y = median)) +
  geom_hline(yintercept = 1, color = "red") +
  geom_vline(xintercept = as.Date("2020-05-28"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-06-14"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-07-15"), lty = "dashed") +
  theme_bw()

# SC2 date range of significant correlations during SAH: mid-to-late April to mid May 2020
all_path_ccf %>%
  filter(!(mobility_metric %in% c(
    "% not wearing masks",
    "Oxford Stringency\nIndex",
    "Within-neighborhood\nmovement",
    "Groceries and\npharmacies"
  ))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2020-07-15") & start_week > as.Date("2020-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min_date = min(start_week),
    max_date = max(start_week)
  ) %>%
  arrange(min_date)
# "2020-04-12" "2020-06-14"

## positive correlations
all_path_ccf %>%
  filter(!(mobility_metric %in% c(
    "% not wearing masks",
    "Oxford Stringency\nIndex",
    "Within-neighborhood\nmovement",
    "Groceries and\npharmacies"
  ))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-21") & start_week >= as.Date("2020-04-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week)
## most mobility indicators, late April to late May/early June 2020

all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-20")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "Child daycare"                     0.656 0.922
# 2 "Elementary and\nhigh schools"      0.764 0.898
# 3 "Colleges"                          0.651 0.827
# 4 "Restaurants"                       0.682 0.815
# 5 "Influx out-of-state\nvisitors"     0.612 0.815
# 6 "% Devices\nleaving home"           0.623 0.807
# 7 "Between-neighborhood\nmovement"    0.650 0.790
# 8 "Infux visitors\nother WA counties" 0.690 0.787
# 9 "Religious\norganizations"          0.594 0.774
# 10 "Transit"                           0.587 0.729

## range of positive correlation values
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks", "Transit", "Oxford Stringency\nIndex", "Within-neighborhood\nmovement", "Groceries and\npharmacies"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-20")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week) %>%
  pull(obs_max_ccf) %>%
  range() #  0.5942475 0.9216414

## negative correlations with OSI, mid-April to late June
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-20")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(start_week)

# range of negative correlation values
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-20")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  pull(obs_max_ccf) %>%
  range() #-0.8505591 -0.6959359

## filter to max correlation values for each mob metric
## late April to mid-May
all_path_ccf %>%
  filter(!(mobility_metric %in% c(
    "% not wearing masks",
    "Oxford Stringency\nIndex",
    "Within-neighborhood\nmovement",
    "Groceries and\npharmacies"
  ))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-20")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  slice_max(obs_max_ccf, n = 3) %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  print(n = 30)

########################################################################################
## winter wave, 2020 - 2021
########################################################################################
ggplot(sc2_rt %>% filter(date > as.Date("2020-09-01") & date < as.Date("2021-02-01"))) +
  geom_line(aes(x = date, y = median)) +
  geom_hline(yintercept = 1, color = "red") +
  geom_vline(xintercept = as.Date("2020-10-05"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-10-28"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-10"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-10"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-03"), lty = "dashed") +
  theme_bw()

## positive correlations
all_path_ccf %>%
  filter(!(mobility_metric %in% c("Elementary and\nhigh schools", "Colleges", "% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week >= as.Date("2020-10-01") & start_week < as.Date("2021-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week)

## early increases
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(!(mobility_metric %in% c("Elementary and\nhigh schools", "Colleges", "% not wearing masks"))) %>%
  filter(start_week >= as.Date("2020-10-01") & start_week < as.Date("2020-11-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric, pathogen) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4)
# mobility_metric              pathogen       n
# 1 "Religious\norganizations" SARS-CoV-2     6
# 2 "Restaurants"              SARS-CoV-2     5

all_path_ccf %>%
  filter(pathogen == "SARS-CoV-2") %>%
  filter(mobility_metric %in% c("Restaurants", "Religious\norganizations")) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week >= as.Date("2020-10-01") & start_week < as.Date("2020-12-10")) %>%
  distinct() %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric             min   max
# <fct>                      <dbl> <dbl>
# 2 "Religious\norganizations" 0.636 0.756
# 3 "Restaurants"              0.569 0.748

all_path_ccf %>%
  filter(mobility_metric %in% c("Infux visitors\nother WA counties", "Influx out-of-state\nvisitors", "% not wearing masks")) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-11-15") & start_week < as.Date("2021-02-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "Influx out-of-state\nvisitors"     0.687 0.828
# 2 "Infux visitors\nother WA counties" 0.653 0.779
# 3 "% not wearing masks"               0.589 0.719

sc2_keep <- all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
           mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2021-03-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric, pathogen) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 3) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

sc2_keep

all_path_ccf %>%
  filter(pathogen == "SARS-CoV-2") %>%
  filter(mobility_metric %in% sc2_keep) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2021-03-01"))

all_path_ccf %>%
  filter(pathogen == "SARS-CoV-2") %>%
  # filter(mobility_metric %in% sc2_keep) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2021-03-01")) %>%
  distinct() %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)

## negative correlations (OSI)
all_path_ccf %>%
  filter(pathogen == "SARS-CoV-2") %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  arrange(start_week)

all_path_ccf %>%
  filter(pathogen == "SARS-CoV-2") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2021-02-01")) %>%
  distinct() %>%
  filter(mobility_metric %in% c("Infux visitors\nother WA counties", " Influx out-of-state\nvisitors", "% Devices\nleaving home")) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)

all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week >= as.Date("2020-11-01") & start_week < as.Date("2020-12-27")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week) %>%
  arrange(obs_max_ccf) %>%
  pull(obs_max_ccf) %>%
  range() #-0.8493392 -0.6147512

########################################################################################
## alpha wave, spring 2021
########################################################################################
ggplot(sc2_rt %>% filter(date > as.Date("2021-02-15") & date < as.Date("2021-04-21"))) +
  geom_line(aes(x = date, y = median)) +
  geom_hline(yintercept = 1, color = "red") +
  geom_vline(xintercept = as.Date("2021-03-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-04-08"), lty = "dashed") +
  theme_bw()

## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2021-04-30") & start_week > as.Date("2021-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)
## only groceries and pharmacies for 3 weeks; visitor inflow 1 week

# negative correlations
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2021-04-21") & start_week > as.Date("2021-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)
## zero weeks

########################################################################################
## delta wave, summer 2021
########################################################################################
ggplot(sc2_rt %>% filter(date > as.Date("2021-06-01") & date < as.Date("2021-08-21"))) +
  geom_line(aes(x = date, y = median)) +
  geom_hline(yintercept = 1, color = "red") +
  geom_vline(xintercept = as.Date("2021-06-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-25"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-08-15"), lty = "dashed") +
  theme_bw()

## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
           mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric, pathogen) %>%
  tally() %>%
  arrange(-n)
# mobility_metric                     pathogen       n
# # <fct>                               <fct>      <int>
# 1 "Infux visitors\nother WA counties" SARS-CoV-2    10
# 2 "Influx out-of-state\nvisitors"     SARS-CoV-2    10
# 3 "% Devices\nleaving home"           SARS-CoV-2    10
# 4 "Restaurants"                       SARS-CoV-2    10
# 5 "% not wearing masks"               SARS-CoV-2    10

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
           mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-01") & start_week > as.Date("2021-06-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(min)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "Restaurants"                       0.661 0.852
# 2 "% not wearing masks"               0.675 0.943
# 3 "% Devices\nleaving home"           0.733 0.936
# 4 "Influx out-of-state\nvisitors"     0.794 0.928
# 5 "Infux visitors\nother WA counties" 0.820 0.889

## negative correlations
all_path_ccf %>%
  filter(sig == "yes" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-05-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)
## OSI early Aug

########################################################################################
## First omicron wave, late 2021 - early 2022
########################################################################################
ggplot(sc2_rt %>% filter(date > as.Date("2021-11-01") & date < as.Date("2022-01-15"))) +
  geom_line(aes(x = date, y = median)) +
  geom_hline(yintercept = 1, color = "red") +
  geom_vline(xintercept = as.Date("2021-11-20"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-12-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2022-01-10"), lty = "dashed") +
  theme_bw()

## check number of statistically significant windows
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "SARS-CoV-2") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4) %>%
  arrange(pathogen, -n) %>%
  print(n = 60)
# 1 SARS-CoV-2 "Groceries and\npharmacies"            13
# 2 SARS-CoV-2 "Within-neighborhood\nmovement"         9
# 3 SARS-CoV-2 "Religious\norganizations"              9
# 4 SARS-CoV-2 "Transit"                               6
# 5 SARS-CoV-2 "Infux visitors\nother WA counties"     5
# 6 SARS-CoV-2 "Restaurants"                           5

## initial increases

sc2_keep <- all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
           mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2022-01-01") & start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric, pathogen) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 3) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

sc2_keep

## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric %in% sc2_keep) %>%
  filter(start_week < as.Date("2021-12-20") & start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric %in% sc2_keep) %>%
  filter(start_week < as.Date("2021-12-20") & start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(min)

## increase and decline
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2022-02-20") & start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, start_week)

## negative correlations with OSI
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0 &
           mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-12-20") & start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)
## none

## length of lags between mobility and Rt
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "SARS-CoV-2") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(mobility_metric, -n) %>%
  print(n = 50)
# mobility_metric                     obs_max_ccf_lag     n
# <fct>                                         <dbl> <int>
# 1 "Within-neighborhood\nmovement"                  -4     9
# 2 "Infux visitors\nother WA counties"              -4     3
# 3 "Infux visitors\nother WA counties"              -3     2
# 4 "Influx out-of-state\nvisitors"                  -1     3
# 5 "Restaurants"                                    -4     5
# 6 "Groceries and\npharmacies"                      -1     7
# 7 "Groceries and\npharmacies"                      -2     6
# 8 "Transit"                                        -4     6
# 9 "Religious\norganizations"                       -3     5
# 10 "Religious\norganizations"                       -4     4
# 11 "Child daycare"                                  -4     2
# 12 "Elementary and\nhigh schools"                   -4     1
