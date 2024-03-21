########################################################################################
## Moving window Spearman cross-correlations between Rt and mobility, Fall 2019 to Winter 2020
## 1. Figure S8. Season 2019-2020 moving window cross-correlations (pre-pandemic)
## 2. Significant cross-correlations in Fall 2019
## 3. Figure S11: late 2019 - early 2020 moving window cross-correlations (early pandemic)
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

keep <- c(
  "Influenza A/H1N1", "Influenza B",
  "RSV A",
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  "hPIV 1 + 2",
  "hPIV 3 + 4",
  "Adenovirus",
  "Enterovirus",
  "Rhinovirus",
  "SARS-CoV-2"
)

color_vec <- c(
  "#332288",
  # "#6699CC",
  "#88CCEE",
  "#44AA99",
  "#117733",
  "#999933",
  "#DDCC77",
  "#661100",
  "#CC6677",
  "#AA4466",
  "#882255",
  "#993377",
  "#AA4499",
  "#ABB065"
)
########################################################################################
## Figure S8. Season 2019-2020 moving window cross-correlations (pre-pandemic)
########################################################################################

## Rt plot
fall_2019_rt <- ggplot() +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism %in% keep) %>%
      filter(epi_date > as.Date("2019-08-20") &
               epi_date < as.Date("2020-02-01")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(epi_date > as.Date("2019-08-20") &
               epi_date < as.Date("2020-02-01")) %>%
      filter(organism %in% keep) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.8
  ) +
  theme_bw(base_size = 12) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec) +
  facet_wrap(~organism, nrow = 1) +
  ylab("Rt") +
  xlab("Epidemic Week")
fall_2019_rt

## cross-correlation plot
fall_2019_plot <- ggplot(
  all_path_ccf %>% filter(obs_max_ccf_lag < 1) %>%
    filter(start_week > as.Date("2019-08-20") & start_week < as.Date("2020-02-01") &
             pathogen != "SARS-CoV-2" &
             !(mobility_metric %in% c("% not wearing masks", "Hospitals", "Oxford Stringency\nIndex"))) %>%
    droplevels(),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.01") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )
fall_2019_plot

# combined Rt and cross-correlations
rt_top <- plot_grid(NULL, fall_2019_rt, NULL, rel_widths = c(0.2, 9, 1), nrow = 1)
rt_mob_fall_2019 <- plot_grid(rt_top, fall_2019_plot, nrow = 2, rel_heights = c(1.5, 12), labels = "AUTO")
rt_mob_fall_2019
unique(all_path_ccf$mobility_metric)

### reduce number of mobility indicators to improve visibility
fall_2019_plot_red <- ggplot(
  all_path_ccf %>% filter(obs_max_ccf_lag < 1) %>%
    filter(start_week > as.Date("2019-08-20") & start_week < as.Date("2020-01-01") &
             pathogen != "SARS-CoV-2" &
             mobility_metric %in% c(
               "% Devices\nleaving home",
               "Between-neighborhood\nmovement",
               # "Influx out-of-state\nvisitors",
               "Religious\norganizations", "Child daycare",
               "Elementary and\nhigh schools", "Colleges"
             )) %>%
    droplevels(),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed", lwd = 1.2) +
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )
fall_2019_plot_red

fall_2019_rt_red <- ggplot() +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism %in% keep) %>%
      filter(epi_date > as.Date("2019-08-20") &
               epi_date < as.Date("2020-01-01")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(epi_date > as.Date("2019-08-20") &
               epi_date < as.Date("2020-01-01")) %>%
      filter(organism %in% keep) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.8
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec) +
  facet_wrap(~organism, nrow = 1) +
  ylab("Rt") +
  xlab("Epidemic Week")
fall_2019_rt_red

# combined Rt and cross-correlations
rt_top <- plot_grid(NULL, fall_2019_rt_red, NULL, rel_widths = c(0.12, 7.5, 0.8), nrow = 1)
rt_mob_fall_2019_red <- plot_grid(rt_top, fall_2019_plot_red, nrow = 2, rel_heights = c(1.5, 5), labels = "AUTO")
rt_mob_fall_2019_red
save_plot(rt_mob_fall_2019_red,
          file = "figures/fig_s8_block_bootstrap_fall_2019_select_mob_metrics_spearman_neg_lags.png",
          base_width = 22, base_height = 10
)

########################################################################################
## 2019-2020 early season increases: August to mid-Nov 2019
## cross-correlations between Rt and mobility
########################################################################################

########################################################################################
## child daycare

all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4)

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, start_week) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hpiv, rsv, adv, flu, hcov, hmpv

## min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

########################################################################################
## schools
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, start_week) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4)
# pathogen         mobility_metric                    n
# <fct>            <fct>                          <int>
# 1 Influenza B      "Elementary and\nhigh schools"    20
# 2 RSV A            "Elementary and\nhigh schools"    20
# 3 hMPV             "Elementary and\nhigh schools"    20
# 4 RSV B            "Elementary and\nhigh schools"    16
# 5 hCoV HKU1 + NL63 "Elementary and\nhigh schools"    15
# 6 Influenza A/H1N1 "Elementary and\nhigh schools"    14
# 7 hPIV 1 + 2       "Elementary and\nhigh schools"    10
# 8 Adenovirus       "Elementary and\nhigh schools"     9
# 9 hCoV 229E + OC43 "Elementary and\nhigh schools"     7
# 10 hPIV 3 + 4       "Elementary and\nhigh schools"     7

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 7) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, rsv, hpiv, hCoV, hmpv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

########################################################################################
## colleges
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, start_week) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 6)

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 6) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  filter(obs_max_ccf >= 0.4) %>%
  arrange(-obs_max_ccf)
# hmpv, hpiv, flu, rsv, hcov

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  filter(obs_max_ccf >= 0.4) %>%
  arrange(-obs_max_ccf)

########################################################################################
## religious orgs
## at least 4 weeks in the season that coefficient is significant
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 3)
# pathogen         mobility_metric                n
# <fct>            <fct>                      <int>
# 1 hMPV             "Religious\norganizations"    16
# 2 Influenza B      "Religious\norganizations"    12
# 3 RSV A            "Religious\norganizations"    11
# 4 hPIV 1 + 2       "Religious\norganizations"    10
# 5 hPIV 3 + 4       "Religious\norganizations"     8
# 6 RSV B            "Religious\norganizations"     7
# 7 hCoV HKU1 + NL63 "Religious\norganizations"     7

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep


## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, rsv, hmpv, hcov, hpiv, adv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

########################################################################################
## restaurants
# all_path_ccf %>%
#   filter(pathogen != "SARS-CoV-2" & mobility_metric == "Restaurants") %>%
#   filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
#   filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
#   distinct() %>%
#   arrange(pathogen, mobility_metric, obs_max_ccf) %>%
#   group_by(pathogen, mobility_metric) %>%
#   tally() %>%
#   arrange(-n) %>%
#   filter(n > 3)
# # pathogen         mobility_metric     n
# # <fct>            <fct>           <int>
# # 1 hPIV 3 + 4       Restaurants        17
# # 2 hCoV HKU1 + NL63 Restaurants         8
# # 3 hCoV 229E + OC43 Restaurants         6
# # 4 RSV B            Restaurants         5
#
# ## max values
# all_path_ccf %>%
#   filter(pathogen != "SARS-CoV-2" & mobility_metric == "Restaurants") %>%
#   filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
#   filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
#   distinct() %>%
#   arrange(pathogen, mobility_metric, obs_max_ccf) %>%
#   group_by(pathogen, mobility_metric) %>%
#   slice_max(obs_max_ccf) %>%
#   filter(obs_max_ccf >= 0.4) %>%
#   arrange(-obs_max_ccf)
# #hpiv, hcov, rsv b
#
# ## min values
# all_path_ccf %>%
#   filter(pathogen != "SARS-CoV-2" & mobility_metric == "Restaurants") %>%
#   filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
#   filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-15")) %>%
#   distinct() %>%
#   arrange(pathogen, mobility_metric, obs_max_ccf) %>%
#   group_by(pathogen, mobility_metric) %>%
#   slice_max(obs_max_ccf) %>%
#   arrange(-obs_max_ccf)

########################################################################################
## b/w neighborhood movement

## at least 4 weeks in the season that coefficient is significant
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.4) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 3)
# pathogen         mobility_metric                      n
# <fct>            <fct>                            <int>
# 1 Adenovirus       "Between-neighborhood\nmovement"    17
# 2 hPIV 1 + 2       "Between-neighborhood\nmovement"    13
# 3 Influenza B      "Between-neighborhood\nmovement"     7
# 4 RSV A            "Between-neighborhood\nmovement"     7
# 5 hMPV             "Between-neighborhood\nmovement"     7
# 6 hCoV HKU1 + NL63 "Between-neighborhood\nmovement"     6
# 7 RSV B            "Between-neighborhood\nmovement"     5

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-08-15")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 5) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep
## strongest correlations later in the season
## max values
all_path_ccf %>%
  filter(pathogen %in% keep & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-02-01") & start_week > as.Date("2019-10-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, start_week) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  # filter(obs_max_ccf>0.5)%>%
  arrange(-obs_max_ccf)
## late nov/dec

all_path_ccf %>%
  filter(pathogen %in% keep & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-02-01") & start_week > as.Date("2019-10-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "Adenovirus" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "hPIV 1 + 2" & mobility_metric == "Between-neighborhood\nmovement") %>% # 3 weeks
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "Influenza B" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) # 5 weeks

########################################################################################
## within-neighborhood movement
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Within-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## few/weak correlations

########################################################################################
## visitor inflow
## out-of-state
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Influx out-of-state\nvisitors") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 3)
# 1 hPIV 1 + 2  "Influx out-of-state\nvisitors"     7
# 2 hMPV        "Influx out-of-state\nvisitors"     5

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Influx out-of-state\nvisitors") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
# hpiv, hmpv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Influx out-of-state\nvisitors") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

## inflow from other WA counties
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Infux visitors\nother WA counties") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(pathogen, obs_max_ccf)
## only one positive, leading relationships

########################################################################################
## percent leaving home
## at least 4 weeks in the season that coefficient is significant
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 5)
# pathogen         mobility_metric               n
# <fct>            <fct>                     <int>
# 1 RSV A            "% Devices\nleaving home"    15
# 2 RSV B            "% Devices\nleaving home"    14
# 3 hCoV HKU1 + NL63 "% Devices\nleaving home"    11
# 4 Influenza A/H1N1 "% Devices\nleaving home"     8
# 5 hMPV             "% Devices\nleaving home"     7
# 6 hPIV 3 + 4       "% Devices\nleaving home"     7

keep <- all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 5) %>%
  pull(pathogen) %>%
  unique() %>%
  as.character()
keep

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-15") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
# hcov, rsv, hpiv, rv, flu

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & pathogen %in% keep & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-15") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

########################################################################################
##### non-env viruses, esp. RV, have fewer relationships with mobility during fall 2019
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  filter(pathogen == "Adenovirus") %>%
  arrange(mobility_metric, -obs_max_ccf)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  filter(pathogen == "Rhinovirus") %>%
  arrange(mobility_metric, -obs_max_ccf)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  filter(n > 4) %>%
  slice_max(n, n = 5) %>%
  arrange(pathogen, -n) %>%
  print(n = 50)

########################################################################################
## Figure S11. Plot early 2020 moving window cross-correlations between Rt and mobility
########################################################################################
## take out masking and OSI because these data aren't available until after March and June 2020
df1 <- all_path_ccf %>%
  filter(start_week < as.Date("2020-06-30") & # non-env viruses have continuous time series
           start_week > as.Date("2019-12-01") &
           !(mobility_metric %in% c("% not wearing masks", "Oxford Stringency\nIndex"))) %>%
  droplevels()

df2 <- df1 %>%
  filter(!(pathogen %in% c("Adenovirus", "Rhinovirus", "Enterovirus", "SARS-CoV-2"))) %>%
  filter(start_week < as.Date("2020-05-01")) %>% # Rt goes to 0 for env viruses
  droplevels()

df3 <- df1 %>%
  filter(pathogen %in% c("Adenovirus", "Rhinovirus", "Enterovirus", "SARS-CoV-2")) %>%
  droplevels()

df4 <- bind_rows(df2, df3)

levels(df4$mobility_metric)
levels(df4$pathogen)

## cross-correlation plot
spring_2020_plot <- ggplot(
  df4 %>% filter(obs_max_ccf_lag < 1),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  annotate("rect",
           xmin = as.Date("2020-03-23"),
           xmax = as.Date("2020-06-05"),
           ymin = -Inf,
           ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 0.7) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed", lwd = 1.2) +
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(begin = 0, end = 0.75) +
  xlab("5-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Weeks)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )
spring_2020_plot

keep <- c(
  "Influenza A/H1N1", "Influenza B",
  "RSV A",
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  "hPIV 1 + 2",
  "hPIV 3 + 4",
  "Adenovirus", "Enterovirus", "Rhinovirus",
  "SARS-CoV-2"
)

# rt plot
spring_2020_rt <- ggplot() +
  annotate("rect",
           xmin = as.Date("2020-03-23"),
           xmax = as.Date("2020-06-05"),
           ymin = -Inf,
           ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 0.7) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism %in% keep) %>%
      filter(epi_date > as.Date("2019-12-01") &
               epi_date < as.Date("2020-06-30")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.2, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(epi_date > as.Date("2019-12-01") &
               epi_date < as.Date("2020-06-30")) %>%
      filter(organism %in% keep) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.9
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec) +
  facet_wrap(~organism, nrow = 1, scales = "free_y") +
  ylab("Rt") +
  xlab("Epidemic Week")
spring_2020_rt

rt_top <- plot_grid(NULL, spring_2020_rt, NULL, rel_widths = c(0.1, 9, 0.9), nrow = 1)
rt_mob_spring_2020 <- plot_grid(rt_top, spring_2020_plot, nrow = 2, rel_heights = c(1.5, 12), labels = "AUTO")
rt_mob_spring_2020
save_plot(rt_mob_spring_2020, file = "figures/fig_s11_block_bootstrap_spring_2020_spearman_neg_lags.png", base_width = 24, base_height = 16)


# keep <- c("Rhinovirus", "Enterovirus", "Adenovirus", "SARS-CoV-2", "RSV B", "hMPV", "hCoV 229E + OC43", "hCoV HKU1 + NL63", "hPIV 3 + 4")
