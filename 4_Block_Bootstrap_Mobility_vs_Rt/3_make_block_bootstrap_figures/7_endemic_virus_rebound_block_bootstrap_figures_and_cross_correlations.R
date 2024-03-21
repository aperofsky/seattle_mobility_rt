########################################################################################
## Moving window Spearman cross-correlations between Rt and mobility, Endemic virus rebound
## Figure S13. Non-enveloped virus moving window cross-correlations (pandemic rebound)
## Significant cross-correlations between non-enveloped virus Rt and mobility during 2020 (post SAH)
## Figure S15. Enveloped virus moving window cross-correlations (pandemic rebound)
## Significant cross-correlations between enveloped virus Rt and mobility during rebound in 2021
## Significant cross-correlations between endemic virus Rt and mobility during Omicron BA.1 wave
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

########################################################################################
######## Figure S13. Non-enveloped virus cross-correlations during pandemic rebound (RV and AdV)
########################################################################################

keep1 <- c("Rhinovirus", "Enterovirus", "Adenovirus")

unique(all_path_ccf$mobility_metric)

all_mob_plot2 <- ggplot(
  all_path_ccf %>%
    filter(obs_max_ccf_lag < 1) %>%
    filter(start_week > as.Date("2020-01-01") & pathogen %in% keep1 &
      !(mobility_metric %in% c("Groceries and\npharmacies", "Transit", "Within-neighborhood\nmovement"))),
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
    xmin = as.Date("2021-11-01"),
    xmax = as.Date("2022-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed", lwd = 1.2) +
  scale_alpha_manual(values = c(0.1, 0.9), name = "p(perm) < 0.01") +
  theme_bw(base_size = 18) +
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
all_mob_plot2

color_vec_lim <- c("#882255", "#993377", "#AA4499")
cont_rt <- ggplot() +
  annotate("rect",
    xmin = as.Date("2020-03-23"),
    xmax = as.Date("2020-06-05"),
    ymin = -Inf,
    ymax = Inf, fill = "orange", alpha = 0.2
  ) +
  annotate("rect",
    xmin = as.Date("2021-11-01"),
    xmax = as.Date("2022-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism %in% keep1) %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date <= as.Date("2022-03-06")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.2, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(organism %in% keep1) %>%
      filter(epi_date > as.Date("2020-01-01") & epi_date <= as.Date("2022-03-06")) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.9
  ) +
  theme_bw(base_size = 18) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank()) +
  scale_color_manual(values = color_vec_lim) +
  scale_fill_manual(values = color_vec_lim) +
  facet_wrap(~organism, nrow = 1) +
  ylab("Rt") +
  xlab("Epidemic Week")
cont_rt

rt_top <- plot_grid(NULL, cont_rt, NULL, rel_widths = c(0.08, 9, 1.1), nrow = 1)
rt_mob_cont <- plot_grid(rt_top, all_mob_plot2, nrow = 2, rel_heights = c(2.5, 12), labels = "AUTO")
rt_mob_cont
save_plot(rt_mob_cont, file = "figures/fig_s13_block_bootstrap_non_env_viruses_5mo_moving_window_spearman_neg_lags.png", base_width = 22, base_height = 16)

########################################################################################
## Non-enveloped virus cross-correlations between Rt and mobility (pandemic rebound during 2020)
########################################################################################

########################################################################################
## Rhinovirus
## statistically significant cross-correlations between Rt and mobility
########################################################################################
non_env_path <- c("RV", "Adenovirus", "EV")
non_env_rt <- read_csv("4_Block_Bootstrap_Mobility_vs_Rt/rt_all_pathogens_15day_mv_avg.csv") %>% filter(organism %in% non_env_path)

ggplot(data = non_env_rt %>%
  filter(date > as.Date("2020-06-01") & date < as.Date("2021-01-01"))) +
  geom_line(aes(x = date, y = median, color = organism)) +
  theme_bw()

## positive correlations
all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric, pathogen) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4)
# mobility_metric                     pathogen       n
# <fct>                               <fct>      <int>
# 1 "Religious\norganizations"          Rhinovirus    30
# 2 "Infux visitors\nother WA counties" Rhinovirus    25
# 3 "Between-neighborhood\nmovement"    Rhinovirus    14
# 4 "% Devices\nleaving home"           Rhinovirus    11
# 5 "Transit"                           Rhinovirus    11
# 6 "Influx out-of-state\nvisitors"     Rhinovirus     9
# 7 "Within-neighborhood\nmovement"     Rhinovirus     8
# 8 "Restaurants"                       Rhinovirus     7
# 9 "Child daycare"                     Rhinovirus     7
# 10 "Colleges"                          Rhinovirus     7

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(mobility_metric %in% c("Infux visitors\nother WA counties", "Religious\norganizations")) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "Infux visitors\nother WA counties" 0.547 0.974
# 2 "Religious\norganizations"          0.580 0.828

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-30")) %>%
  distinct() %>%
  arrange(start_week) # June 2020 to early Jan 2021

keep <- all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-08-15")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 5) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(mobility_metric %in% keep) %>%
  filter(!(mobility_metric %in% c("Infux visitors\nother WA counties", "Religious\norganizations"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-08-15")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# 1 "Influx out-of-state\nvisitors"  0.665 0.961
# 2 "Between-neighborhood\nmovement" 0.615 0.920
# 3 "Transit"                        0.653 0.881
# 4 "Child daycare"                  0.599 0.848
# 5 "% Devices\nleaving home"        0.612 0.843
# 6 "Restaurants"                    0.585 0.831
# 7 "Within-neighborhood\nmovement"  0.535 0.771

# negative correlations
all_path_ccf %>%
  filter(pathogen == "Rhinovirus" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)
#            mobility_metric start_week obs_max_ccf_lag obs_max_ccf sig   pathogen
# 1 Oxford Stringency\nIndex 2020-08-23              -1  -0.5823216 yes Rhinovirus
# 2 Oxford Stringency\nIndex 2020-08-30              -1  -0.6070566 yes Rhinovirus
# 3 Oxford Stringency\nIndex 2020-09-06              -1  -0.6918868 yes Rhinovirus
# 4 Oxford Stringency\nIndex 2020-09-13               0  -0.8204624 yes Rhinovirus
# 5 Oxford Stringency\nIndex 2020-09-20               0  -0.8232969 yes Rhinovirus
# 6 Oxford Stringency\nIndex 2020-09-27               0  -0.7875638 yes Rhinovirus
# 7 Oxford Stringency\nIndex 2020-10-04               0  -0.7216602 yes Rhinovirus
########################################################################################
## Adenovirus
## statistically significant cross-correlations between Rt and mobility
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)

all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4)
# mobility_metric                    n
# <fct>                          <int>
# 1 "Elementary and\nhigh schools"         13
# 2 "Religious\norganizations"             12
# 3 "Colleges"                              8
# 4 "Infux visitors\nother WA counties"     7
# 5 "Within-neighborhood\nmovement"         5
# 6 "Influx out-of-state\nvisitors"         5

all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric %in%
    c(
      "Elementary and\nhigh schools", "Religious\norganizations", "Colleges",
      "Infux visitors\nother WA counties"
    )) %>%
  filter(start_week >= as.Date("2020-06-01") & start_week <= as.Date("2020-10-15")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                   min   max
# <fct>                           <dbl> <dbl>
# 1 "Religious\norganizations"          0.603 0.881
# 2 "Elementary and\nhigh schools"      0.552 0.840
# 3 "Colleges"                          0.580 0.791
# 4 "Infux visitors\nother WA counties" 0.615 0.649
## rel. orgs: early June to late Aug
## schools: early June to late Aug
## colleges: early june and late july

all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf)

# negative correlations with OSI
all_path_ccf %>%
  filter(pathogen == "Adenovirus" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)
## June to late July

########################################################################################
## Enterovirus
## statistically significant cross-correlations between Rt and mobility
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(pathogen == "Enterovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-10-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)

keep <- all_path_ccf %>%
  filter(pathogen == "Enterovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-10-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 1) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

all_path_ccf %>%
  filter(pathogen == "Enterovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-10-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n)
# mobility_metric                           n
# 1 "Within-neighborhood\nmovement"         7
# 2 "Religious\norganizations"              7
# 3 "Influx out-of-state\nvisitors"         6
# 4 "Transit"                               6
# 5 "Infux visitors\nother WA counties"     5
# 6 "Elementary and\nhigh schools"          5
# 7 "% Devices\nleaving home"               3
# 8 "Restaurants"                           3
# 9 "Colleges"                              3

all_path_ccf %>%
  filter(pathogen == "Enterovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric %in% keep) %>%
  filter(start_week >= as.Date("2020-06-01") & start_week <= as.Date("2020-10-15")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                   min   max
# <fct>                           <dbl> <dbl>
# 1 "Religious\norganizations"          0.607 0.806
# 2 "Influx out-of-state\nvisitors"     0.576 0.806
# 3 "Infux visitors\nother WA counties" 0.587 0.737
# 4 "Elementary and\nhigh schools"      0.613 0.735
# 5 "Colleges"                          0.689 0.727
# 6 "Transit"                           0.622 0.708
# 7 "% Devices\nleaving home"           0.569 0.630
# 8 "Within-neighborhood\nmovement"     0.556 0.618
# 9 "Restaurants"                       0.533 0.605

all_path_ccf %>%
  filter(pathogen == "Enterovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf)

# negative correlations with OSI
all_path_ccf %>%
  filter(pathogen == "Enterovirus" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  # filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)
## no significant correlations, except 2 weeks in Oct 2021

########################################################################################
## Figure S15. Enveloped virus cross-correlations (pandemic rebound)
########################################################################################

keep2 <- c(
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  # "hPIV 1 + 2",
  "hPIV 3 + 4"
)

env_virus_rebound_plot <- ggplot(
  all_path_ccf %>% filter(obs_max_ccf_lag < 1) %>%
    filter(start_week > as.Date("2020-11-01") & start_week < as.Date("2022-06-01") & pathogen %in% keep2 &
      !(mobility_metric %in% c("Groceries and\npharmacies", "Transit", "Within-neighborhood\nmovement"))),
  aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  annotate("rect",
    xmin = as.Date("2021-11-01"),
    xmax = as.Date("2022-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  geom_vline(xintercept = as.Date("2021-04-19"), lty = "dashed", color = "darkgreen", lwd = 1) +
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
    plot.title = element_text(hjust = 0.5),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )

env_virus_rebound_plot

keep2 <- c(
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  # "hPIV 1 + 2",
  "hPIV 3 + 4"
)


color_vec3 <- c("#117733", "#999933", "#DDCC77", "#661100", "#AA4466")


env_rt <- ggplot() +
  geom_vline(xintercept = as.Date("2021-04-19"), lty = "dashed", color = "darkgreen", lwd = 1) +
  annotate("rect",
    xmin = as.Date("2021-11-01"),
    xmax = as.Date("2022-02-01"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = rt_weekly %>%
      filter(organism %in% keep2) %>%
      filter(epi_date > as.Date("2021-01-16") & epi_date < as.Date("2022-03-07")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.2, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(organism %in% keep2) %>%
      filter(epi_date > as.Date("2021-01-16") & epi_date < as.Date("2022-03-07")) %>%
      dplyr::select(epi_date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, ymin = lower, ymax = upper, fill = organism), alpha = 0.9
  ) +
  theme_bw(base_size = 18) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank()) +
  scale_color_manual(values = color_vec3) +
  scale_fill_manual(values = color_vec3) +
  facet_wrap(~organism, nrow = 1) +
  ylab("Rt") +
  xlab("Epidemic Week")
env_rt

rt_top <- plot_grid(NULL, env_rt, NULL, rel_widths = c(0.1, 9, 1.1), nrow = 1)
rt_mob_env <- plot_grid(rt_top, env_virus_rebound_plot, nrow = 2, rel_heights = c(2, 12), labels = "AUTO")
rt_mob_env
save_plot(rt_mob_env, file = "figures/fig_s15_block_bootstrap_env_viruses_rebound_5mo_moving_window_spearman_neg_lags.png", base_width = 22, base_height = 14)


########################################################################################
## Enveloped virus cross-correlations between Rt and mobility (pandemic rebound during 2021)
########################################################################################

########################################################################################
## RSV B, summer 2021
## statistically significant cross-correlations between Rt and mobility
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
    mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, start_week)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n)
# mobility_metric                      n
# # <fct>                            <int>
# 1 "Religious\norganizations"          15
# 2 "Between-neighborhood\nmovement"    13
# 3 "% not wearing masks"               11
# 4 "% Devices\nleaving home"           10
# 5 "Influx out-of-state\nvisitors"      9
# 6 "Restaurants"                        9
# 7 "Child daycare"                      9
# 8 "Elementary and\nhigh schools"       9

keep <- all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 &
    mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()
keep

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  filter(mobility_metric %in% keep) %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                    min   max
# <fct>                            <dbl> <dbl>
# 1 "Child daycare"                  0.641 0.883
# 2 "Between-neighborhood\nmovement" 0.582 0.861
# 3 "% Devices\nleaving home"        0.590 0.860
# 4 "Elementary and\nhigh schools"   0.579 0.854
# 5 "% not wearing masks"            0.570 0.837
# 6 "Religious\norganizations"       0.626 0.830
# 7 "Restaurants"                    0.596 0.827
# 8 "Influx out-of-state\nvisitors"  0.574 0.818

all_path_ccf %>%
  filter(mobility_metric == "Between-neighborhood\nmovement" & pathogen == "RSV B") %>%
  filter(start_week > as.Date("2021-03-01") & start_week < as.Date("2021-09-01")) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)

## max correlation values: early to mid April
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.5) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)

# min correlation values
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.3) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(start_week)
# early May to late June


## Fall 2019
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric %in% keep) %>%
  filter(start_week < as.Date("2020-03-01") & start_week > as.Date("2019-07-15")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, start_week) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n)
# 1 "Elementary and\nhigh schools"      22
# 2 "Religious\norganizations"          17
# 3 "% Devices\nleaving home"           16
# 4 "Child daycare"                     15
# 5 "Between-neighborhood\nmovement"    13
# 6 "Restaurants"                       13
# 7 "Influx out-of-state\nvisitors"      5
# 3 - 5.5 mo
########################################################################################
## hPIV 3+4, hMPV, hCoV; Spring 2021
## statistically significant cross-correlations between Rt and mobility
########################################################################################
# hpiv 3+4
## check number of statistically significant windows
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(pathogen == "hPIV 3 + 4") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, start_week, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(pathogen, -n) %>%
  filter(n >= 4) %>%
  print(n = 10)
# pathogen   mobility_metric                      n
# <fct>      <fct>                            <int>
# 1 hPIV 3 + 4 "Restaurants"                          15
# 2 hPIV 3 + 4 "Elementary and\nhigh schools"         15
# 3 hPIV 3 + 4 "% not wearing masks"                  13
# 4 hPIV 3 + 4 "Between-neighborhood\nmovement"       11
# 5 hPIV 3 + 4 "% Devices\nleaving home"              10
# 6 hPIV 3 + 4 "Infux visitors\nother WA counties"     7
# 7 hPIV 3 + 4 "Influx out-of-state\nvisitors"         7
# 8 hPIV 3 + 4 "Child daycare"                         7
# 9 hPIV 3 + 4 "Religious\norganizations"              6

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.5) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  filter(pathogen == "hPIV 3 + 4") %>%
  arrange(mobility_metric, start_week, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  filter(pathogen == "hPIV 3 + 4") %>%
  arrange(mobility_metric, start_week, obs_max_ccf)

# hMPV
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  filter(pathogen == "hMPV") %>%
  arrange(mobility_metric, start_week, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4) %>%
  print(n = 10)
# mobility_metric                         n
# <fct>                               <int>
# 1 "% Devices\nleaving home"              11
# 2 "Between-neighborhood\nmovement"       10
# 3 "Restaurants"                          10
# 4 "Elementary and\nhigh schools"         10
# 5 "% not wearing masks"                  10
# 6 "Within-neighborhood\nmovement"         9
# 7 "Influx out-of-state\nvisitors"         8
# 8 "Child daycare"                         8
# 9 "Infux visitors\nother WA counties"     6
# 10 "Religious\norganizations"              6

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  filter(pathogen == "hMPV") %>%
  arrange(mobility_metric, start_week, obs_max_ccf)

# "hCoV 229E + OC43"
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  distinct() %>%
  filter(pathogen == "hCoV 229E + OC43") %>%
  arrange(mobility_metric, start_week)

# "hCoV HKU1 + NL63"
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  distinct() %>%
  filter(pathogen == "hCoV HKU1 + NL63") %>%
  arrange(mobility_metric, -obs_max_ccf)

range(all_path_ccf$start_week)

########################################################################################
## Endemic pathogen declines during Omicron BA.1 wave late 2021
## statistically significant cross-correlations between Rt and mobility
########################################################################################
## check number of statistically significant windows
## endemic env viruses
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(!(pathogen %in% c("SARS-CoV-2"))) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n > 4) %>%
  arrange(pathogen, -n) %>%
  print(n = 60)

# rhinovirus
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "Rhinovirus") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(-n) %>%
  print(n = 50)
# mobility_metric                     obs_max_ccf_lag     n
# <fct>                                         <dbl> <int>
# 1 "Transit"                                         0    13
# 2 "Colleges"                                        0     9
# 3 "Restaurants"                                     0     5
# 4 "Infux visitors\nother WA counties"               0     4
# 5 "% Devices\nleaving home"                         0     4
# 6 "Religious\norganizations"                        0     4
# 7 "Child daycare"                                   0     4
# 8 "Colleges"                                       -1     4
# 9 "Influx out-of-state\nvisitors"                   0     2
# 10 "% not wearing masks"                             0     2
# 11 "Elementary and\nhigh schools"                    0     1

# adenovirus
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "Adenovirus") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(mobility_metric, -n) %>%
  print(n = 50)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "Enterovirus") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(mobility_metric, -n) %>%
  print(n = 50)

## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  arrange(pathogen, mobility_metric)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  arrange(pathogen, mobility_metric) %>%
  mutate(sync = if_else(obs_max_ccf_lag < 0, "yes", "no")) %>%
  group_by(pathogen, mobility_metric, sync) %>%
  tally() %>%
  arrange(pathogen, mobility_metric, -n) %>%
  filter(n > 2) %>%
  print(n = 80)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(pathogen, mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(pathogen, -obs_max_ccf_lag) %>%
  filter(obs_max_ccf_lag != 0)

# % devices leaving home
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.5) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  filter(mobility_metric == "% Devices\nleaving home") %>%
  arrange(pathogen, mobility_metric)
# RSV B, hMPV, hCoV, RV, AdV
