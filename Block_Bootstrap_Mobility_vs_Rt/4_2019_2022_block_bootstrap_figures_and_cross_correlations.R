########################################################################################
## Moving window cross-correlations between Rt and mobility, Fall 2019 to Winter 2022
## 1. Figure S8. Season 2019-2020 moving window cross-correlations (pre-pandemic)
## 2. Significant cross-correlations in Fall 2019
## 3. Figure S10: late 2019 - early 2020 moving window cross-correlations (early pandemic)
## 4. Figure 5: SARS-CoV-2 moving window cross-correlations
## 5. Figure S12. Non-enveloped virus moving window cross-correlations (pandemic rebound)
## 6. Figure S14. Enveloped virus moving window cross-correlations (pandemic rebound)
## 7. Significant cross-correlations between SARS-CoV-2 Rt and mobility during each wave
## 8. Significant cross-correlations between endemic virus Rt and mobility during pandemic rebound
## 8a. non-enveloped virus
## 8b. enveloped viruses
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
## Import block bootstrap results
########################################################################################

dir <- "Block_Bootstrap_Mobility_vs_Rt/biowulf_block_bootstrap_5mo_rolling_window_exp_decay/biowulf_block_bootstrap_5mo_rolling_window_exp_decay_output/"

## adenovirus
load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
head(actual_data_and_perm)
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "adeno"

load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))
actual_data_and_perm2$pathogen <- "adeno"

load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))
actual_data_and_perm3$pathogen <- "adeno"

all_adeno <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
range(all_adeno$start_week)
unique(all_adeno$mobility_metric)

########################################################################################
### covid
load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "covid"

load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))
actual_data_and_perm2$pathogen <- "covid"

load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))
actual_data_and_perm3$pathogen <- "covid"

all_covid <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
unique(all_covid$mobility_metric)

########################################################################################
## seasonal CoV
load(paste0(dir, "hCoV_229E_OC43_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "scov_229E_OC43"

load(paste0(dir, "scov_229E_OC43_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021.RData"))
actual_data_and_perm2 <- actual_data_and_perm
actual_data_and_perm2$pathogen <- "scov_229E_OC43"

all_scov_229E_OC43 <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_scov_229E_OC43$mobility_metric)

load(paste0(dir, "hCoV_HKU1_NL63_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "scov_HKU1_NL63"

load(paste0(dir, "scov_HKU1_NL63_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021.RData"))
actual_data_and_perm2 <- actual_data_and_perm
actual_data_and_perm2$pathogen <- "scov_HKU1_NL63"

all_scov_HKU1_NL63 <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_scov_HKU1_NL63$mobility_metric)

########################################################################################
## hmpv
load(paste0(dir, "hmpv_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "hmpv"

load(paste0(dir, "hmpv_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021.RData"))
actual_data_and_perm2 <- actual_data_and_perm
actual_data_and_perm2$pathogen <- "hmpv"

all_hmpv <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_hmpv$mobility_metric)

########################################################################################
## HPIV
load(paste0(dir, "hpiv_3_4_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "hpiv_3_4"

load(paste0(dir, "hpiv_3_4_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021.RData"))
actual_data_and_perm2 <- actual_data_and_perm
actual_data_and_perm2$pathogen <- "hpiv_3_4"

all_hpiv <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
range(all_hpiv$start_week)
unique(all_hpiv$mobility_metric)

load(paste0(dir, "hpiv_1_2_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
all_hpiv_1_2 <- actual_data_and_perm
all_hpiv_1_2$pathogen <- "hpiv_1_2"
unique(all_hpiv_1_2$mobility_metric)

########################################################################################
## rhino
load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "rhino"

load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))
actual_data_and_perm2$pathogen <- "rhino"

load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))
actual_data_and_perm3$pathogen <- "rhino"

all_rhino <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
sort(unique(all_rhino$mobility_metric))

########################################################################################
## rsv a
load(paste0(dir, "rsv_a_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
all_rsv_a <- actual_data_and_perm
all_rsv_a$pathogen <- "rsv_a"

########################################################################################
## rsv b
load(paste0(dir, "rsv_b_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
actual_data_and_perm1 <- actual_data_and_perm
actual_data_and_perm1$pathogen <- "rsv_b"

load(paste0(dir, "rsv_b_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021.RData"))
actual_data_and_perm2 <- actual_data_and_perm
actual_data_and_perm2$pathogen <- "rsv_b"

all_rsv_b <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_rsv_b$mobility_metric)

########################################################################################
### h1n1
load(paste0(dir, "h1n1_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
all_h1n1 <- actual_data_and_perm
all_h1n1$pathogen <- "h1n1"

########################################################################################
# ivb
load(paste0(dir, "ivb_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData"))
all_ivb <- actual_data_and_perm
all_ivb$pathogen <- "ivb"

########################################################################################
## combine data for all pathogens
########################################################################################
all_path_ccf <- bind_rows(
  all_adeno, all_covid, all_hmpv, all_h1n1, all_ivb,
  all_hpiv, all_hpiv_1_2, all_rhino, all_rsv_a,
  all_rsv_b, all_scov_229E_OC43, all_scov_HKU1_NL63
)

all_path_ccf$pathogen <- as.factor(all_path_ccf$pathogen)
levels(all_path_ccf$pathogen)
levels(all_path_ccf$mobility_metric)
unique(all_path_ccf$pathogen)
all_path_ccf$pathogen <- factor(all_path_ccf$pathogen, levels = c(
  "h1n1", "ivb", "rsv_a",
  "rsv_b", "hmpv",
  "scov_229E_OC43",
  "scov_HKU1_NL63",
  "hpiv_1_2",
  "hpiv_3_4",
  "rhino", "adeno", "covid"
))

levels(all_path_ccf$pathogen)

levels(all_path_ccf$pathogen) <- c(
  "Influenza A/H1N1", "Influenza B", "RSV A", "RSV B", "hMPV",
  "hCoV 229E + OC43", "hCoV HKU1 + NL63",
  "hPIV 1 + 2", "hPIV 3 + 4", "Rhinovirus",
  "Adenovirus", "SARS-CoV-2"
)

levels(all_path_ccf$mobility_metric)

all_path_ccf <- all_path_ccf %>%
  filter(!(mobility_metric %in% c("performing_arts_or_sports_events"))) %>%
  droplevels()
unique(all_path_ccf$mobility_metric)

all_path_ccf$mobility_metric <- factor(all_path_ccf$mobility_metric,
  levels = c(
    "within_neighborhood_movement",
    "within_city_movement",
    "within_state_movement",
    "out_of_state_movement",
    "fb_leaving_home_custom",
    "full_service_restaurants",
    "groceries_and_pharmacies",
    "transit",
    "religious_orgs", "child_day_care",
    "elementary_and_secondary_schools",
    "colleges", "not_masking",
    "oxford_stringency_index"
  )
)
levels(all_path_ccf$mobility_metric)

levels(all_path_ccf$mobility_metric) <- c(
  "Within-neighborhood\nmovement",
  "Between-neighborhood\nmovement",
  "Infux visitors\nother WA counties",
  "Influx out-of-state\nvisitors",
  "% Devices\nleaving home",
  "Restaurants",
  "Groceries and\npharmacies",
  "Transit",
  "Religious\norganizations",
  "Child daycare",
  "Elementary and\nhigh schools",
  "Colleges",
  "% not wearing masks",
  "Oxford Stringency\nIndex"
)

all_endemic_results_avg <- read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")
all_endemic_results_avg$organism <- as.factor(all_endemic_results_avg$organism)
levels(all_endemic_results_avg$organism)

levels(all_endemic_results_avg$organism) <- c(
  "Adenovirus",
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
# color_vec = c("#332288","#6699CC","#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77","#661100", "#CC6677", "#AA4466","#882255","#AA4499")

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
  "#AA4499",
  "#ABB065"
)

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
  "Adenovirus", "Rhinovirus",
  "SARS-CoV-2"
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
        epi_date < as.Date("2020-01-15")) %>%
      dplyr::select(epi_date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = epi_date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = rt_weekly %>%
      filter(epi_date > as.Date("2019-08-20") &
        epi_date < as.Date("2020-01-15")) %>%
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
  all_path_ccf %>% filter(start_week > as.Date("2019-08-20") &
    start_week < as.Date("2020-01-15") &
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
  all_path_ccf %>%
    filter(start_week > as.Date("2019-08-20") & start_week < as.Date("2020-01-01") &
      pathogen != "SARS-CoV-2" &
      mobility_metric %in% c(
        "% Devices\nleaving home", "Between-neighborhood\nmovement",
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
save_plot(rt_mob_fall_2019_red, file = "figures/fig_s8_block_bootstrap_fall_2019_select_mob_metrics.png", base_width = 22, base_height = 10)

########################################################################################
## 2019-2020 early season increases: August to mid-Nov 2019
## cross-correlations between Rt and mobility
########################################################################################

########################################################################################
## child daycare
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hpiv, flu, rsv, adv, hmpv, hcov

## min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Child daycare") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hpiv, flu, rsv, adv, hmpv, hcov

########################################################################################
## schools
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, hmpv, rsv, hCoV, hPIV

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Elementary and\nhigh schools") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, hmpv, rsv, hCoV, hPIV

########################################################################################
## colleges
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, hmpv, rsv, hcov, hpiv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Colleges") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, hmpv, rsv, hcov, hpiv

########################################################################################
## religious orgs
## at least 4 weeks in the season that coefficient is significant
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally()%>%
  arrange(-n)%>%
  filter(n>3)

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## flu, rsv, hmpv, hcov, hpiv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hmpv, rsv, hpiv, flu, hcov

########################################################################################
## restaurants
## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Restaurants") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  filter(obs_max_ccf>0.5)%>%
  arrange(-obs_max_ccf)
## hpiv 3+4, hCoV

## min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Restaurants") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hpiv 3+4, hCoV

########################################################################################
## b/w neighborhood movement

## at least 4 weeks in the season that coefficient is significant
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally()%>%
  arrange(-n)%>%
  filter(n>3)
# 1 Adenovirus       "Between-neighborhood\nmovement"    10
# 2 Influenza B      "Between-neighborhood\nmovement"     8


keep = all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  tally()%>%
  arrange(-n)%>%
  filter(n>4)%>%
  pull(pathogen)%>%
  unique()%>%
  as.character()

## max values
all_path_ccf %>%
  filter(pathogen %in% keep & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  # filter(obs_max_ccf>0.5)%>%
  arrange(-obs_max_ccf)

all_path_ccf %>%
  filter(pathogen %in% keep & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-01") & start_week > as.Date("2019-09-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "Adenovirus" & mobility_metric == "Between-neighborhood\nmovement") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "hPIV 1 + 2" & mobility_metric == "Between-neighborhood\nmovement") %>% #3 weeks
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0)

all_path_ccf %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-09-01")) %>%
  filter(pathogen == "Influenza B" & mobility_metric == "Between-neighborhood\nmovement")%>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) #5 weeks

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
  tally()%>%
  arrange(-n)
# 1 hPIV 1 + 2  "Influx out-of-state\nvisitors"     7
# 2 Rhinovirus  "Influx out-of-state\nvisitors"     2
# 3 Influenza B "Influx out-of-state\nvisitors"     1
# 4 Adenovirus  "Influx out-of-state\nvisitors"     1

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Influx out-of-state\nvisitors") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  filter(obs_max_ccf>0.5)%>%
  arrange(-obs_max_ccf)
## rhinovirus and hpiv 1+2

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "Influx out-of-state\nvisitors") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-11-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## rhinovirus and hpiv 1+2

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
## no positive, leading relationshps

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
  tally()%>%
  arrange(-n)
# pathogen         mobility_metric               n
# <fct>            <fct>                     <int>
# 1 RSV B            "% Devices\nleaving home"    11
# 2 RSV A            "% Devices\nleaving home"     9
# 3 hCoV HKU1 + NL63 "% Devices\nleaving home"     6
# 4 hPIV 3 + 4       "% Devices\nleaving home"     5
# 5 Rhinovirus       "% Devices\nleaving home"     5
# 6 hCoV 229E + OC43 "% Devices\nleaving home"     4

## max values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-01-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_max(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hcov, rsv, hpiv, hcov, rv

# min values
all_path_ccf %>%
  filter(pathogen != "SARS-CoV-2" & mobility_metric == "% Devices\nleaving home") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2019-12-10") & start_week > as.Date("2019-08-01")) %>%
  distinct() %>%
  arrange(pathogen, mobility_metric, obs_max_ccf) %>%
  group_by(pathogen, mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(-obs_max_ccf)
## hcov, rsv, hpiv, hcov, rv

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
  filter(n > 3) %>%
  arrange(pathogen, -n) %>%
  print(n = 50)

########################################################################################
## Figure S10. Plot early 2020 moving window cross-correlations between Rt and mobility
########################################################################################
## take out masking and OSI because these data aren't available until after March and June 2020
df1 <- all_path_ccf %>%
  filter(start_week < as.Date("2020-06-30") & # non-env viruses have continuous time series
    start_week > as.Date("2019-12-01") &
    !(mobility_metric %in% c("% not wearing masks", "Oxford Stringency\nIndex"))) %>%
  droplevels()

df2 <- df1 %>%
  filter(!(pathogen %in% c("Adenovirus", "Rhinovirus", "SARS-CoV-2"))) %>%
  filter(start_week < as.Date("2020-05-01")) %>% # Rt goes to 0 for env viruses
  droplevels()

df3 <- df1 %>%
  filter(pathogen %in% c("Adenovirus", "Rhinovirus", "SARS-CoV-2")) %>%
  droplevels()

df4 <- bind_rows(df2, df3)

levels(df4$mobility_metric)
levels(df4$pathogen)

## cross-correlation plot
spring_2020_plot <- ggplot(
  df4,
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
  "Adenovirus", "Rhinovirus",
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
    labels = scales::number_format(accuracy = 0.1))+
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec) +
  facet_wrap(~organism, nrow = 1,scales="free_y") +
  ylab("Rt") +
  xlab("Epidemic Week")
spring_2020_rt

rt_top <- plot_grid(NULL, spring_2020_rt, NULL, rel_widths = c(0.1, 9, 0.9), nrow = 1)
rt_mob_spring_2020 <- plot_grid(rt_top, spring_2020_plot, nrow = 2, rel_heights = c(1.5, 12), labels = "AUTO")
rt_mob_spring_2020
save_plot(rt_mob_spring_2020, file = "figures/fig_s10_block_bootstrap_spring_2020.png", base_width = 24, base_height = 16)


keep <- c("Rhinovirus", "Adenovirus", "SARS-CoV-2", "RSV B", "hMPV", "hCoV 229E + OC43", "hCoV HKU1 + NL63", "hPIV 3 + 4")

########################################################################################
######## Figure 5: Cross-correlations between Rt and mobility across SARS-CoV-2 waves
########################################################################################

all_mob_plot2 <- ggplot(
  all_path_ccf %>% filter(start_week >= as.Date("2020-02-16") & start_week < as.Date("2022-02-15") &
    pathogen == "SARS-CoV-2"),
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
      filter(epi_date > as.Date("2020-01-01") & epi_date < as.Date("2022-02-15")) %>%
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
save_plot(rt_mob_cont, file = "figures/fig_5_block_bootstrap_covid_5mo_moving_window.png", base_width = 14, base_height = 16)

# covid_mob_red <- ggplot(
#   all_path_ccf %>% filter(start_week >= as.Date("2020-02-16") & start_week < as.Date("2022-02-15") &
#     pathogen == "SARS-CoV-2" &
#     mobility_metric %in% c(
#       "Between-neighborhood\nmovement",
#       "Infux visitors\nother WA counties",
#       "Influx out-of-state\nvisitors",
#       "% Devices\nleaving home",
#       "Restaurants",
#       "% not wearing masks",
#       "Oxford Stringency\nIndex"
#     )),
#   aes(x = start_week, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
# ) +
#   scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = c(0.02, 0.02)) +
#   annotate("rect",
#     xmin = as.Date("2020-03-23"),
#     xmax = as.Date("2020-06-05"),
#     ymin = -Inf,
#     ymax = Inf, fill = "orange", alpha = 0.2
#   ) +
#   annotate("rect",
#     xmin = as.Date("2020-09-25"),
#     xmax = as.Date("2021-02-01"),
#     ymin = -Inf,
#     ymax = Inf, fill = "blue", alpha = 0.1
#   ) +
#   annotate("rect",
#     xmin = as.Date("2021-03-01"),
#     xmax = as.Date("2021-05-31"),
#     ymin = -Inf,
#     ymax = Inf, fill = "blue", alpha = 0.1
#   ) +
#   annotate("rect",
#     xmin = as.Date("2021-07-01"),
#     xmax = as.Date("2021-11-01"),
#     ymin = -Inf,
#     ymax = Inf, fill = "blue", alpha = 0.1
#   ) +
#   annotate("rect",
#     xmin = as.Date("2021-12-01"),
#     xmax = as.Date("2022-03-01"),
#     ymin = -Inf,
#     ymax = Inf, fill = "blue", alpha = 0.1
#   ) +
#   geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen", lwd = 1) +
#   geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
#   geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
#   facet_grid(mobility_metric ~ pathogen) +
#   geom_hline(aes(yintercept = 0), lty = "dashed", lwd = 1.2) +
#   scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
#   theme_bw(base_size = 18) +
#   ylab("Cross-Correlation Coefficient") +
#   scale_fill_viridis_c() +
#   xlab("5-Month Window Mid-Point") +
#   labs(fill = "Optimal Lag\n(Weeks)") +
#   theme(
#     strip.background.x = element_blank(),
#     strip.text.x = element_blank(),
#     strip.background.y = element_rect(fill = "white"),
#     legend.position = "bottom",
#     strip.text.y.right = element_text(angle = 0)
#   )
# covid_mob_red
#
# rt_top <- plot_grid(NULL, cont_rt, NULL, rel_widths = c(0.25, 10, 2), nrow = 1)
# rt_mob_cont <- plot_grid(rt_top, covid_mob_red, nrow = 2, rel_heights = c(2, 7), labels = NULL)
# rt_mob_cont
# save_plot(rt_mob_cont, file = "figures/fig_5_block_bootstrap_covid_5mo_moving_window_select_mobility_indicators.png", base_width = 14, base_height = 12)
########################################################################################
######## Fig S12. Non-enveloped virus cross-correlations during pandemic rebound (RV and AdV)
########################################################################################

keep1 <- c("Rhinovirus", "Adenovirus")

all_mob_plot2 <- ggplot(
  all_path_ccf %>% filter(start_week > as.Date("2020-01-01") &
    pathogen %in% keep1 & mobility_metric != "Hospitals"),
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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

color_vec_lim <- c("#882255", "#AA4499")
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
save_plot(rt_mob_cont, file = "figures/fig_s12_block_bootstrap_non_env_viruses_5mo_moving_window.png", base_width = 22, base_height = 16)

########################################################################################
## Figure S14. Enveloped virus cross-correlations (pandemic rebound)
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
  all_path_ccf %>%
    filter(start_week > as.Date("2020-11-01") &
      start_week < as.Date("2022-06-01") &
      pathogen %in% keep2 & mobility_metric != "Hospitals"),
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 18) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
save_plot(rt_mob_env, file = "figures/fig_s14_block_bootstrap_env_viruses_rebound_5mo_moving_window.png", base_width = 22, base_height = 16)

########################################################################################
## SARS-CoV-2 waves 2020 - 2022
## statistically significant cross-correlations between Rt and mobility
########################################################################################

########################################################################################
## SAH period
########################################################################################

# SC2 date range of significant correlations during SAH: late April to late May 2020
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-07-15") & start_week > as.Date("2020-02-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  pull(start_week) %>%
  range() # "2020-04-26" "2020-05-24"

## positive correlations
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)
## most mobility indicators, late april - may

all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-01")) %>%
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
# 1 "Colleges"                          0.814 0.969
# 2 "Between-neighborhood\nmovement"    0.731 0.961
# 3 "Child daycare"                     0.696 0.960
# 4 "Elementary and\nhigh schools"      0.778 0.959
# 5 "Religious\norganizations"          0.662 0.944
# 6 "Transit"                           0.763 0.927
# 7 "Restaurants"                       0.800 0.918
# 8 "% Devices\nleaving home"           0.737 0.841
# 9 "Infux visitors\nother WA counties" 0.661 0.794
# 10 "Influx out-of-state\nvisitors"     0.749 0.783

## range of positive correlation values
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week >= as.Date("2020-04-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  pull(obs_max_ccf) %>%
  range() # 0.6608960 0.9689244

## negative correlations with OSI
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2020-06-30") & start_week >= as.Date("2020-04-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)

# range of negative correlation values
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2020-06-30") & start_week >= as.Date("2020-04-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf) %>%
  pull(obs_max_ccf) %>%
  range() #-0.9565565 -0.8127755

## filter to max correlation values for each mob metric
## late April to mid-May
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2020-06-01") & start_week > as.Date("2020-04-01")) %>%
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
## positive correlations
all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)

all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric,pathogen)%>%
  tally()%>%
  filter(n>3)
# mobility_metric                     pathogen       n
# <fct>                               <fct>      <int>
# 1 "Infux visitors\nother WA counties" SARS-CoV-2     5
# 2 "Influx out-of-state\nvisitors"     SARS-CoV-2     6
# 3 "% Devices\nleaving home"           SARS-CoV-2     6
# 4 "Restaurants"                       SARS-CoV-2     8

sc2_keep = all_path_ccf %>%
  filter(!(mobility_metric %in% c("% not wearing masks"))) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric,pathogen)%>%
  tally()%>%
  filter(n>3) %>%
  pull(mobility_metric)%>%
  unique()

all_path_ccf %>%
  filter(mobility_metric %in% sc2_keep) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  filter(obs_max_ccf > 0.5) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "% Devices\nleaving home"           0.701 0.827
# 2 "Restaurants"                       0.648 0.826
# 3 "Infux visitors\nother WA counties" 0.702 0.737
# 5 "Influx out-of-state\nvisitors"     0.660 0.690

## negative correlations (OSI)
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  arrange(obs_max_ccf)

all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-10-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  arrange(obs_max_ccf) %>%
  pull(obs_max_ccf) %>%
  range() #-0.8488734 -0.7556842

########################################################################################
## alpha wave, spring 2021
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week < as.Date("2021-05-01") & start_week > as.Date("2021-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)
## only groceries and pharmacies for 2 weeks

# negative correlations
all_path_ccf %>%
  filter(mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week < as.Date("2021-05-01") & start_week > as.Date("2021-02-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, -obs_max_ccf)
## no correlations with OSI

########################################################################################
## delta wave, summer 2021
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric,pathogen)%>%
  tally()%>%
  arrange(-n)
# mobility_metric                     pathogen       n
# <fct>                               <fct>      <int>
# 1 "Influx out-of-state\nvisitors"     SARS-CoV-2     8
# 2 "% Devices\nleaving home"           SARS-CoV-2     6
# 3 "% not wearing masks"               SARS-CoV-2     6
# 4 "Infux visitors\nother WA counties" SARS-CoV-2     1

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  filter(obs_max_ccf > 0.5) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "% not wearing masks"               0.888 0.935
# 2 "Influx out-of-state\nvisitors"     0.836 0.908
# 3 "% Devices\nleaving home"           0.752 0.835

## negative correlations
all_path_ccf %>%
  filter(sig == "yes" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2021-08-15") & start_week > as.Date("2021-06-15")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)
## significant correlation during only 1 week

########################################################################################
## First omicron wave, late 2021 - early 2022
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2022-03-01") & start_week > as.Date("2021-10-30")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)
## only groceries and pharmacies and w/i neighborhood movement have pos correlation with Rt

## negative correlations with OSI
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0 & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(start_week < as.Date("2022-03-01") & start_week > as.Date("2021-10-30")) %>%
  distinct() %>%
  filter(pathogen == "SARS-CoV-2") %>%
  arrange(mobility_metric, obs_max_ccf)
## none

########################################################################################
## Non-enveloped virus cross-correlations between Rt and mobility (pandemic rebound during 2020)
########################################################################################

########################################################################################
## Rhinovirus
## statistically significant cross-correlations between Rt and mobility
########################################################################################

## positive correlations
all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)%>%
  group_by(mobility_metric,pathogen)%>%
  tally()%>%
  arrange(-n)%>%
  filter(n>4)
# mobility_metric                     pathogen       n
# <fct>                               <fct>      <int>
# 1 "Religious\norganizations"          Rhinovirus    18
# 2 "Infux visitors\nother WA counties" Rhinovirus    10
# 3 "Influx out-of-state\nvisitors"     Rhinovirus    10
# 4 "Between-neighborhood\nmovement"    Rhinovirus     9
# 5 "% Devices\nleaving home"           Rhinovirus     9
# 6 "Restaurants"                       Rhinovirus     9

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week)

keep = all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>% 
  group_by(mobility_metric)%>%
  tally()%>%
  filter(n>4)%>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(mobility_metric %in% keep)%>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>%
filter(obs_max_ccf > 0.5) %>%
  group_by(mobility_metric) %>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                       min   max
# <fct>                               <dbl> <dbl>
# 1 "Influx out-of-state\nvisitors"     0.785 0.955
# 2 "Between-neighborhood\nmovement"    0.786 0.925
# 3 "Infux visitors\nother WA counties" 0.761 0.921
# 4 "Restaurants"                       0.758 0.918
# 5 "Religious\norganizations"          0.635 0.910
# 6 "% Devices\nleaving home"           0.770 0.884

all_path_ccf %>%
  filter(pathogen == "Rhinovirus") %>%
  filter(mobility_metric == "Religious\norganizations") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2020-12-15")) %>%
  distinct() %>%
  arrange(-obs_max_ccf)

# negative correlations
all_path_ccf %>%
  filter(pathogen == "Rhinovirus" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf)
#           mobility_metric start_week obs_max_ccf_lag obs_max_ccf sig   pathogen
# 1 Oxford Stringency\nIndex 2020-09-13               0  -0.7903236 yes Rhinovirus
# 2 Oxford Stringency\nIndex 2020-09-20               0  -0.7567030 yes Rhinovirus
# 3 Oxford Stringency\nIndex 2020-09-06               0  -0.7346442 yes Rhinovirus
# 4 Oxford Stringency\nIndex 2020-09-27               0  -0.6919610 yes Rhinovirus
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
  group_by(mobility_metric)%>%
  tally()%>%
  filter(n>4)%>%
  arrange(-n)
# mobility_metric                    n
# <fct>                          <int>
# 1 "Religious\norganizations"         8
# 2 "Elementary and\nhigh schools"     7

all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, start_week) %>% 
  group_by(mobility_metric)%>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                   min   max
# <fct>                           <dbl> <dbl>
# 1 "Religious\norganizations"      0.767 0.872
# 2 "Elementary and\nhigh schools"  0.703 0.823

all_path_ccf %>%
  filter(pathogen == "Adenovirus") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf)

# negative correlations
all_path_ccf %>%
  filter(pathogen == "Adenovirus" & mobility_metric == "Oxford Stringency\nIndex") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf < 0) %>%
  filter(start_week > as.Date("2020-06-01") & start_week < as.Date("2021-01-01")) %>%
  distinct() %>%
  arrange(mobility_metric, -obs_max_ccf)
## no sig. correlations with OSI

########################################################################################
## Enveloped virus cross-correlations between Rt and mobility (pandemic rebound during 2021)
########################################################################################

########################################################################################
## RSV B, summer 2021
## statistically significant cross-correlations between Rt and mobility
########################################################################################
## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric)%>%
  tally()%>%
  filter(n>3)%>%
  arrange(-n)
# mobility_metric                      n
# # <fct>                            <int>
# 1 "Religious\norganizations"           9
# 2 "Between-neighborhood\nmovement"     8
# 3 "Child daycare"                      7
# 4 "Elementary and\nhigh schools"       7
# 5 "% not wearing masks"                5
# 6 "Influx out-of-state\nvisitors"      4

keep = all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric)%>%
  tally()%>%
  filter(n>3)%>%
  arrange(-n) %>%
  pull(mobility_metric)%>%unique()%>%as.character()

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & mobility_metric != "Oxford Stringency\nIndex") %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-08-01")) %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  filter(obs_max_ccf>0.5 & mobility_metric %in% keep)%>%
  arrange(mobility_metric, obs_max_ccf)%>%
  group_by(mobility_metric)%>%
  summarize(
    min = min(obs_max_ccf),
    max = max(obs_max_ccf)
  ) %>%
  arrange(-max)
# mobility_metric                    min   max
# <fct>                            <dbl> <dbl>
# 1 "% not wearing masks"            0.862 0.912
# 2 "Child daycare"                  0.732 0.894
# 3 "Elementary and\nhigh schools"   0.716 0.872
# 4 "Influx out-of-state\nvisitors"  0.726 0.837
# 5 "Between-neighborhood\nmovement" 0.699 0.812
# 6 "Religious\norganizations"       0.697 0.775

all_path_ccf %>%
  filter(mobility_metric == "Between-neighborhood\nmovement" & pathogen == "RSV B") %>%
  filter(start_week > as.Date("2021-03-01") & start_week < as.Date("2021-09-01")) %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf)

## max correlation values
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
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  filter(mobility_metric != "Oxford Stringency\nIndex") %>%
  distinct() %>%
  filter(pathogen == "RSV B") %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  group_by(mobility_metric) %>%
  slice_min(obs_max_ccf) %>%
  arrange(obs_max_ccf)

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
# 1 hPIV 3 + 4 "Elementary and\nhigh schools"      10
# 2 hPIV 3 + 4 "Between-neighborhood\nmovement"     7
# 3 hPIV 3 + 4 "Restaurants"                        7
# 4 hPIV 3 + 4 "% not wearing masks"                6
# 5 hPIV 3 + 4 "Child daycare"                      5

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-05-01")) %>%
  distinct() %>%
  filter(pathogen == "hPIV 3 + 4") %>%
  arrange(mobility_metric, start_week, obs_max_ccf)%>%
  group_by(mobility_metric)%>%
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
# 1 "Between-neighborhood\nmovement"       10
# 2 "% Devices\nleaving home"              10
# 3 "Within-neighborhood\nmovement"         7
# 4 "Infux visitors\nother WA counties"     7
# 5 "Child daycare"                         7
# 6 "Influx out-of-state\nvisitors"         6
# 7 "% not wearing masks"                   5

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
  filter(start_week > as.Date("2020-12-01") & start_week < as.Date("2021-09-01")) %>%
  distinct() %>%
  filter(pathogen == "hCoV HKU1 + NL63") %>%
  arrange(mobility_metric, -obs_max_ccf)

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
  group_by(pathogen,mobility_metric) %>%
  tally() %>%
  arrange(-n) %>%
  filter(n>2) %>%
  arrange(pathogen)%>%
  print(n = 50)

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

# adenovirus
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.5) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  filter(!(mobility_metric %in% c("Oxford Stringency\nIndex"))) %>%
  filter(pathogen == "Adenovirus") %>%
  distinct() %>%
  arrange(pathogen, mobility_metric) %>%
  group_by(mobility_metric, obs_max_ccf_lag) %>%
  tally() %>%
  arrange(-n) %>%
  print(n = 50)

## positive correlations
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  arrange(pathogen, mobility_metric)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  arrange(pathogen, mobility_metric) %>%
  distinct(pathogen, mobility_metric, obs_max_ccf_lag)

# % devices leaving home
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0.5) %>%
  filter(start_week > as.Date("2021-11-01") & start_week < as.Date("2022-02-01")) %>%
  distinct() %>%
  filter(pathogen != "SARS-CoV-2") %>%
  filter(mobility_metric == "% Devices\nleaving home") %>%
  arrange(pathogen, mobility_metric)

