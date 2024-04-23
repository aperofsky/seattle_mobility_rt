########################################################################################
## Moving window cross correlations between Rt and mobility
## Figure S7: February 2019 snowstorm; plotting correlations for each window
########################################################################################

library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
library(gghighlight)

########################################################################################
## Import data
########################################################################################

dir <- "4_Block_Bootstrap_Mobility_vs_Rt/1_snowstorm_rolling_window_block_bootstrap/biowulf_block_bootstrap_1mo_rolling_window_output_spearman/"

## adenovirus
load(paste0(dir, "adeno_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## seasonal CoV
load(paste0(dir, "scov_229E_OC43_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

load(paste0(dir, "scov_HKU1_NL63_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## hmpv
load(paste0(dir, "hmpv_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## HPIV
load(paste0(dir, "hpiv_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## rhino
load(paste0(dir, "rhino_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## entero
load(paste0(dir, "entero_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

## rsv a
load(paste0(dir, "rsv_a_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

# rsv b
load(paste0(dir, "rsv_b_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

### h1n1
load(paste0(dir, "h1n1_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

# h3n2
load(paste0(dir, "h3n2_snowstorm_sliding_window_actual_and_null_output_spearman.RData"))

all_path_ccf <- bind_rows(
  adeno_block_bs,
  entero_block_bs,
  h1n1_block_bs,
  h3n2_block_bs,
  hmpv_block_bs,
  hpiv_block_bs,
  rhino_block_bs,
  rsva_block_bs,
  rsvb_block_bs,
  scov_229E_OC43_block_bs,
  scov_HKU1_NL63_block_bs
)
head(all_path_ccf)

all_path_ccf$pathogen <- as.factor(all_path_ccf$pathogen)
all_path_ccf$month <- zoo::as.yearmon(all_path_ccf$start_date, "%Y %m")
levels(all_path_ccf$pathogen)
levels(all_path_ccf$mobility_metric)
unique(all_path_ccf$pathogen)

all_path_ccf$pathogen <- factor(all_path_ccf$pathogen, levels = c(
  "h1n1", "h3n2", "rsv_a",
  "rsv_b", "hmpv",
  "scov_229E_OC43",
  "scov_HKU1_NL63",
  "hpiv_3_4",
  "rhino", "entero", "adeno"
))

levels(all_path_ccf$pathogen)
levels(all_path_ccf$pathogen) <- c(
  "Influenza A/H1N1",
  "Influenza A/H3N2",
  "RSV A",
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  "hPIV 3 + 4",
  "Rhinovirus",
  "Enterovirus",
  "Adenovirus"
)

levels(all_path_ccf$mobility_metric)
levels(all_path_ccf$pathogen)
unique(all_path_ccf$mobility_metric)
all_path_ccf$mobility_metric <- factor(all_path_ccf$mobility_metric,
  levels = c(
    "within_neighborhood_movement",
    "within_city_movement",
    "within_state_movement",
    "out_of_state_movement",
    "full_service_restaurants",
    "performing_arts_or_sports_events",
    "groceries_and_pharmacies",
    "transit",
    "religious_orgs", "child_day_care",
    "elementary_and_secondary_schools",
    "colleges"
  )
)
levels(all_path_ccf$mobility_metric)

levels(all_path_ccf$mobility_metric) <- c(
  "Within-neighborhood\nmovement",
  "Between-neighborhood\nmovement",
  "Infux visitors\nother WA counties",
  "Influx out-of-state\nvisitors",
  "Restaurants",
  "Performing arts or\nsports events",
  "Groceries and\npharmacies",
  "Transit",
  "Religious\norganizations",
  "Child daycare",
  "Elementary and\nhigh schools",
  "Colleges"
)

all_endemic_results_avg <- read_csv("2_Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv") %>%
  filter(date >= as.Date("2019-01-15") & date < as.Date("2019-08-01"))
unique(all_endemic_results_avg$organism)
all_endemic_results_avg$organism <- as.factor(all_endemic_results_avg$organism)

levels(all_endemic_results_avg$organism)
levels(all_endemic_results_avg$organism) <- c(
  "Adenovirus",
  "Enterovirus",
  "Human metapneumovirus",
  # "Human parainfluenza 1 + 2", #not circulating
  "Human parainfluenza 3 + 4",
  "Influenza A/H1N1",
  "Influenza A/H3N2",
  # "Influenza B",# not circulating
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Rhinovirus",
  "Seasonal CoV 229E + OC43",
  "Seasonal CoV HKU1 + NL63"
)

unique(all_endemic_results_avg$organism)

########################################################################################
## Figure S7. Plot snowstorm block bootstrap
########################################################################################

## only include colors for pathogens circulating during Feb 2019
color_vec <- c(
  "#332288",
  "#6699CC",
  # "#88CCEE",
  "#44AA99",
  "#117733",
  "#999933",
  "#DDCC77",
  "#661100",
  # "#CC6677",
  "#AA4466",
  "#882255",
  "#993377",
  "#AA4499"
)

unique(all_endemic_results_avg$organism)
all_endemic_results_avg$organism <- factor(all_endemic_results_avg$organism,
  levels = c(
    "Influenza A/H1N1",
    "Influenza A/H3N2",
    "Respiratory syncytial virus (RSV) A",
    "Respiratory syncytial virus (RSV) B",
    "Human metapneumovirus",
    "Seasonal CoV 229E + OC43",
    "Seasonal CoV HKU1 + NL63",
    "Human parainfluenza 3 + 4",
    "Rhinovirus", "Enterovirus", "Adenovirus"
  )
)
levels(all_endemic_results_avg$organism)
range(all_path_ccf$start_date)
df1 <- all_path_ccf %>% filter(start_date >= as.Date("2019-01-20") & start_date < as.Date("2019-03-01"))

all_mob_plot <- ggplot(
  data = df1 %>% filter(obs_max_ccf_lag < 1),
  aes(x = start_date, group = sig, y = as.numeric(obs_max_ccf), label = obs_max_ccf_lag)
) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0.02, 0.02)) +
  annotate("rect",
    xmin = as.Date("2019-02-03"), # snowstorm period
    xmax = as.Date("2019-02-15"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_point(pch = 21, size = 5, aes(fill = as.numeric(obs_max_ccf_lag), alpha = sig)) +
  geom_text(aes(alpha = sig), hjust = 0.5, vjust = 0.5, size = 3, color = "white") +
  facet_grid(mobility_metric ~ pathogen) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.01") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c(begin = 0, end = 0.75, breaks = c(0, -7, -14), labels = c(0, -7, -14)) +
  xlab("1-Month Window Mid-Point") +
  labs(fill = "Optimal Lag\n(Days)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white"),
    legend.position = "bottom",
    strip.text.y.right = element_text(angle = 0)
  )
all_mob_plot

levels(all_endemic_results_avg$organism) <- levels(all_path_ccf$pathogen)

keep <- c(
  "Influenza A/H1N1",
  "Influenza A/H3N2",
  "RSV A",
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  "hPIV 3 + 4",
  "Rhinovirus", "Enterovirus", "Adenovirus"
)

winter_2019_rt <- ggplot() +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  annotate("rect",
    xmin = as.Date("2019-02-03"), # snowstorm period
    xmax = as.Date("2019-02-15"),
    ymin = -Inf,
    ymax = Inf, fill = "blue", alpha = 0.1
  ) +
  geom_line(
    data = all_endemic_results_avg %>%
      filter(organism %in% keep) %>%
      filter(date >= as.Date("2019-01-20") & date < as.Date("2019-03-01")) %>%
      dplyr::select(date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = all_endemic_results_avg %>%
      filter(date >= as.Date("2019-01-20") & date < as.Date("2019-03-01")) %>%
      filter(organism %in% keep) %>%
      dplyr::select(date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = date, ymin = lower, ymax = upper, fill = organism), alpha = 0.5
  ) +
  theme_bw(base_size = 12) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", expand = c(0.02, 0.02)) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec) +
  facet_wrap(~organism, nrow = 1) +
  ylab("Rt") +
  xlab("Date")
winter_2019_rt

rt_top <- plot_grid(NULL, winter_2019_rt, NULL, rel_widths = c(0.15, 9, 1), nrow = 1)
rt_mob_winter_2019 <- plot_grid(rt_top, all_mob_plot, nrow = 2, rel_heights = c(1.5, 12), labels = "AUTO")
rt_mob_winter_2019
save_plot(rt_mob_winter_2019, file = "figures/fig_s7_block_bootstrap_snowstorm_and_rt_spearman_neg_lags.png", base_width = 22, base_height = 16)

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0 & month == "Feb 2019") %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  filter(start_date >= as.Date("2019-02-01") & start_date <= as.Date("2019-02-19")) %>%
  # arrange(-obs_max_ccf,pathogen,mobility_metric)%>%
  group_by(pathogen, mobility_metric) %>%
  summarize(mean = mean(obs_max_ccf)) %>%
  arrange(-mean, pathogen, mobility_metric) %>%
  print(n = 50)

all_path_ccf %>%
  filter(pathogen %in% c("RSV A", "RSV B", "Adenovirus", "Enterovirus") & month == "Feb 2019") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & month == "Feb 2019" & obs_max_ccf > 0) %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  filter(start_date >= as.Date("2019-02-01") & start_date <= as.Date("2019-02-19")) %>%
  # arrange(-obs_max_ccf,pathogen,mobility_metric)%>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(pathogen, -n) %>%
  slice_max(n, n = 5) %>%
  print(n = 40)

mob_keep <- all_path_ccf %>%
  filter(pathogen %in% c("RSV A", "RSV B", "Adenovirus", "Enterovirus") & month == "Feb 2019") %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & month == "Feb 2019" & obs_max_ccf > 0) %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  # filter(start_date >= as.Date("2019-02-01") & start_date <= as.Date("2019-02-19"))%>%
  # filter(start_date >= as.Date("2019-02-01") & start_date <= as.Date("2019-02-19"))%>%
  # arrange(-obs_max_ccf,pathogen,mobility_metric)%>%
  group_by(pathogen, mobility_metric) %>%
  tally() %>%
  arrange(pathogen, -n) %>%
  slice_max(n, n = 5) %>%
  pull(mobility_metric) %>%
  unique() %>%
  as.character()

mob_keep
# [1] "Elementary and\nhigh schools"      "Infux visitors\nother WA counties" "Transit"
# [4] "Child daycare"                     "Colleges"                          "Restaurants"
# [7] "Religious\norganizations"          "Influx out-of-state\nvisitors"

## avg ccf in Feb 2019 for metrics with strongest relationships with RSV and AdV
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1 & obs_max_ccf > 0) %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  arrange(-obs_max_ccf, pathogen, mobility_metric) %>%
  filter(pathogen %in% c("RSV A", "RSV B", "Adenovirus", "Enterovirus") & month == "Feb 2019") %>%
  filter(mobility_metric %in% mob_keep) %>%
  filter(mobility_metric != "Transit") %>%
  # filter(start_date >= as.Date("2019-02-03") & start_date <= as.Date("2019-02-19"))%>%
  group_by(pathogen, mobility_metric) %>%
  summarize(mean = mean(obs_max_ccf)) %>%
  group_by(pathogen) %>%
  # slice_max(mean, n = 6) %>%
  arrange(mobility_metric, -mean) %>%
  print(n = 40)
# pathogen    mobility_metric                      mean
# 1 Adenovirus  "Infux visitors\nother WA counties" 0.966
# 2 Enterovirus "Infux visitors\nother WA counties" 0.929
# 3 RSV B       "Infux visitors\nother WA counties" 0.850
# 4 RSV A       "Infux visitors\nother WA counties" 0.833
# 5 Enterovirus "Influx out-of-state\nvisitors"     0.977
# 6 Adenovirus  "Influx out-of-state\nvisitors"     0.949
# 7 RSV B       "Influx out-of-state\nvisitors"     0.852
# 8 RSV A       "Influx out-of-state\nvisitors"     0.823
# 9 Adenovirus  "Restaurants"                       0.900
# 10 Enterovirus "Restaurants"                       0.852
# 11 RSV B       "Restaurants"                       0.811
# 12 RSV A       "Restaurants"                       0.789
# 13 Adenovirus  "Religious\norganizations"          0.968
# 14 Enterovirus "Religious\norganizations"          0.919
# 15 RSV B       "Religious\norganizations"          0.860
# 16 RSV A       "Religious\norganizations"          0.844
# 17 Adenovirus  "Child daycare"                     0.945
# 18 Enterovirus "Child daycare"                     0.909
# 19 RSV B       "Child daycare"                     0.869
# 20 RSV A       "Child daycare"                     0.843
# 21 RSV B       "Elementary and\nhigh schools"      0.928
# 22 RSV A       "Elementary and\nhigh schools"      0.922
# 23 Adenovirus  "Elementary and\nhigh schools"      0.866
# 24 Adenovirus  "Colleges"                          0.938
# 25 Enterovirus "Colleges"                          0.903
# 26 RSV B       "Colleges"                          0.829
# 27 RSV A       "Colleges"                          0.806

all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1) %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  filter(pathogen %in% c("hPIV 3 + 4", "hMPV", "Rhinovirus") & month == "Feb 2019") %>%
  group_by(pathogen, mobility_metric) %>%
  summarize(mean = mean(obs_max_ccf)) %>%
  group_by(pathogen) %>%
  slice_max(mean, n = 5) %>%
  print(n = 30)
# pathogen   mobility_metric                       mean
# <fct>      <fct>                                <dbl>
# 1 hMPV       "Transit"                           -0.750
# 2 hMPV       "Colleges"                          -0.758
# 3 hMPV       "Restaurants"                       -0.769
# 4 hMPV       "Child daycare"                     -0.781
# 5 hMPV       "Groceries and\npharmacies"         -0.790
# 6 hPIV 3 + 4 "Within-neighborhood\nmovement"      0.217
# 7 hPIV 3 + 4 "Influx out-of-state\nvisitors"      0.207
# 8 hPIV 3 + 4 "Elementary and\nhigh schools"       0.159
# 9 hPIV 3 + 4 "Colleges"                           0.131
# 10 hPIV 3 + 4 "Child daycare"                      0.0791
# 11 Rhinovirus "Infux visitors\nother WA counties"  0.117
# 12 Rhinovirus "Religious\norganizations"           0.108
# 13 Rhinovirus "Child daycare"                      0.0808
# 14 Rhinovirus "Colleges"                           0.0666
# 15 Rhinovirus "Influx out-of-state\nvisitors"      0.0647
