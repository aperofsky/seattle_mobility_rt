########################################################################################
## Moving window cross correlations between Rt and mobility
## Figure S7: February 2019 snowstorm; plotting correlations for each window
## Mean cross-correlations between Rt and mobility during Feb 2019
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

dir <- "Block_Bootstrap_Mobility_vs_Rt/biowulf_block_bootstrap_snowstorm_rolling_window_exp_decay/biowulf_block_bootstrap_1mo_rolling_window_exp_decay_output/"

## adenovirus
load(paste0(dir, "adeno_snowstorm_sliding_window_actual_and_null_output.RData"))

## seasonal CoV
load(paste0(dir, "scov_229E_OC43_snowstorm_sliding_window_actual_and_null_output.RData"))

load(paste0(dir, "scov_HKU1_NL63_snowstorm_sliding_window_actual_and_null_output.RData"))

## hmpv
load(paste0(dir, "hmpv_snowstorm_sliding_window_actual_and_null_output.RData"))

## HPIV
load(paste0(dir, "hpiv_snowstorm_sliding_window_actual_and_null_output.RData"))

## rhino
load(paste0(dir, "rhino_snowstorm_sliding_window_actual_and_null_output.RData"))

## rsv a
load(paste0(dir, "rsv_a_snowstorm_sliding_window_actual_and_null_output.RData"))

# rsv b
load(paste0(dir, "rsv_b_snowstorm_sliding_window_actual_and_null_output.RData"))

### h1n1
load(paste0(dir, "h1n1_snowstorm_sliding_window_actual_and_null_output.RData"))

# h3n2
load(paste0(dir, "h3n2_snowstorm_sliding_window_actual_and_null_output.RData"))

all_path_ccf <- bind_rows(
  adeno_block_bs,
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
all_path_ccf$pathogen <- factor(all_path_ccf$pathogen, levels = c(
  "h1n1", "h3n2", "rsv_a",
  "rsv_b", "hmpv",
  "scov_229E_OC43",
  "scov_HKU1_NL63",
  "hpiv_3_4",
  "rhino", "adeno"
))

levels(all_path_ccf$pathogen) <- c(
  "Influenza A/H1N1", "Influenza A/H3N2",
  "RSV A", "RSV B", "hMPV",
  "hCoV 229E + OC43", "hCoV HKU1 + NL63",
  "hPIV 3 + 4", "Rhinovirus",
  "Adenovirus"
)

levels(all_path_ccf$mobility_metric)
levels(all_path_ccf$pathogen)
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
  "Religious\norganizations", "Child daycare",
  "Elementary and\nhigh schools",
  "Colleges"
)


all_endemic_results_avg <- read_csv("Block_Bootstrap_Mobility_vs_Rt/rt_all_pathogens_15day_mv_avg.csv") %>%
  filter(date >= as.Date("2019-01-15") & date < as.Date("2019-08-01"))
unique(all_endemic_results_avg$organism)
all_endemic_results_avg$organism <- as.factor(all_endemic_results_avg$organism)

levels(all_endemic_results_avg$organism) <- c(
  "Adenovirus",
  "Human metapneumovirus",
  # "Human parainfluenza 1 + 2",
  "Human parainfluenza 3 + 4",
  "Influenza A/H1N1", "Influenza A/H3N2",
  # "Influenza B",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Rhinovirus",
  "Seasonal CoV 229E + OC43",
  "Seasonal CoV HKU1 + NL63"
)



unique(all_endemic_results_avg$organism)
# color_vec <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

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
  # "#661100",
  "#CC6677",
  "#AA4466",
  "#882255",
  "#AA4499"
  # "#ABB065"
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
    "Rhinovirus", "Adenovirus"
  )
)
levels(all_endemic_results_avg$organism)
range(all_path_ccf$start_date)
df1 <- all_path_ccf %>% filter(start_date >= as.Date("2019-01-15") & start_date < as.Date("2019-03-01"))

names(df1)
all_mob_plot <- ggplot(
  data = df1,
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
  scale_alpha_manual(values = c(0.2, 0.9), name = "p(perm)<0.05") +
  theme_bw(base_size = 16) +
  ylab("Cross-Correlation Coefficient") +
  scale_fill_viridis_c() +
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
# save_plot(all_mob_plot,file="figures/block_bootstrap_snowstorm_period_1mo_moving_window.png",base_width = 26,base_height = 14)

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
  "Rhinovirus", "Adenovirus"
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
      filter(date >= as.Date("2019-01-15") & date < as.Date("2019-03-01")) %>%
      dplyr::select(date, median, organism) %>%
      distinct() %>% droplevels(),
    aes(x = date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = all_endemic_results_avg %>%
      filter(date >= as.Date("2019-01-15") & date < as.Date("2019-03-01")) %>%
      filter(organism %in% keep) %>%
      dplyr::select(date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% droplevels(),
    aes(x = date, ymin = lower, ymax = upper, fill = organism), alpha = 0.8
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
rt_mob_fall_2019 <- plot_grid(rt_top, all_mob_plot, nrow = 2, rel_heights = c(1.5, 12), labels = "AUTO")
rt_mob_fall_2019
save_plot(rt_mob_fall_2019, file = "figures/fig_s7_block_bootstrap_snowstorm_and_rt.png", base_width = 22, base_height = 16)

## avg ccf in Feb 2019 for metrics with strongest relationships with RSV and AdV
all_path_ccf %>%
  filter(sig == "yes" & obs_max_ccf_lag < 1) %>% # only include dates when there's a significant correlation and synchronous/leading relationship
  distinct() %>%
  arrange(mobility_metric, obs_max_ccf) %>%
  filter(pathogen %in% c("RSV A", "RSV B", "Adenovirus") & month == "Feb 2019") %>%
  group_by(pathogen, mobility_metric) %>%
  summarize(mean = mean(obs_max_ccf)) %>%
  group_by(pathogen) %>%
  slice_max(mean, n = 5) %>%
  arrange(pathogen,-mean)%>%
  print(n = 20)
# pathogen   mobility_metric                      mean
# <fct>      <fct>                               <dbl>
# 1 RSV A      "Between-neighborhood\nmovement"    0.934
# 2 RSV A      "Elementary and\nhigh schools"      0.918
# 3 RSV A      "Influx out-of-state\nvisitors"     0.913
# 4 RSV A      "Child daycare"                     0.911
# 5 RSV A      "Religious\norganizations"          0.906
# 6 RSV B      "Elementary and\nhigh schools"      0.943
# 7 RSV B      "Between-neighborhood\nmovement"    0.943
# 8 RSV B      "Influx out-of-state\nvisitors"     0.932
# 9 RSV B      "Religious\norganizations"          0.924
# 10 RSV B      "Child daycare"                     0.907
# 11 Adenovirus "Religious\norganizations"          0.981
# 12 Adenovirus "Infux visitors\nother WA counties" 0.973
# 13 Adenovirus "Child daycare"                     0.972
# 14 Adenovirus "Influx out-of-state\nvisitors"     0.962
# 15 Adenovirus "Within-neighborhood\nmovement"     0.950

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
# 1 hMPV       "Colleges"                          -0.859
# 2 hMPV       "Within-neighborhood\nmovement"     -0.866
# 3 hMPV       "Restaurants"                       -0.872
# 4 hMPV       "Infux visitors\nother WA counties" -0.874
# 5 hMPV       "Transit"                           -0.875
# 6 hPIV 3 + 4 "Between-neighborhood\nmovement"     0.582
# 7 hPIV 3 + 4 "Child daycare"                      0.391
# 8 hPIV 3 + 4 "Influx out-of-state\nvisitors"      0.240
# 9 hPIV 3 + 4 "Colleges"                           0.189
# 10 hPIV 3 + 4 "Infux visitors\nother WA counties"  0.138
# 11 Rhinovirus "Influx out-of-state\nvisitors"      0.956
# 12 Rhinovirus "Religious\norganizations"           0.895
# 13 Rhinovirus "Child daycare"                      0.816
# 14 Rhinovirus "Infux visitors\nother WA counties"  0.702
# 15 Rhinovirus "Between-neighborhood\nmovement"     0.690
