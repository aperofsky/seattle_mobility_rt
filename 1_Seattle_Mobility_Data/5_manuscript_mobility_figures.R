## start here if you don't have access to the raw SafeGraph data
####################################################################################
## Mobility figures
## Figure S4: % change in baseline for individual POI categories
## Figure 2: % change in baseline for large scale movements and different POI categories; % staying home and % masking
## Figure S17: Omicron period % change in baseline for large scale movements and different POI categories; % staying home and % masking
####################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
library(jcolors)
library(ggpubfigs)
library(pals)
library(zoo)

combined_mob <- read_rds("1_Seattle_Mobility_Data/mobility_data/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)
range(combined_mob$date)
names(combined_mob)

## add in time periods for covid waves
combined_long <- combined_mob %>%
  dplyr::select(
    -mask_wearing_final,
    # -fb_stay_put_custom
  ) %>%
  mutate(time_period = case_when(
    date < as.Date("2020-01-01") ~ "pre-pandemic",
    date >= as.Date("2020-01-01") & date < as.Date("2021-03-01") ~ "ancestral",
    date >= as.Date("2021-03-01") & date < as.Date("2021-11-01") ~ "delta",
    date >= as.Date("2021-11-01") ~ "omicron"
  )) %>%
  tidyr::pivot_longer(cols = c(within_neighborhood_movement:fb_stay_put_custom), names_to = "mobility_metric", values_to = "value")

combined_long$mobility_metric <- as.factor(combined_long$mobility_metric)
levels(combined_long$mobility_metric)
levels(combined_long$mobility_metric) <- c(
  "child day care",
  "colleges",
  "elementary and high schools",
  "% devices staying home",
  "full service restaurants",
  "groceries and pharmacies",
  "influx out-of-state visitors",
  "performing arts and sports events",
  "religious organizations",
  "transit",
  "between-neighborhood movement",
  "within-neighborhood movement",
  "influx other WA counties"
)

combined_long$mobility_metric <- factor(combined_long$mobility_metric, levels = c(
  "within-neighborhood movement",
  "between-neighborhood movement",
  "influx other WA counties",
  "influx out-of-state visitors",
  "% devices staying home",
  "full service restaurants",
  "groceries and pharmacies",
  "transit",
  "performing arts and sports events",
  "religious organizations",
  "child day care",
  "elementary and high schools",
  "colleges"
))
levels(combined_long$mobility_metric)

unique(combined_long$epi_date)
####################################################################################
## declines in mobility during snowstorm
####################################################################################
combined_long %>%
  filter(!is.na(value)) %>%
  dplyr::select(date, mobility_metric, value) %>%
  filter(date <= as.Date("2019-02-24") & date >= as.Date("2019-02-03")) %>%
  arrange(mobility_metric, date) %>%
  group_by(mobility_metric) %>%
  slice_min(value, n = 1) %>%
  slice_head(n = 1) %>%
  arrange(value) %>%
  print(n = 15)
# date       mobility_metric                    value
# <date>     <fct>                              <dbl>
# 1 2019-02-05 elementary and high schools       -84.0
# 2 2019-02-11 performing arts and sports events -81.6
# 3 2019-02-09 colleges                          -80.9
# 4 2019-02-17 transit                           -75.6
# 5 2019-02-04 child day care                    -60.3
# 6 2019-02-04 influx out-of-state visitors      -53.0
# 7 2019-02-04 groceries and pharmacies          -49.1
# 8 2019-02-10 full service restaurants          -48.2
# 9 2019-02-04 influx other WA counties          -38.2
# 10 2019-02-04 between-neighborhood movement     -32.2
# 11 2019-02-18 within-neighborhood movement      -29.3
# 12 2019-02-09 religious organizations           -26.8

####################################################################################
## declines in mobility after State of Emergency and by the time of business closures on March 16
####################################################################################

combined_long %>%
  filter(!is.na(value)) %>%
  dplyr::select(date, mobility_metric, value) %>%
  filter(date <= as.Date("2020-03-17") & date >= as.Date("2020-02-27")) %>%
  arrange(mobility_metric) %>%
  group_by(mobility_metric) %>%
  slice_min(value, n = 1) %>%
  slice_head(n = 1) %>%
  arrange(value) %>%
  print(n = 15)
# date       mobility_metric                   value
# <date>     <fct>                             <dbl>
# 1 2020-03-15 transit                           -94.6
# 2 2020-03-15 performing arts and sports events -88.8
# 3 2020-03-14 elementary and high schools       -86.1
# 4 2020-03-15 colleges                          -82.9
# 5 2020-03-16 influx out-of-state visitors      -68.9
# 6 2020-03-15 child day care                    -67.3
# 7 2020-03-16 within-neighborhood movement      -62.7
# 8 2020-03-16 between-neighborhood movement     -62.2
# 9 2020-03-17 full service restaurants          -59.6
# 10 2020-03-14 religious organizations           -54.6
# 11 2020-03-16 influx other WA counties          -43.3
# 12 2020-03-08 groceries and pharmacies          -28.8

####################################################################################
## declines in mobility after SAH orders
####################################################################################
combined_long %>%
  filter(!is.na(value)) %>%
  dplyr::select(date, mobility_metric, value) %>%
  filter(date <= as.Date("2020-06-03") & date >= as.Date("2020-03-23")) %>%
  arrange(mobility_metric) %>%
  group_by(mobility_metric) %>%
  slice_min(value, n = 1) %>%
  slice_head(n = 1) %>%
  arrange(value) %>%
  print(n = 15)
# date       mobility_metric                   value
# <date>     <fct>                             <dbl>
# 1 2020-05-09 transit                           -95.4
# 2 2020-05-25 performing arts and sports events -92.3
# 3 2020-05-10 colleges                          -90.4
# 4 2020-05-25 elementary and high schools       -90.4
# 5 2020-04-06 influx out-of-state visitors      -84.3
# 6 2020-04-12 full service restaurants          -78.4
# 7 2020-04-12 child day care                    -76.0
# 8 2020-04-20 between-neighborhood movement     -73.7
# 9 2020-04-20 within-neighborhood movement      -70.6
# 10 2020-05-04 religious organizations           -65.1
# 11 2020-04-06 influx other WA counties          -55.2
# 12 2020-04-12 groceries and pharmacies          -47.6

####################################################################################
## max % staying at home during SAH orders
####################################################################################

combined_long %>%
  filter(!is.na(value)) %>%
  dplyr::select(date, mobility_metric, value) %>%
  filter(date <= as.Date("2020-06-03") & date >= as.Date("2020-03-23") & mobility_metric == "% devices staying home") %>%
  slice_max(value, n = 1)
# date       mobility_metric                    value
# <date>     <fct>                              <dbl>
# 13 2020-03-29 % devices staying home             53.4
####################################################################################
## Figure S4: % change in baseline for individual POI categories
####################################################################################
## use 2-week moving average to reduce noise
## 2-week moving average to reduce noise
combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = NA)))

combined_long <- combined_mob_avg %>%
  dplyr::select(
    -mask_wearing_final,
    -fb_stay_put_custom
  ) %>%
  mutate(time_period = case_when(
    date < as.Date("2020-01-01") ~ "pre-pandemic",
    date >= as.Date("2020-01-01") & date < as.Date("2021-03-01") ~ "ancestral",
    date >= as.Date("2021-03-01") & date < as.Date("2021-11-01") ~ "delta",
    date >= as.Date("2021-11-01") ~ "omicron"
  )) %>%
  tidyr::pivot_longer(cols = c(within_neighborhood_movement:colleges), names_to = "mobility_metric", values_to = "value")

combined_long$mobility_metric <- as.factor(combined_long$mobility_metric)
levels(combined_long$mobility_metric)
levels(combined_long$mobility_metric) <- c(
  "child day care",
  "colleges",
  "elementary and high schools",
  "full service restaurants",
  "groceries and pharmacies",
  "influx out-of-state visitors",
  "performing arts and sports events",
  "religious organizations",
  "transit",
  "between-neighborhood movement",
  "within-neighborhood movement",
  "influx other WA counties"
)

combined_long$mobility_metric <- factor(combined_long$mobility_metric, levels = c(
  "within-neighborhood movement",
  "between-neighborhood movement",
  "influx other WA counties",
  "influx out-of-state visitors",
  "full service restaurants",
  "groceries and pharmacies",
  "transit",
  "performing arts and sports events",
  "religious organizations",
  "child day care",
  "elementary and high schools",
  "colleges"
))
levels(combined_long$mobility_metric)


col_vec <- c("blue", "#843C39", "#DE9ED6", "darkgreen", "#f46d43", "#9C9EDE", "#E7969C", "#A55194", "#6DDE88", "#5254A3", "#BD9E39", "#637939", "skyblue")

p <- ggplot() +
  geom_rect(
    data = data.frame(x = 0, y = 0), aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_rect(
    data = data.frame(x = 0, y = 0), aes(xmin = as.Date("2020-03-16"), xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "orange"
  ) +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen") +
  geom_line(
    data = combined_long %>% filter(date >= as.Date("2018-11-01") & date < as.Date("2022-07-01")),
    aes(x = date, y = value, color = mobility_metric), lwd = 1
  ) +
  facet_wrap(~mobility_metric, scales = "free_y", dir = "v") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_color_manual(values = col_vec) +
  scale_x_date(expand = c(0, 0)) +
  ylab("Percent change from baseline") +
  ylim(-90, 90)
p
save_plot(p, filename = "1_Seattle_Mobility_Data/mobility_figures/fig_s4_individual_mobility_metrics.png", base_width = 16, base_height = 10)
save_plot(p, filename = "figures/fig_s4_individual_mobility_metrics.png", base_width = 16, base_height = 10)

####################################################################################
## Large Scale Population Movements: Figure 2A and Figure 16A
####################################################################################

flow_plot <- ggplot() +
  geom_rect(aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_rect(aes(xmin = as.Date("2020-03-16"), xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "orange"
  ) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(date >= as.Date("2018-11-18") & date < as.Date("2022-07-01")) %>%
      filter(mobility_metric %in% c(
        "influx out-of-state visitors",
        "influx other WA counties",
        "between-neighborhood movement",
        "within-neighborhood movement"
      )),
    aes(x = date, y = value, color = mobility_metric), lty = "solid", lwd = 1, alpha = 0.8
  ) +
  theme_bw(base_size = 16) +
  xlab("Date") +
  ylab("Percent Change from Baseline") +
  scale_color_manual(
    values = friendly_pal("muted_nine"),
    labels = c("Within-Neighborhood", "Between-Neighborhood", "Within-State", "Out-of-State"),
    name = NULL
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b-%y") +
  ggtitle("Within-City Movement and In-Flows") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ylim(c(-90, 90))
flow_plot

## declines during omicron
omicron_flow_plot <- ggplot() +
  geom_rect(aes(xmin = as.Date("2021-11-01"), xmax = as.Date("2022-02-01"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(mobility_metric %in% c(
        "influx out-of-state visitors",
        "influx other WA counties",
        "between-neighborhood movement",
        "within-neighborhood movement"
      ) & date > as.Date("2021-09-01") & date < as.Date("2022-04-01")),
    aes(x = date, y = value, color = mobility_metric), lty = "solid", lwd = 1, alpha = 0.8
  ) +
  theme_bw(base_size = 16) +
  xlab("Date") +
  ylab("Percent Change from Baseline") +
  scale_color_manual(
    values = friendly_pal("muted_nine"),
    labels = c("Within-Neighborhood", "Between-Neighborhood", "Within-State", "Out-of-State"), name = NULL
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 months", date_labels = "%b-%y", limits = c(as.Date("2021-09-01"), as.Date("2022-04-01"))) +
  ggtitle("Within-City Movement and In-Flows") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ylim(c(-60, 60))

####################################################################################
## Foot traffic to different categories of POIs: Figure 2B and Figure 16B
####################################################################################

poi_plot <- ggplot() +
  geom_rect(aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_rect(aes(xmin = as.Date("2020-03-16"), xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "orange"
  ) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(date >= as.Date("2018-11-18") & date < as.Date("2022-07-01")) %>%
      filter(mobility_metric %in% c(
        "transit", "elementary and high schools", "colleges",
        "full service restaurants", "religious organizations",
        "groceries and pharmacies"
      )),
    aes(x = date, y = value, color = mobility_metric), lwd = 0.7, alpha = 0.8
  ) +
  theme_bw(base_size = 16) +
  xlab("Date") +
  ylab("Percent Change from Baseline") +
  scale_color_manual(
    values = c("#a45ccc", "#b88f06", "#4faa45", "#e6408e", "#afc84b", "#176cc0", "#332288", "#CC6677"),
    breaks = c(
      "transit",
      "full service restaurants", "religious organizations",
      "groceries and pharmacies", "colleges", "elementary and high schools"
    ),
    name = NULL
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b-%y") +
  ggtitle("Visits to Points of Interest (POIs)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ylim(c(-90, 90))
poi_plot

## declines during omicron
omicron_poi_plot <- ggplot() +
  geom_rect(aes(xmin = as.Date("2021-11-01"), xmax = as.Date("2022-02-01"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(mobility_metric %in% c(
        "transit", "elementary and high schools", "colleges",
        "full service restaurants", "religious organizations",
        "groceries and pharmacies"
      ) & date > as.Date("2021-09-01") & date < as.Date("2022-04-01")),
    aes(x = date, y = value, color = mobility_metric), lwd = 0.7, alpha = 0.8
  ) +
  theme_bw(base_size = 16) +
  xlab("Date") +
  ylab("Percent Change from Baseline") +
  scale_color_manual(
    values = c("#a45ccc", "#b88f06", "#4faa45", "#e6408e", "#afc84b", "#176cc0", "#332288", "#CC6677"),
    breaks = c(
      "transit",
      "full service restaurants", "religious organizations",
      "groceries and pharmacies", "colleges", "elementary and high schools"
    ),
    name = NULL
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 months", date_labels = "%b-%y", limits = c(as.Date("2021-09-01"), as.Date("2022-04-01"))) +
  ggtitle("Visits to Points of Interest (POIs)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ylim(c(-90, 90))
omicron_poi_plot

####################################################################################
## % Devices Staying Home and % Masking in Public: Figure 2C and Figure 16C
####################################################################################

combined_mob2 <-
  combined_mob %>%
  dplyr::select(date, fb_stay_put_custom, mask_wearing_final) %>%
  rename(start_date = date) %>%
  filter(!is.na(fb_stay_put_custom)) %>%
  mutate(fb_stay_put_custom = zoo::rollmean(fb_stay_put_custom, k = 7, align = "center", fill = NA)) %>%
  rename(date = start_date) %>%
  filter(date < as.Date("2022-07-01"))

coeff <- 2
stay_home_plot <- ggplot(combined_mob2) +
  geom_rect(aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue", data = data.frame(x = 0, y = 0)
  ) +
  geom_rect(aes(xmin = as.Date("2020-03-16"), xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "orange", data = data.frame(x = 0, y = 0)
  ) +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 1) +
  geom_line(aes(x = date, y = fb_stay_put_custom, color = "% Staying Home"), lwd = 1) +
  geom_line(aes(x = date, y = mask_wearing_final / coeff, color = "% Wearing Masks"), lwd = 1) +
  theme_bw(base_size = 16) +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b-%y", limits = c(as.Date("2018-11-18"), as.Date("2022-06-30"))) +
  ggtitle("Percentage Staying Home or Masking") +
  ylab("% Staying Home") +
  xlab("Date") +
  scale_y_continuous(
    # Features of the first axis
    name = "% Staying Home",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * coeff, name = "% Wearing Masks")
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title = element_blank()) +
  scale_color_jcolors(palette = "pal8")
stay_home_plot

## max devices staying at home during omicron wave = 47.6%
combined_mob2 %>%
  filter(date > as.Date("2021-09-01") & date < as.Date("2022-04-01")) %>%
  slice_max(fb_stay_put_custom)
# date       fb_stay_put_custom mask_wearing_final
# <date>                  <dbl>              <dbl>
# 1 2021-12-29               47.6               92.4

# increases during omicron
omicron_stay_home_plot <- ggplot() +
  geom_rect(aes(xmin = as.Date("2021-11-01"), xmax = as.Date("2022-02-01"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_line(aes(x = date, y = mask_wearing_final, color = "% Wearing Masks"), lwd = 1, data = combined_mob2 %>% filter(date > as.Date("2021-09-01") & date < as.Date("2022-04-01"))) +
  geom_line(aes(x = date, y = fb_stay_put_custom / 0.5, color = "% Staying Home"), lwd = 1, data = combined_mob2 %>% filter(date > as.Date("2021-09-01") & date < as.Date("2022-04-01"))) +
  theme_bw(base_size = 16) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 months", date_labels = "%b-%y", limits = c(as.Date("2021-09-01"), as.Date("2022-04-01"))) +
  ggtitle("Percentage Staying Home or Masking") +
  ylab("% Staying Home") +
  xlab("Date") +
  scale_y_continuous(

    # Features of the first axis
    name = "% Staying Home",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * 0.5, name = "% Wearing Masks")
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title = element_blank()) +
  scale_color_jcolors(palette = "pal8")
omicron_stay_home_plot

####################################################################################
## Figure 2
####################################################################################
combined <- plot_grid(flow_plot, poi_plot, stay_home_plot, nrow = 3, labels = "AUTO", rel_heights = c(0.75, 0.75, 0.75), align = "hv")
combined
save_plot(combined, filename = "1_Seattle_Mobility_Data/mobility_figures/fig_2_mobility_flow_and_poi_visits_time_series.png", base_width = 14, base_height = 14)
save_plot(combined, filename = "figures/fig_2_mobility_flow_and_poi_visits_time_series.png", base_width = 14, base_height = 14)
####################################################################################
## Figure S17
####################################################################################
omicron_combined <- plot_grid(omicron_flow_plot, omicron_poi_plot, omicron_stay_home_plot, nrow = 3, labels = "AUTO", rel_heights = c(0.75, 0.75, 0.75), align = "hv")
omicron_combined
save_plot(omicron_combined, filename = "1_Seattle_Mobility_Data/mobility_figures/fig_s17_omicron_period_mobility_flow_and_poi_visits_time_series.png", base_width = 10, base_height = 14)
save_plot(omicron_combined, filename = "figures/fig_s17_omicron_period_mobility_flow_and_poi_visits_time_series.png", base_width = 10, base_height = 14)
