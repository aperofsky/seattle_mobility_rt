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
  geom_hline(aes(yintercept = 0), lty = "dashed",alpha=0.5) +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(date >= as.Date("2018-11-18") & date < as.Date("2022-07-01")) %>%
      filter(mobility_metric %in% c(
        "influx out-of-state visitors",
        "influx other WA counties",
        "between-neighborhood movement",
        "within-neighborhood movement"
      )),
    aes(x = date, y = value, color = mobility_metric), lty = "solid", lwd = 0.7, alpha = 0.8
  ) +
  theme_bw(base_size = 7) +
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
  geom_hline(aes(yintercept = 0), lty = "dashed",alpha=0.5) +
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_line(
    data = combined_long %>% filter(!is.na(value)) %>%
      filter(date >= as.Date("2018-11-18") & date < as.Date("2022-07-01")) %>%
      filter(mobility_metric %in% c(
        "transit", "elementary and high schools", "colleges",
        "full service restaurants", "religious organizations",
        "groceries and pharmacies"
      )),
    aes(x = date, y = value, color = mobility_metric), lwd = 0.5, alpha = 0.8
  ) +
  theme_bw(base_size = 7) +
  xlab("Date") +
  ylab("Percent Change from Baseline") +
  scale_color_manual(
    values = c("#a45ccc", "#b88f06", "#4faa45", "#e6408e", "#afc84b", "#176cc0", "#332288", "#CC6677"),
    breaks = c(
      "transit",
      "full service restaurants", "religious organizations",
      "groceries and pharmacies", "colleges", "elementary and high schools"
    ),
    labels=c(
      "Transit",
      "Full service restaurants", "Religious organizations",
      "Groceries and pharmacies", "Colleges", "Elementary and high schools"
    ),
    name = NULL
  ) +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b-%y") +
  ggtitle("Visits to Points of Interest (POIs)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 1))+
  ylim(c(-90, 90))
poi_plot

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
  geom_vline(xintercept = as.Date("2020-02-28"), lty = "dashed", color = "darkgreen", lwd = 0.5) +
  geom_line(aes(x = date, y = fb_stay_put_custom, color = "% Staying Home"), lwd = 0.7) +
  geom_line(aes(x = date, y = mask_wearing_final / coeff, color = "% Wearing Masks"), lwd = 0.7) +
  theme_bw(base_size = 7) +
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

####################################################################################
## Figure 2
####################################################################################
combined <- plot_grid(flow_plot, poi_plot, stay_home_plot, nrow = 3, labels = "AUTO",
                      label_size = 7, rel_heights = c(0.75, 0.75, 0.75), align = "hv")
combined
# save_plot(combined, filename = "1_Seattle_Mobility_Data/mobility_figures/fig_2_mobility_flow_and_poi_visits_time_series.png", base_width = 14, base_height = 14)
save_plot(combined, filename = "figures/fig_2_mobility_flow_and_poi_visits_time_series.pdf",  units="mm",base_width = 180, base_height = 180, dpi=300)

