## Note: code is for the raw SafeGraph visits, but Dewey doesn't permit sharing of the raw mobility data
## Requires outputs from 1_import_SafeGraph_Data.R but these can't be shared publicly
## Inputs and outputs are not publicly available, but you can run the code if you have access to the SafeGraph data
########################################################################################################
## 1. Combine SafeGraph Weekly Patterns dataset with Home and Visitor Panel data and census data
## 2. Adjust foot traffic indices by SG device panel size
## 3. Figure 3: Seattle Mobility Network Maps
## 4. Large scale movement mobility metrics
## 5. Foot traffic to different categories of POIs
########################################################################################################

library(SafeGraphR)
library(dplyr)
library(data.table)
library(readr)
library(tidyr)
library(padr)
library(igraph)
library(forecast)
library(censusapi)
library(tigris)
library(tidycensus)
Sys.setenv(PROJ_LIB = "")
library(sf) #  use older version of sf; see: https://github.com/r-spatial/sf/issues/2298#issuecomment-1867563910
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
# library(reshape) #reshape is needed but don't load whole package

## get API key here: https://api.census.gov/data/key_signup.html
# Add key to .Renviron
# Sys.setenv(CENSUS_KEY = "your key")
# Reload .Renviron
readRenviron("~/.Renviron")
# Check to see that the expected key is output in your R console
Sys.getenv("CENSUS_KEY")

dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Seattle_Flu_Study/Seattle_SG_Mobility/"
# ####################################################
# ## Seattle Mobility Network Maps: Figure 3 in main text
# ####################################################
# # downloaded from: https://data-seattlecitygis.opendata.arcgis.com/datasets/SeattleCityGIS::census-block-groups-2010/about
# ## Seattle CBGs
cbg_gps <- read.csv("1_Seattle_Mobility_Data/mobility_data/Seattle_Census_Block_Groups_2010.csv", colClasses = c("GEOID10" = "character"))
head(cbg_gps)
unique(nchar(cbg_gps$GEOID10))

# Map of King County, WA
# requires API key
king <- get_acs(
  state = "WA", county = "King", geography = "block group",
  variables = "B01003_001", geometry = TRUE, year = 2019
)

# filter to CBGs in Seattle proper
seattle_map <- king %>% filter(GEOID %in% cbg_gps$GEOID10)

weekly_visitors_cbg_within_seattle <- read_rds(paste0(dir, "SG_data/weekly_visitors_cbg_within_seattle_2018_2022.rds"))

range(weekly_visitors_cbg_within_seattle$scaled_visits_adj)

seattle_map$estimate <- as.numeric(seattle_map$estimate)

## check for duplicates
weekly_visitors_cbg_within_seattle %>%
  group_by(visitor_home_cbg, poi_cbg, start_date) %>%
  tally() %>%
  filter(n != 1)

range(weekly_visitors_cbg_within_seattle$start_date)

weekly_visitors_cbg_within_seattle %>%
  filter(start_date %in% c(
    as.Date("2019-02-11"),
    as.Date("2019-07-15"),
    as.Date("2019-03-18"),
    as.Date("2021-07-12"),
    as.Date("2021-11-24"),
    as.Date("2022-01-17")
  )) %>%
  pull(scaled_visits_county) %>%
  range() # 105.5147 23393.8629

# color on log scale
my_breaks <- reshape::round_any(exp(seq(log(100), log(24000), length = 5)), 5)

####################################################
## february 2019 snowstorm
####################################################
weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2019-02-11") & home_location == "non_resident") %>%
  summarize(net_movement = sum(scaled_visits_county)) #  984572

feb_2019_map <-
  ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-02-11") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-02-11") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6),
    axis.title.x = element_text(size = 5, hjust = 0.5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Feb 2019 Snowstorm") +
  xlab("Net Movement = 984,572") +
  ylab(NULL)
feb_2019_map

### histogram of degree values
network_df <- weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2019-02-11") & home_location == "non_resident") %>%
  dplyr::select(visitor_home_cbg, poi_cbg, scaled_visits_county) %>%
  as.data.frame() %>%
  rename(weight = scaled_visits_county)
range(unique(network_df$weight))

mygraph <- graph_from_data_frame(network_df, directed = F)
degree_dist <- degree_distribution(mygraph) %>% as.data.frame()
cbg_degree <- degree(mygraph, loops = F, normalized = F) %>% as.data.frame()
colnames(cbg_degree) <- "degree"

mean(cbg_degree$degree)
median(cbg_degree$degree)
range(cbg_degree$degree)
range(log10(cbg_degree$degree))

# log scale on x-axis
feb_deg_log <- ggpubr::gghistogram(cbg_degree,
  size = 0.2,
  font.label = list(size = 7),
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  theme_bw(base_size = 7) +
  theme(axis.title = element_text(size = 5)) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(1, 190))
feb_deg_log

feb_2019_map_w_deg_log <- plot_grid(feb_2019_map, feb_deg_log, nrow = 2, rel_heights = c(6, 2), labels = c("A", "B"), label_size = 7)
feb_2019_map_w_deg_log

feb_2019_map_w_deg_log <- plot_grid(feb_2019_map, feb_deg_log, nrow = 2, rel_heights = c(6, 2), labels = NA, label_size = 7)
feb_2019_map_w_deg_log
####################################################
## july 2019 baseline
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2019-07-22") & home_location == "non_resident") %>%
  summarize(net_movement = sum(scaled_visits_county)) # 1139033

jul_2019_map <- ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-07-22") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-07-22") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6),
    axis.title.x = element_text(size = 5, hjust = 0.5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("July 2019") +
  xlab("Net Movement =  1,139,033") +
  ylab(NULL)

network_df <- weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2019-07-22") & home_location == "non_resident") %>%
  dplyr::select(visitor_home_cbg, poi_cbg, scaled_visits_county) %>%
  as.data.frame() %>%
  rename(weight = scaled_visits_county)

mygraph <- graph_from_data_frame(network_df, directed = F)
degree_dist <- degree_distribution(mygraph) %>% as.data.frame()
cbg_degree <- degree(mygraph, loops = F, normalized = F) %>% as.data.frame()
colnames(cbg_degree) <- "degree"
mean(cbg_degree$degree)
median(cbg_degree$degree)
range(cbg_degree$degree)

# log scale on x-axis
jul_2019_log <- ggpubr::gghistogram(cbg_degree,
  size = 0.2,
  font.label = list(size = 7),
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  theme_bw(base_size = 7) +
  theme(axis.title = element_text(size = 5)) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(1, 190))
jul_2019_log


jul_2019_map_w_deg_log <- plot_grid(jul_2019_map, jul_2019_log, nrow = 2, rel_heights = c(6, 2))
jul_2019_map_w_deg_log

####################################################
## march 2020 lockdown
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "non_resident" & start_date == as.Date("2020-03-16")) %>%
  summarize(net_movement = sum(scaled_visits_county)) # 495167

march_2020_map <-
  ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2020-03-16") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2020-03-16") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6),
    axis.title.x = element_text(size = 5, hjust = 0.5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("March 2020 Lockdown") +
  xlab("Net Movement = 495,167") +
  ylab(NULL)
march_2020_map

network_df <- weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2020-03-16") & home_location == "non_resident") %>%
  dplyr::select(visitor_home_cbg, poi_cbg, scaled_visits_county) %>%
  as.data.frame() %>%
  rename(weight = scaled_visits_county)

mygraph <- graph_from_data_frame(network_df, directed = F)
degree_dist <- degree_distribution(mygraph) %>% as.data.frame()
cbg_degree <- degree(mygraph, loops = F, normalized = F) %>% as.data.frame()
colnames(cbg_degree) <- "degree"
mean(cbg_degree$degree)
median(cbg_degree$degree)
range(cbg_degree$degree)

## degree histogram
mar_2020_log <- ggpubr::gghistogram(cbg_degree,
  size = 0.2,
  font.label = list(size = 7),
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  theme_bw(base_size = 7) +
  theme(axis.title = element_text(size = 5)) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(1, 190))
mar_2020_log

mar_2020_map_w_deg_log <- plot_grid(march_2020_map, mar_2020_log, nrow = 2, rel_heights = c(6, 2))
mar_2020_map_w_deg_log

####################################################
## Delta wave: July 2021
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2021-07-05") & home_location == "non_resident") %>%
  summarize(net_movement = sum(scaled_visits_county)) # 721159

july_2021_map <-
  ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") &
        home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visits", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6),
    axis.title.x = element_text(size = 5, hjust = 0.5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Delta Wave (July 2021)") +
  xlab("Net Movement = 721,159") +
  ylab(NULL)
july_2021_map

network_df <- weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2021-07-05") & home_location == "non_resident") %>%
  dplyr::select(visitor_home_cbg, poi_cbg, scaled_visits_county) %>%
  as.data.frame() %>%
  rename(weight = scaled_visits_county)

mygraph <- graph_from_data_frame(network_df, directed = F)
degree_dist <- degree_distribution(mygraph) %>% as.data.frame()
cbg_degree <- degree(mygraph, loops = F, normalized = F) %>% as.data.frame()
colnames(cbg_degree) <- "degree"
mean(cbg_degree$degree)
median(cbg_degree$degree)
range(cbg_degree$degree)

## degree histogram
jul_2021_log <- ggpubr::gghistogram(cbg_degree,
  size = 0.2,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  theme_bw(base_size = 7) +
  theme(axis.title = element_text(size = 5)) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(1, 190))
jul_2021_log

jul_2021_map_w_deg_log <- plot_grid(july_2021_map, jul_2021_log, nrow = 2, rel_heights = c(6, 2))
jul_2021_map_w_deg_log

####################################################
## Omicron wave Jan 2022 (avoid potential effects of Christmas/NYE on movement)
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "non_resident" &
    start_date == as.Date("2022-01-17")) %>%
  summarize(net_movement = sum(scaled_visits_county)) # 703591

weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "non_resident" &
    start_date == as.Date("2022-01-17")) %>%
  pull(scaled_visits_county) %>%
  range()

jan_2022_map <-
  ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2022-01-17") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2022-01-17") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6),
    axis.title.x = element_text(size = 5, hjust = 0.5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Omicron Wave (January 2022)") +
  xlab("Net Movement = 703,591") +
  ylab(NULL)

network_df <- weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2022-01-17") & home_location == "non_resident") %>%
  dplyr::select(visitor_home_cbg, poi_cbg, scaled_visits_county) %>%
  as.data.frame() %>%
  rename(weight = scaled_visits_county)

mygraph <- graph_from_data_frame(network_df, directed = F)
degree_dist <- degree_distribution(mygraph) %>% as.data.frame()
cbg_degree <- degree(mygraph, loops = F, normalized = F) %>% as.data.frame()
colnames(cbg_degree) <- "degree"
mean(cbg_degree$degree)
median(cbg_degree$degree)
range(cbg_degree$degree)

## degree histogram
jan_2022_log <- ggpubr::gghistogram(cbg_degree,
  size = 0.2,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  theme_bw(base_size = 7) +
  theme(axis.title = element_text(size = 5)) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(1, 190))
jan_2022_log

jan_2022_map_w_deg_log <- plot_grid(jan_2022_map, jan_2022_log, nrow = 2, rel_heights = c(6, 2))
jan_2022_map_w_deg_log

# remake to get legend for figure with all maps
july_2021_map <-
  ggplot(seattle_map) +
  geom_sf(lwd = 0.1, alpha = 0.2) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") & home_location == "non_resident" & scaled_visits_county < 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county,
      alpha = scaled_visits_county
    ), lwd = 0.2
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") &
        home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.4
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of\nVisits", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.justification = "center",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Delta Wave (July 2021)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Net Movement = 721,159") +
  ylab(NULL)
july_2021_map

## maps with degree histogram on log scale = Figure 3
hist_combined <- plot_grid(
  feb_2019_map_w_deg_log + theme(legend.position = "none"),
  jul_2019_map_w_deg_log + theme(legend.position = "none"),
  mar_2020_map_w_deg_log + theme(legend.position = "none"),
  jul_2021_map_w_deg_log + theme(legend.position = "none"),
  jan_2022_map_w_deg_log + theme(legend.position = "none"),
  nrow = 1
)
hist_combined
leg <- get_legend(july_2021_map + theme(legend.direction = "vertical"))

map_combined_and_hist <- plot_grid(hist_combined, leg, nrow = 1, rel_widths = c(1, 0.1))
map_combined_and_hist

## Figure 3
save_plot(map_combined_and_hist,
  filename = "figures/fig_3_within_seattle_movement_up_to_Jan_2022_with_log_trans_deg_histogram.pdf",
  units = "mm", base_width = 180, base_height = 90, dpi = 300
)

#presentation
save_plot(map_combined_and_hist,
          filename = "figures/presentation_within_seattle_movement_up_to_Jan_2022_with_log_trans_deg_histogram.png",
          units = "mm", base_width = 180, base_height = 90, dpi = 300
)
