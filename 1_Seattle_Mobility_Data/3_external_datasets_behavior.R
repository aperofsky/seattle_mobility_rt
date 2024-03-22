####################################################################################
# Import and combine external data sources on % staying home, % masking, NPI stringency
####################################################################################
# Note: these outputs are already available in the Seattle_Mobility_Data/mobility_data folder
# Safegraph % staying home (via Delphi covidcast)
# Facebook % masking in public (via Delphi covidcast)
# Outbreaks Near You (ONY) % masking in public (via Rader et al 2021 github repo, obtained from Leech 2020)
# Oxford stringency index
# Facebook % staying put (movement range maps from humanitarian data exchange)

library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(timetk)
library(lubridate)
library(covidcast)
library(cdcfluview)
library(ggplot2)
library(cowplot)
library(jcolors)
####################################################################################
## COVIDCAST DATA: SG % staying home and FB % masking in public
####################################################################################
## note: you will need an API key from COVIDcast if you want to pull a lot of data
## options(covidcast.auth = "your_key")
seattle_counties <- "53033" # fips code for King County, WA

##########################################
## Safegraph % staying home
SG <- covidcast_signals(
  data_source = "safegraph",
  signal = c("completely_home_prop_7dav", "median_home_dwell_time_7dav"),
  geo_type = "county",
  geo_values = seattle_counties
)

SG_agg <- aggregate_signals(SG)
head(SG_agg)
range(SG_agg$time_value) # "2019-01-01" "2021-04-16"

names(SG_agg) <- c("fips", "date", "completely_home_prop", "median_home_dwell_time")

ggplot(SG_agg) +
  geom_line(aes(x = date, y = completely_home_prop))
today <- Sys.Date()
today <- gsub("-", "_", today)
write_csv(SG_agg, file = paste0("1_Seattle_Mobility_Data/mobility_data/safegraph_percent_staying_home_king_co_wa_", today, ".csv"))

##########################################
## pull FB survey data on proportion of individuals masking in King County

mask_wearing <- covidcast_signals(
  data_source = "fb-survey",
  signal = c("smoothed_wwearing_mask_7d", "smoothed_wwearing_mask"),
  geo_type = "county",
  geo_values = seattle_counties
)

mask_wearing_agg <- aggregate_signals(mask_wearing)
head(mask_wearing_agg)
names(mask_wearing_agg) <- c("fips", "date", "mask_wearing_new", "mask_wearing_old")

mask_wearing_agg <-
  mask_wearing_agg %>%
  mutate(mask_wearing_final = ifelse(date > as.Date("2021-02-10"), mask_wearing_new, mask_wearing_old))

ggplot(mask_wearing_agg) +
  geom_line(aes(x = date, y = mask_wearing_new)) +
  geom_line(aes(x = date, y = mask_wearing_old)) +
  geom_line(aes(x = date, y = mask_wearing_final), color = "blue") +
  facet_wrap(~fips)

range(mask_wearing_agg$date) # "2020-09-08" "2022-06-25"

mask_wearing_agg <- mask_wearing_agg %>% arrange(date)
mask_wearing_agg$date <- as.Date(mask_wearing_agg$date)
today <- Sys.Date()
today <- gsub("-", "_", today)
write_csv(mask_wearing_agg, file = paste0("1_Seattle_Mobility_Data/mobility_data/covidcast_masking_king_co_wa_", today, ".csv"))

## supplement FB masking data with Outbreaks Near You daily masking data for WA state because it starts in June 2020
# wa_mask <- read_csv("https://raw.githubusercontent.com/g-leech/masks_v_mandates/main/data/raw/rader/rader_us_wearing_aggregated_mean_shop_and_work.csv") %>%
#   filter(label == "Washington") %>%
#   dplyr::select(-label)
# head(wa_mask)
# write_csv(wa_mask, file = "Seattle_Mobility_Data/mobility_data/ONY_wa_state_mask_data.csv")
wa_mask_rader <- read_csv("1_Seattle_Mobility_Data/mobility_data/ONY_wa_state_mask_data.csv")

head(mask_wearing_agg)

## combine masking data
skc_and_wa_masking <- full_join(mask_wearing_agg, wa_mask_rader, by = c("date" = "response_date")) %>%
  mutate(percent_mc = (100 * percent_mc) + 5) %>% # trends are the same but wa state was slightly lower than skc
  arrange(date) %>%
  mutate(mask_wearing_final = if_else(date > as.Date("2020-09-08"), mask_wearing_final, percent_mc)) %>%
  dplyr::select(-fips, -mask_wearing_new, -mask_wearing_old, -percent_mc)
range(skc_and_wa_masking$date) # "2020-06-02" "2022-06-25"
write_csv(skc_and_wa_masking, file = "1_Seattle_Mobility_Data/mobility_data/skc_masking_june_2020_to_june_2022.csv")

## covidcast indicators
covidcast_ind <- full_join(SG_agg, skc_and_wa_masking, by = c("date"))
covidcast_ind <- covidcast_ind %>% arrange(date)
head(covidcast_ind)
range(covidcast_ind$date)

###########################################
## Oxford Stringency Indicator
###########################################
ox <- url("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/United%20States/OxCGRT_USA_latest.csv")
ox_df <- read.csv(ox)
head(ox_df)
unique(ox_df$RegionName)
unique(ox_df$Jurisdiction)

wa_ox_df <- ox_df %>% filter(RegionName == "Washington")
wa_ox_df$Date <- lubridate::as_date(as.character(wa_ox_df$Date), format = "%Y%m%d", tz = NULL)
head(wa_ox_df)
tail(wa_ox_df)
range(wa_ox_df$Date)
names(wa_ox_df)

ggplot(wa_ox_df) +
  geom_line(aes(x = Date, y = ContainmentHealthIndex_Average_ForDisplay)) +
  geom_line(aes(x = Date, y = EconomicSupportIndex_ForDisplay), color = "blue") +
  geom_line(aes(x = Date, y = GovernmentResponseIndex_Average_ForDisplay), color = "green") +
  geom_line(aes(x = Date, y = StringencyIndex_Average_ForDisplay), color = "red")

wa_ox_df_lim <- wa_ox_df %>%
  dplyr::select(Date, GovernmentResponseIndex_Average_ForDisplay, StringencyIndex_Average_ForDisplay) %>%
  rename(
    oxford_stringency_index = StringencyIndex_Average_ForDisplay,
    oxford_gov_response_index = GovernmentResponseIndex_Average_ForDisplay,
    date = Date
  )
head(wa_ox_df_lim)
today <- Sys.Date()
today <- gsub("-", "_", today)
write_csv(wa_ox_df_lim, file = paste0("1_Seattle_Mobility_Data/mobility_data/oxford_stringency_index_washington_state_", today, ".csv"))

###########################################
## Facebook mobility data
## Source: https://data.humdata.org/dataset/movement-range-maps; discontinued May 22 2022
## Note the full dataset is not included in this GH repo because the files are too large
## If you want to work with the original dataset, download data to "Seattle_Mobility_Data/mobility_data/" folder and uncomment first section
###########################################
## original fb dataset including all locations (large files)
# fb_movement1 = fread("1_Seattle_Mobility_Data/mobility_data/movement-range-data-2020-03-01-2020-12-31/movement-range-data-2020-03-01--2020-12-31.txt",header = T,fill = T)
# fb_movement2 = fread("1_Seattle_Mobility_Data/mobility_data/movement-range-data-2022-05-22/movement-range-2022-05-22.txt",header = T,fill = T)
# fb_movement = bind_rows(fb_movement1,fb_movement2)

# all_day_bing_tiles_visited_relative_change: Positive or negative change in movement relative to baseline
# all_day_ratio_single_tile_users: Positive proportion of users staying put within a single location

### filter to King County
# skc_fb <- fb_movement %>%
#   filter(country == "USA" & polygon_source == "FIPS" & polygon_id == "53033") %>%
#   as_tibble()
# write_csv(skc_fb, file = "1_Seattle_Mobility_Data/mobility_data/skc_facebook_stay_at_home_raw.csv")

## easier to start with fb data filtered to king county
skc_fb <- read_csv("1_Seattle_Mobility_Data/mobility_data/skc_facebook_stay_at_home_raw.csv", col_types = c("D", "c", "c", "c", "c", "d", "d", "c", "c"))
skc_fb$all_day_bing_tiles_visited_relative_change <- as.numeric(skc_fb$all_day_bing_tiles_visited_relative_change)
skc_fb$all_day_ratio_single_tile_users <- as.numeric(skc_fb$all_day_ratio_single_tile_users)
skc_fb$ds <- as.Date(skc_fb$ds)
skc_fb <- skc_fb[complete.cases(skc_fb), ]

ggplot(skc_fb) +
  geom_line(aes(x = ds, y = all_day_bing_tiles_visited_relative_change))

ggplot(skc_fb) +
  geom_line(aes(x = ds, y = all_day_ratio_single_tile_users))

## check for duplicates
skc_fb %>%
  group_by(ds) %>%
  tally() %>%
  filter(n != 1)

skc_fb <-
  skc_fb %>%
  rename(
    fb_movement = all_day_bing_tiles_visited_relative_change,
    fb_stay_put = all_day_ratio_single_tile_users,
    date = ds
  ) %>%
  dplyr::select(date, polygon_id, polygon_name, fb_movement, fb_stay_put)

skc_fb <-
  skc_fb %>%
  mutate(fb_movement = fb_movement * 100, fb_stay_put = fb_stay_put * 100)
head(skc_fb)

today <- Sys.Date()
today <- gsub("-", "_", today)
write_csv(skc_fb, file = paste0("1_Seattle_Mobility_Data/mobility_data/fb_skc_mobility_as_of_", today, ".csv"))


rolling_avg_7 = slidify(.f = mean, .period = 7, .align = "center", .partial = F)

external_dataset <- plyr::join_all(
  list(
    data.frame(covidcast_ind %>% dplyr::select(-median_home_dwell_time)),
    data.frame(wa_ox_df_lim %>% dplyr::select(-oxford_gov_response_index)),
    data.frame(skc_fb %>% dplyr::select(date, fb_stay_put))
  ),
  type = "full",
  by = c("date")
) %>%
  mutate(fb_stay_put_custom = fb_stay_put + 14) %>%
  mutate(fb_stay_put_custom = if_else(date < as.Date("2020-03-01"), 100 * completely_home_prop, fb_stay_put_custom)) %>%
  mutate(fb_stay_put_custom = rolling_avg_7(fb_stay_put_custom))%>% #fb data is much noiser than SG data
  dplyr::select(date,completely_home_prop,fb_stay_put, fb_stay_put_custom, mask_wearing_final, oxford_stringency_index)

external_dataset %>% filter(!is.na(fb_stay_put)) %>% last()

a = ggplot(external_dataset %>% filter(date<as.Date("2022-05-23")))+
  geom_line(aes(x=date,y=fb_stay_put/100,color="Meta % Stay Put"),linewidth=1,alpha=0.8)+
  geom_line(aes(x=date,y=completely_home_prop,color="SafeGraph % Completely Staying at Home"),linewidth=1,alpha=0.8)+
  ylab("Percent Staying Home")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits=c(0.15,0.6))+
  theme_bw(base_size = 18)+
  scale_color_manual(name=NULL,values=as.vector(jcolors("pal8")),
                     breaks=c("SafeGraph % Completely Staying at Home", "Meta % Stay Put"))+
  scale_x_date("Date",expand=c(0,0),date_breaks = "6 months",date_labels = "%b-%Y")+
  theme(legend.position = "bottom")
a

b = ggplot(external_dataset %>% filter(date<as.Date("2022-05-20")))+
  geom_line(aes(x=date,y=(fb_stay_put + 14)/100,color="Adjusted Meta % Stay Put"),linewidth=1,alpha=0.8)+
  geom_line(aes(x=date,y=completely_home_prop,color="SafeGraph % Completely Staying at Home"),linewidth=1)+
  ylab("Percent Staying Home")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits=c(0.15,0.6))+
  theme_bw(base_size = 18)+
  scale_color_manual(name=NULL,values=as.vector(jcolors("pal8")),
                     breaks=c("SafeGraph % Completely Staying at Home", "Adjusted Meta % Stay Put"))+
  scale_x_date("Date",expand=c(0,0),date_breaks = "6 months",date_labels = "%b-%Y")+
  theme(legend.position = "bottom")
b

c = ggplot(external_dataset %>% filter(date<as.Date("2022-05-20")))+
  geom_line(aes(x=date,y=fb_stay_put_custom/100,color="SafeGraph + Meta Custom % Staying at Home"),linewidth=1)+
  ylab("Percent Staying Home")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits=c(0.15,0.6))+
  theme_bw(base_size = 18)+
  scale_color_manual(name=NULL,values=as.vector(jcolors("pal2")),
                     breaks=c("SafeGraph + Meta Custom % Staying at Home"))+
  scale_x_date("Date",expand=c(0,0),date_breaks = "6 months",date_labels = "%b-%Y")+
  theme(legend.position = "bottom")
c

combined = plot_grid(a,b,c,nrow=3,labels="AUTO")
combined
save_plot(plot=combined,file="figures/fig_s24_custom_stay_home_metric.png",base_width = 16,base_height = 16)

external_dataset <- external_dataset %>% dplyr::select(date, fb_stay_put_custom, mask_wearing_final, oxford_stringency_index)
head(external_dataset)

write_csv(external_dataset,file="1_Seattle_Mobility_Data/mobility_data/covidcast_fb_osi_combined_dataset.csv")