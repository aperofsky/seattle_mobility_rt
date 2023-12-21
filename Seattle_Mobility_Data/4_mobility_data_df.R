
## Note: code is for the raw SafeGraph visits, but Dewey doesn't permit sharing of the raw mobility data
## requires outputs from 1_import_SafeGraph_Data.R and 2_process_Safegraph_data_and_city_maps.R but these can't be shared publicly

##########################################################################
## Combine SG data with other population behavior datasets
##########################################################################

library(dplyr)
library(tidyr)
library(forecast)
library(cdcfluview)
library(data.table)
library(readr)
library(timeDate)
library(lubridate)

# dir <- "Seattle_SG_Mobility/SG_data/"

## Safegraph large scale pop movements
movement <- read_rds(paste0(dir, "within_and_between_cbg_movement_indicators_up_to_2022_09_26_red.rds"))
range(movement$start_date) # "2018-06-04" "2022-09-26"
head(movement) # raw counts adjusted by panel size
names(movement)

## Safegraph visits to categories of POIs
seattle_daily_visits_naics_all_years <- read_rds(paste0(dir, "seattle_daily_visits_naics.rds"))
head(seattle_daily_visits_naics_all_years) # visits to categories of pois adjusted by panel size
range(seattle_daily_visits_naics_all_years$date) # "2018-06-04" "2022-10-02"
unique(seattle_daily_visits_naics_all_years$industry) # 283 categories
names(seattle_daily_visits_naics_all_years)

## External data sources: King Co. masking survey data, Oxford stringency stringency, FB stay put
ext_behavior_daily_data <- read_csv("Seattle_Mobility_Data/mobility_data/covidcast_fb_osi_combined_dataset.csv")
head(ext_behavior_daily_data)

##########################################################################
## Filter SG POI categories
##########################################################################

## keep SG categories with most visits
keep <- seattle_daily_visits_naics_all_years %>%
  mutate(year = lubridate::year(date)) %>%
  filter(!(industry %in% c(
    "Nature Parks and Other Similar Institutions", # parks are "sinks" in Seattle (inflated pings)
    "Other Airport Operations", "Assisted Living Facilities for the Elderly",
    "Pet and Pet Supplies Stores", "Florists",
    "Hardware Stores", "Sporting Goods Stores", "Used Car Dealers",
    "Gift, Novelty, and Souvenir Stores", "Kidney Dialysis Centers",
    "Used Merchandise Stores", "Golf Courses and Country Clubs",
    "Automotive Parts and Accessories Stores", "Historical Sites",
    "New Car Dealers", "Book Stores", "Tobacco Stores",
    "Home Health Care Services", "Florists",
    "Hobby, Toy, and Game Stores", "Sports and Recreation Instruction",
    "Casinos (except Casino Hotels)", "Department Stores",
    "Automotive Repair and Maintenance", "Furniture Stores",
    "All Other General Merchandise Stores"
  ))) %>%
  group_by(industry) %>%
  mutate(visits_fixed = tsclean(scaled_visits_state)) %>%
  group_by(year, industry) %>%
  summarize(total_visits = sum(visits_fixed, na.rm = T)) %>%
  arrange(year, -total_visits) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  distinct(industry) %>%
  pull(industry)
keep

## add in additional categories
keep2 <- c(
  keep, "Breweries", "Offices of Dentists", "Transit and Ground Passenger Transportation",
  "All Other Transit and Ground Passenger Transportation",
  "Colleges, Universities, and Professional Schools", "General Medical and Surgical Hospitals",
  "Sports Teams and Clubs", "Zoos and Botanical Gardens"
) %>% unique()
keep2

## aggregate by date and POI category
all_bus_scaled <- seattle_daily_visits_naics_all_years %>%
  group_by(date, start_date, industry) %>%
  summarize(visits_by_day = sum(scaled_visits_state, na.rm = T)) %>%
  ungroup()

## filter to select categories
naics_filt <- all_bus_scaled %>%
  filter(industry %in% keep2) %>%
  ungroup() %>%
  distinct()

unique(naics_filt$industry) # 30 categories

head(naics_filt)
names(naics_filt)

## clean time series
naics_filt <-
  naics_filt %>%
  group_by(industry) %>%
  mutate(visits_fixed = tsclean(visits_by_day)) %>%
  ungroup()

naics_filt$industry <- as.factor(naics_filt$industry)
levels(naics_filt$industry)

## shorten names of POI categories
levels(naics_filt$industry) <- c(
  "amusement_and_recreation",
  "other_transit",
  "amusement_and_theme_parks",
  "breweries",
  "child_day_care",
  "colleges",
  "convenience_stores",
  "bars",
  "elementary_and_secondary_schools",
  "fitness",
  "full_service_restaurants",
  "gas_stations",
  "hospitals",
  "hotels",
  "non_residential_buildings",
  "libraries",
  "limited_service_restaurants",
  "dentist_offices",
  "optometrist_offices",
  "physician_offices",
  "pharmacies",
  "performing_arts_or_sports_events",
  "religious_orgs",
  "snack_and_beverage_bars",
  "sports_teams_and_clubs",
  "supermarkets",
  "transit",
  "warehouse_clubs",
  "zoos"
)

levels(naics_filt$industry)
unique(naics_filt$industry)

## make POI categories columns
naics_wide <- naics_filt %>%
  dplyr::select(-visits_by_day, -start_date) %>%
  pivot_wider(names_from = industry, values_from = visits_fixed) %>%
  filter(date >= as.Date("2018-06-01"))

##########################################################################
## Join SG datasets
##########################################################################
## join daily data on POI visits with weekly data on large scale movements
combined_SG_mob <- left_join(naics_wide, movement, by = c("date" = "start_date")) %>%
  arrange(date)%>%
  fill(c(within_neighborhood_movement:out_of_state_movement), .direction = "down") # fill in movement data (weekly); POI visits are daily
head(combined_SG_mob)
combined_SG_mob %>% dplyr::select(date,within_neighborhood_movement:out_of_state_movement)

## combine some categories
combined_SG_mob$schools_and_daycare <- combined_SG_mob$colleges + combined_SG_mob$elementary_and_secondary_schools + combined_SG_mob$child_day_care
combined_SG_mob$all_restaurants <- combined_SG_mob$full_service_restaurants + combined_SG_mob$limited_service_restaurants + combined_SG_mob$snack_and_beverage_bars
combined_SG_mob$groceries_and_pharmacies <- combined_SG_mob$supermarkets + combined_SG_mob$pharmacies + combined_SG_mob$warehouse_clubs

combined_SG_mob$all_amusement_and_recreation <- combined_SG_mob$amusement_and_recreation +
  combined_SG_mob$amusement_and_theme_parks +
  combined_SG_mob$performing_arts_or_sports_events +
  combined_SG_mob$zoos +
  combined_SG_mob$sports_teams_and_clubs

combined_SG_mob$doctor_and_dentist_offices <- combined_SG_mob$dentist_offices + combined_SG_mob$physician_offices + combined_SG_mob$optometrist_offices

combined_SG_mob$transit <- combined_SG_mob$transit + combined_SG_mob$other_transit

combined_SG_mob$bars_and_breweries <- combined_SG_mob$bars + combined_SG_mob$breweries

##########################################################################
## % change from baseline for SG metrics
##########################################################################

combined_SG_mob$wday <- lubridate::wday(combined_SG_mob$date, label = T)

HolidayTable <- rbindlist(lapply(2018:2022, function(y) {
  data.frame(Holiday = as.Date(c(USNewYearsDay(y), USMemorialDay(y), USIndependenceDay(y), USLaborDay(y), USThanksgivingDay(y), USChristmasDay(y))))
}))

combined_SG_mob <- combined_SG_mob %>% mutate(holiday = if_else(date %in% HolidayTable$Holiday, "yes", "no"))

## mobility metrics
cols2 <- names(combined_SG_mob)[!(names(combined_SG_mob) %in% c("date", "wday", "holiday", "other_transit"))]

## feb 2019 snowstorm dates
snow_dates <- seq.Date(from = as.Date("2019-02-03"), to = as.Date("2019-02-15"), by = "day")

scale_to_date_custom <- function(data, date1, date2, columns = c(
                                   "out_of_state_movement_king", "within_state_movement_king",
                                   "within_city_movement_king", "within_neighborhood_movement_king"
                                 )) {
  mean <- data %>%
    filter(!(start_date %in% snow_dates)) %>%
    filter(holiday == "no") %>%
    filter(start_date >= as.Date(date1), start_date <= as.Date(date2)) %>%
    summarize_at(.vars = columns, mean) %>%
    pivot_longer(cols = all_of(columns), names_to = "metric", values_to = "mean")
  # summarize_at(all_of(columns), .funs = "mean") %>% #all_of() outside of select is now deprecated
  # pivot_longer(all_of(columns), names_to = "metric", values_to = "mean")

  # ((after value – before value) / before value) * 100 = % change
  df_long <- data %>%
    pivot_longer(cols = all_of(columns), names_to = "metric", values_to = "value")

  df_long <- left_join(df_long, mean, by = "metric")

  df2 <- df_long %>%
    mutate(per_change = ((value - mean) / mean) * 100) %>%
    dplyr::select(start_date, metric, per_change) %>%
    pivot_wider(names_from = metric, values_from = per_change)
  return(df2)
}


## scale mobility data relative to 2019 (% change in baseline)
combined_SG_mob_rel <- scale_to_date_custom(
  data = combined_SG_mob %>%
    rename(start_date = date) %>%
    dplyr::select(-contains(c("fips", "home_panel", "pop", "SG_visits", "multiplier", "scaled"))),
  date1 = "2019-01-01", date2 = "2019-12-31", columns = cols2
)
head(combined_SG_mob_rel)

##########################################################################
## combine SG data with external datasets
##########################################################################

combined_mob <- full_join(combined_SG_mob_rel, ext_behavior_daily_data, by = c("start_date" = "date")) %>% rename(date = start_date)
combined_mob <- bind_cols(combined_mob, as_tibble(mmwr_week(combined_mob$date)[, 1:2])) %>% as.data.table()
combined_mob$epi_date <- mmwr_week_to_date(combined_mob$mmwr_year, combined_mob$mmwr_week)
combined_mob <- combined_mob %>% as_tibble()

# further filter POIs to select categories
combined_mob2 <- combined_mob %>% dplyr::select(
  date,epi_date,
  contains("movement"),
  full_service_restaurants,
  groceries_and_pharmacies,
  full_service_restaurants,
  performing_arts_or_sports_events,
  transit,
  religious_orgs,
  child_day_care,
  elementary_and_secondary_schools,
  colleges,
  fb_stay_put_custom,
  oxford_stringency_index,
  mask_wearing_final
)
write_rds(combined_mob2, file = "Seattle_Mobility_Data/mobility_data/mobility_metrics_for_epidemia.rds")
write_rds(combined_mob2, file = "Epidemia_Models/mobility_metrics_for_epidemia.rds")
write_rds(combined_mob2, file = "Block_Bootstrap_Mobility_vs_Rt//mobility_metrics_for_epidemia.rds")

