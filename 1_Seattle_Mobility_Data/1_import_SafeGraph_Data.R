## Note: code is for the raw SafeGraph visits, but Dewey doesn't permit sharing of the raw mobility data.
## Inputs and outputs are not publicly available, but you can run the code if you have access to the raw SafeGraph data.
####################################################
## 1. Import and organize/clean SafeGraph data:
### - Home panel
### - Visitor panel
### - Normalization stats
### - Weekly patterns
## 2. Create mobility indices:
### - Daily visits to various categories of POIs
### - Weekly visitors by visitor home CBG
####################################################
#### Note: SafeGraph data are for WA state only

library(SafeGraphR)
library(dplyr)
library(data.table)
library(ggplot2)
library(readr)
library(fasttime)
library(lubridate)
library(dtplyr)
library(tigris)

data("fips_codes") # tigris package

####################################################
## Home Summary Data
####################################################
panel_dir <- "/Volumes/My Passport for Mac/SafeGraph_Dewey/home_panel_summary/" # path to directory of your home panel data

# Make sure that census_block_group is read as a character
# note colClasses is a data.table::fread() argument
all_panel <- read_many_csvs(
  dir = panel_dir, recursive = T,
  colClasses = c(census_block_group = "character")
)
all_panel[, census_block_group := as.character(as.numeric(census_block_group))]
head(all_panel)
range(all_panel$date_range_start)
range(all_panel$date_range_end)
all_panel %>% filter(iso_country_code == "US")

all_panel <- all_panel[complete.cases(all_panel), ]

## check for duplicates
all_panel[duplicated(all_panel[, c("date_range_start", "census_block_group", "region", "iso_country_code")]), ] %>%
  arrange(region, census_block_group, date_range_start)

all_panel <- unique(all_panel)

all_panel[, start_date := fasttime::fastPOSIXct(substr(date_range_start, 0, 10), tz = "UTC")]

## filter out Canada
all_panel_US <- all_panel %>% filter(iso_country_code == "US")
all_panel_US <- all_panel_US %>% as.data.table()
all_panel_US[, date_range_start := NULL]
all_panel_US[, date_range_end := NULL]
all_panel_US[, c("state_fips", "county_fips") := fips_from_cbg(census_block_group)]
all_panel_US$week_date <- lubridate::floor_date(all_panel_US$start_date, unit = "week")
all_panel_US[, start_date := as.Date(week_date) + 1]
all_panel_US <- all_panel_US %>% filter(!(region %in% c("gu", "mp")))
all_panel_US <- all_panel_US %>% as.data.table()

## check for duplicates
all_panel_US[duplicated(all_panel_US[, c("start_date", "census_block_group", "region", "iso_country_code")]), ] %>%
  arrange(region, census_block_group, start_date) %>%
  pull(region) %>%
  unique()
head(all_panel_US)

## sanity check
all_panel_US %>%
  filter(region == "wa") %>%
  group_by(census_block_group, start_date) %>%
  tally() %>%
  filter(n != 1)

all_panel_US %>%
  filter(grepl("^53033", census_block_group)) %>%
  group_by(census_block_group) %>%
  slice_max(number_devices_primary_daytime) %>%
  arrange(-number_devices_primary_daytime)

ggplot(all_panel_US %>% as_tibble() %>% filter(census_block_group == "530330224006")) +
  geom_line(aes(x = start_date, y = number_devices_primary_daytime))

range(all_panel_US$start_date) # "2018-01-01" "2022-04-18"

all_panel_US_dt <- lazy_dt(all_panel_US)
## check for duplicates
all_panel_US_dt %>%
  group_by(region, iso_country_code, census_block_group, start_date) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

# directory for storing raw mobility data outputs
dir <- "~/Seattle_Flu_Study/Seattle_SG_Mobility/"
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Seattle_Flu_Study/Seattle_SG_Mobility/"
write_rds(all_panel_US, file = paste0(dir, "SG_data/Safegraph_daily_panel_devices_by_CBG_2018_to_2022.rds")) # daily panel size by CBG
# all_panel_US = read_rds(paste0(dir,"SG_data/Safegraph_daily_panel_devices_by_CBG_2018_to_2022.rds"))

## aggregate by county
all_panel_US_by_county <-
  all_panel_US %>%
  group_by(state_fips, county_fips, start_date) %>%
  summarize(
    county_devices_residing = sum(number_devices_residing),
    county_devices_daytime = sum(number_devices_primary_daytime)
  ) %>%
  ungroup()

ggplot(all_panel_US_by_county %>% as_tibble() %>% filter(state_fips == "53" & county_fips == "033")) +
  geom_line(aes(x = start_date, y = county_devices_daytime)) +
  geom_vline(aes(xintercept = as.Date("2018-05-01")))

all_panel_US_dt <- lazy_dt(all_panel_US_by_county)

## check for duplicates
all_panel_US_dt %>%
  group_by(state_fips, county_fips, start_date) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

all_panel_US_by_county <- as_tibble(all_panel_US_by_county)
write_rds(all_panel_US_by_county, file = paste0(dir,"SG_data/Safegraph_daily_panel_devices_by_county_2018_to_2022.rds")) # daily panel size by county
# all_panel_US_by_county <- read_rds(paste0(dir,"SG_data/Safegraph_daily_panel_devices_by_county_2018_to_2022.rds"))

ggplot(all_panel_US_by_county %>% filter(state_fips == "53")) +
  geom_line(aes(x = start_date, y = county_devices_residing)) +
  facet_wrap(~county_fips, scale = "free_y")

ggplot(all_panel_US %>% as_tibble() %>% filter(state_fips == "53" & county_fips == "033" & census_block_group == "530330294061")) +
  geom_line(aes(x = start_date, y = number_devices_residing)) +
  facet_wrap(~census_block_group, scale = "free_y")

# remove big files we don't need anymore
# rm(all_panel_US)
# rm(panel_current)
rm(all_panel)

####################################################
## Normalization Stats
####################################################
norm_dir <- "/Volumes/My Passport for Mac/SafeGraph_Dewey/normalization_stats/" # path to directory of your normalization data

norm <- read_many_csvs(
  dir = norm_dir,
  recursive = TRUE,
  makedate = TRUE
)
range(norm$date) # "2018-01-01" "2022-10-02"

norm <- norm %>% distinct()

## check for duplicates
norm %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1) %>%
  pull(date) %>%
  range()
## two entries for each date in April 2022

norm %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1) %>%
  tail()

norm_stats_US <- norm %>%
  filter(iso_country_code == "US") %>%
  arrange(date) %>%
  as_tibble()

head(norm_stats_US)
range(norm_stats_US$date)

norm_stats_US_dt <- lazy_dt(norm_stats_US)

## check for duplicates
norm_stats_US_dt %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

## get rid of duplicates
norm_stats_US_dt <-
  norm_stats_US_dt %>%
  group_by(region, iso_country_code, date) %>%
  slice_max(total_visits) %>%
  slice_max(total_home_visits) %>%
  slice_max(total_devices_seen) %>%
  ungroup() %>%
  distinct()

## check for duplicates
norm_stats_US_dt %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

norm_stats_US <- as_tibble(norm_stats_US_dt)

write_rds(norm_stats_US, file = paste0(dir,"SG_data/Safegraph_daily_normalization_stats_2018_to_2022.rds"))
# norm_stats_US <- read_rds(paste0(dir,"SG_data/Safegraph_daily_normalization_stats_2018_to_2022.rds"))
rm(norm)

####################################################
## Visit Panel Data
####################################################

visit_panel_dir <- "/Volumes/My Passport for Mac/SafeGraph_Dewey/visit_panel_summary/" # path to directory of your visit panel data
visit_panel <- read_many_csvs(dir = visit_panel_dir)
head(visit_panel)

visit_panel$date <- as.Date(paste(visit_panel$starting_year, visit_panel$starting_month, visit_panel$starting_day, sep = "-"), format = "%Y-%m-%d")
range(visit_panel$date) # "2018-01-01" "2022-09-26"

visit_panel <- visit_panel %>% unique()

all_visit_panel <- visit_panel %>% arrange(date)

## filter out Canada
all_visit_panel_US <- all_visit_panel %>%
  filter(iso_country_code == "US") %>%
  as_tibble()

head(all_visit_panel_US)
range(all_visit_panel_US$starting_year)
range(all_visit_panel_US$date)

## april 2022 has two entries for each date
all_visit_panel_US %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1)

ggplot(all_visit_panel_US) +
  geom_line(aes(x = date, y = num_visits)) +
  facet_wrap(~region, scales = "free_y")

ggplot(all_visit_panel_US %>% filter(region == "wa")) +
  geom_line(aes(x = date, y = num_visits)) +
  geom_vline(aes(xintercept = as.Date("2018-05-01")))

all_visit_panel_US_dt <- lazy_dt(all_visit_panel_US)

# check for duplicates
all_visit_panel_US_dt %>%
  group_by(date, region, iso_country_code) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

## get rid of duplicates
all_visit_panel_US_dt <-
  all_visit_panel_US_dt %>%
  group_by(region, iso_country_code, date) %>%
  slice_max(num_visits) %>%
  slice_max(num_unique_visitors) %>%
  ungroup() %>%
  distinct()

## check for duplicates
all_visit_panel_US_dt %>%
  group_by(region, iso_country_code, date) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

all_visit_panel_US <- as_tibble(all_visit_panel_US_dt)

write_rds(all_visit_panel_US, file = paste0(dir, "SG_data/Safegraph_daily_visit_panel_2018_to_2022.rds"))
# all_visit_panel_US <- read_rds(paste0(dir,"SG_data/Safegraph_daily_visit_panel_2018_to_2022.rds"))
# rm(visit_panel)
unique(all_visit_panel_US$region)

####################################################
## WEEKLY PATTERNS
####################################################

patterns_dir <- "/Volumes/My Passport for Mac/SafeGraph_Dewey/patterns_wa_state" # path to directory of your weekly patterns data

## read patterns data for king county
df <- read_many_patterns(
  dir = patterns_dir,
  silent = F,
  gen_fips = T,
  select = c(
    "placekey", "parent_placekey", "naics_code", "category_tags",
    "poi_cbg", "date_range_start", "date_range_end", "visits_by_day",
    "visitor_home_cbgs",
    "visitor_home_aggregation", "visitor_country_of_origin",
    "raw_visit_counts", "raw_visitor_counts"
  ),
  filter = 'state_fips==53 & county_fips =="033"'
)

# check for duplicates
df %>%
  group_by(date_range_start, placekey) %>%
  tally() %>%
  ungroup() %>%
  filter(n != 1) %>%
  distinct(date_range_start)

df_filt <- df %>%
  filter(!is.na(poi_cbg)) %>%
  as.data.table()
df_filt[, start_date := fasttime::fastPOSIXct(substr(date_range_start, 0, 10), tz = "UTC")]
df_filt <- df_filt %>% distinct()

# check for duplicates
df_filt %>%
  group_by(start_date, placekey) %>%
  tally() %>%
  ungroup() %>%
  filter(n != 1) %>%
  distinct(start_date)

df_filt <- df_filt %>% as_tibble()
range(df_filt$start_date)
head(df_filt)

# check for duplicates
df_filt %>%
  group_by(start_date, placekey) %>%
  tally() %>%
  ungroup() %>%
  filter(n != 1) %>%
  distinct(start_date)

write_rds(x = df_filt, paste0(dir,"SG_data/SG_patterns_raw_2018_01_01_to_2022_09_26.rds"))
# df_filt <- read_rds(paste0(dir,"SG_data/SG_patterns_raw_2018_01_01_to_2022_09_26.rds")) %>% as_tibble()

####################################################
# Filter out POIs that only occasionally appear in the dataset
####################################################
range(df_filt$start_date) # "2018-01-01 UTC" "2022-09-26 UTC"

df_filt_dt <- lazy_dt(df_filt)

sort(unique(df_filt$start_date))

## initially limit to late Oct 2018 onwards
all_df_filt <- df_filt_dt %>%
  unique() %>%
  # filter(start_date >= as.Date("2018-06-01")) %>%
  filter(start_date >= as.Date("2018-10-29")) %>%
  dplyr::select(-contains("date_range")) %>%
  as_tibble()

head(all_df_filt)
range(as.Date(all_df_filt$start_date)) # "2018-10-29" "2022-09-26"
length(unique(all_df_filt$start_date)) # 205

## number of POIs
all_df_filt %>%
  pull(placekey) %>%
  unique() %>%
  length() # 55249

# sanity check; should be zero
rem <- all_df_filt %>%
  group_by(placekey, start_date) %>%
  tally() %>%
  filter(n != 1)
rem

## check for duplicates
rem %>%
  pull(placekey) %>%
  unique() %>%
  length()
# no double counted

appear_df <- all_df_filt %>%
  group_by(placekey) %>%
  tally()
length(unique(appear_df$placekey))
length(unique(all_df_filt$start_date)) # 205

ggplot(appear_df) +
  geom_density(aes(x = n)) +
  geom_vline(xintercept = 195, lty = "dashed")

range(appear_df$n) # 1 205
total <- 205
100 * 200 / total # 98% of dates
100 * 190 / total # 93%
100 * 185 / total # 90%

appear_df %>%
  filter(n > 200) %>%
  pull(placekey) %>%
  unique() %>%
  length() # 14733

appear_df %>%
  filter(n > 190) %>%
  pull(placekey) %>%
  unique() %>%
  length() # 18295

appear_df %>%
  filter(n >= 185) %>%
  pull(placekey) %>%
  unique() %>%
  length() # 19817

## limit POIs to those present during 90% of dates since June 2018
keep <- appear_df %>%
  filter(n >= 185) %>%
  pull(placekey)
length(keep) # 19817

df_all_years <- all_df_filt %>% filter(placekey %in% keep)
range(unique(df_all_years$start_date))

####################################################
# Daily visits to cbgs in Seattle
####################################################
## note:  visits by day is originally in json format and collapsed by week
daily_visits_cbg_all_years <- expand_integer_json(df_all_years, "visits_by_day", index = "day", by = c("poi_cbg", "state_fips", "county_fips", "start_date"))[]
daily_visits_cbg_all_years$date <- as.Date(daily_visits_cbg_all_years$start_date) + (daily_visits_cbg_all_years$day - 1)

# check for duplicates
daily_visits_cbg_all_years %>%
  distinct(date) %>%
  group_by(date) %>%
  tally() %>%
  filter(n != 1)

write_rds(daily_visits_cbg_all_years, file = paste0(dir,"SG_data/SG_patterns_cbg_visits_by_day_2018_to_2022.rds"))

####################################################
# Daily visits to POIs by naics category
####################################################

visits <- df_all_years %>%
  dplyr::select(placekey, visits_by_day, naics_code, state_fips, county_fips, start_date) %>%
  filter(county_fips == "033" & !is.na(naics_code))

# check for duplicates
visits %>%
  group_by(placekey, start_date) %>%
  tally() %>%
  filter(n != 1)

daily_visits_naics_all_years <- expand_integer_json(visits, expand = "visits_by_day", index = "day", by = c("naics_code", "start_date"))[]
daily_visits_naics_all_years$date <- as.Date(daily_visits_naics_all_years$start_date) + (daily_visits_naics_all_years$day - 1)
daily_visits_naics_all_years$date <- as.Date(daily_visits_naics_all_years$date)
length(unique(daily_visits_naics_all_years$naics_code)) # number of categories in King Co

# check for duplicates
daily_visits_naics_all_years %>%
  group_by(naics_code, date) %>%
  tally() %>%
  filter(n != 1)

# import names of NAICS codes
naics <- openxlsx::read.xlsx("1_Seattle_Mobility_Data/mobility_data/2-6 digit_2017_Codes.xlsx")[, 1:2]
head(naics)
names(naics) <- c("naics_code", "industry")
naics$naics_code <- as.integer(naics$naics_code)
naics <- naics %>% filter(!is.na(naics_code))

daily_visits_naics_all_years$naics_code <- as.integer(daily_visits_naics_all_years$naics_code)

daily_visits_naics_all_years <- left_join(daily_visits_naics_all_years, naics %>% distinct(), by = "naics_code") %>%
  mutate(industry = ifelse(naics_code == 532299, "All Other Consumer Goods Rental", industry))

sort(unique(daily_visits_naics_all_years$industry))

daily_visits_naics_all_years$industry <- trimws(daily_visits_naics_all_years$industry, which = "right")

# check for duplicates
daily_visits_naics_all_years %>%
  group_by(date, naics_code, industry) %>%
  tally() %>%
  filter(n != 1)

daily_visits_naics_all_years %>%
  distinct(naics_code, industry) %>%
  group_by(naics_code) %>%
  tally() %>%
  filter(n != 1)

daily_visits_naics_all_years %>%
  distinct(naics_code, industry) %>%
  group_by(industry) %>%
  tally() %>%
  filter(n != 1)

sort(unique(daily_visits_naics_all_years$industry))
sort(unique(daily_visits_naics_all_years$naics_code))

daily_visits_naics_all_years %>%
  distinct(naics_code, industry) %>%
  group_by(industry) %>%
  tally() %>%
  filter(n != 1)

daily_visits_naics_all_years %>%
  filter(industry == "All Other Consumer Goods Rental") %>%
  distinct(naics_code)

## spot fix
daily_visits_naics_all_years <- daily_visits_naics_all_years %>% mutate(naics_code = ifelse(naics_code == 532289, 532299, naics_code))

daily_visits_naics_all_years %>%
  group_by(naics_code, industry, date) %>%
  tally() %>%
  filter(n != 1)

## aggregate daily visits by naics code
daily_visits_naics_all_years <- daily_visits_naics_all_years %>%
  group_by(industry, naics_code, start_date, date, day) %>%
  summarize(visits_by_day = sum(visits_by_day)) %>%
  ungroup()

## check for duplicates
daily_visits_naics_all_years %>%
  group_by(naics_code, industry, date) %>%
  tally() %>%
  filter(n != 1)

daily_visits_naics_all_years <- daily_visits_naics_all_years %>% as_tibble()
write_rds(daily_visits_naics_all_years, file = paste0(dir,"SG_data/SG_patterns_naics_visits_by_day_2018_to_2022.rds"))

####################################################
# Weekly visitors by visitor home cbg
####################################################

weekly_visitors_cbg_all_years <- expand_cat_json(df_all_years, "visitor_home_cbgs", index = "home_cbg", by = c("poi_cbg", "start_date"))[]
head(weekly_visitors_cbg_all_years)
names(weekly_visitors_cbg_all_years)[4] <- "weekly_visits"

write_rds(weekly_visitors_cbg_all_years, file = paste0(dir,"SG_data/SG_patterns_weekly_visitor_home_cbg_all_years_2018_to_2022.rds"))

####################################################
# Datasets for Scaling Visits
####################################################

####################################################
## home panel sizes by cbg, county, state, and national
####################################################
head(all_panel_US_by_county)
names(all_panel_US_by_county)[4] <- "county_home_panel_devices"

home_panel_state <-
  all_panel_US_by_county %>%
  group_by(state_fips, start_date) %>%
  summarize(sample_size_state_month = sum(county_home_panel_devices)) %>%
  ungroup()

home_panel_US <-
  all_panel_US_by_county %>%
  group_by(start_date) %>%
  summarize(sample_size_US_month = sum(county_home_panel_devices)) %>%
  ungroup()

home_panel_pop_size <-
  left_join(all_panel_US_by_county %>% dplyr::select(-county_devices_daytime), home_panel_state, by = c("start_date", "state_fips")) %>%
  left_join(home_panel_US, by = "start_date") %>%
  rename(sample_size_county_month = county_home_panel_devices)

home_panel_cbg <- all_panel_US %>%
  as_tibble() %>%
  dplyr::select(
    start_date, census_block_group,
    state_fips, county_fips, number_devices_residing
  )

head(home_panel_cbg)
names(home_panel_cbg)[5] <- c("home_panel_poi_cbg")

head(home_panel_pop_size)
names(home_panel_pop_size)[4:6] <- c("home_panel_poi_county", "home_panel_poi_state", "home_panel_poi_US")

home_panel_pop_size <- left_join(home_panel_cbg %>% filter(start_date >= as.Date("2018-01-01")),
  home_panel_pop_size,
  by = c("state_fips", "county_fips", "start_date")
)

####################################################
## add fips code to visit panel data
####################################################
all_visit_panel_US$start_date <- as.Date(paste(all_visit_panel_US$starting_year, all_visit_panel_US$starting_month, all_visit_panel_US$starting_day, sep = "-"), format = "%Y-%m-%d")

ggplot(all_visit_panel_US %>% filter(!(region %in% c("ALL", "ALL_US", "ALL_CA")))) +
  geom_line(aes(x = date, y = num_visits, color = region)) +
  facet_wrap(~region, scales = "free_y")

all_visit_panel_US_summary <-
  all_visit_panel_US %>%
  filter(!(region %in% c("ALL_US", "ALL", "ALL_CA"))) %>%
  mutate(visits_per_visitor_state = num_visits / num_unique_visitors) %>%
  mutate(region = toupper(region)) %>%
  ungroup()

all_visit_panel_US_summary <- left_join(all_visit_panel_US_summary,
  fips_codes[, c("state", "state_code", "state_name")] %>%
    distinct(),
  by = c("region" = "state")
)

####################################################
## combine home panel and visitor panel data
####################################################

home_and_visitor_panel <- left_join(home_panel_pop_size, all_visit_panel_US_summary,
  by = c("state_fips" = "state_code", "start_date")
) %>%
  dplyr::select(
    census_block_group, state_fips, county_fips,
    region, iso_country_code, start_date,
    home_panel_poi_cbg, home_panel_poi_county,
    home_panel_poi_state, home_panel_poi_US,
    num_visits, num_unique_visitors, visits_per_visitor_state
  )

## remove big files we don't need anymore
rm(all_panel_US)
rm(all_panel_US_by_county)
rm(df_all_years)
rm(df_filt)
rm(all_df_filt)

home_and_visitor_panel$census_block_group <- stringr::str_pad(home_and_visitor_panel$census_block_group, width = 12, pad = 0, side = "left")
home_and_visitor_panel <- home_and_visitor_panel %>% distinct()
range(home_and_visitor_panel$start_date)

home_and_visitor_panel_dt <- lazy_dt(home_and_visitor_panel)

# check for duplicates
home_and_visitor_panel_dt %>%
  group_by(start_date, census_block_group) %>%
  tally() %>%
  filter(n != 1) %>%
  as_tibble()

home_and_visitor_panel <- home_and_visitor_panel %>% as_tibble()
write_rds(home_and_visitor_panel, file = paste0(dir, "SG_data/home_and_visitor_panel_all_years.rds"))
