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
Sys.setenv(CENSUS_KEY = "your key")
# Reload .Renviron
readRenviron("~/.Renviron")
# Check to see that the expected key is output in your R console
Sys.getenv("CENSUS_KEY")
####################################################
## Census Data
####################################################

# The first two digits of a Census Block Group represent the state_fips code and the next three digits represent a county_fips code
census_data <- read.csv("1_Seattle_Mobility_Data/mobility_data/cbg_fips_codes_2019.csv") # from SafeGraph's open census data
head(census_data)

census_data <- census_data %>%
  mutate(
    state_fips = stringr::str_pad(state_fips, width = 2, side = "left", pad = "0"),
    county_fips = stringr::str_pad(county_fips, width = 3, side = "left", pad = "0")
  ) %>%
  mutate(fips = paste(state_fips, county_fips, sep = ""))
head(census_data)
census_data %>% filter(state == "WA")

## get census block group population sizes
pop_csv <- read.csv("1_Seattle_Mobility_Data/mobility_data/cbg_b00.csv") # from SafeGraph's open census data

cbg_pop_sizes <- pop_csv %>%
  dplyr::select(census_block_group, B00001e1) %>%
  rename(cbg_pop = B00001e1)
cbg_pop_sizes[, c("state_fips", "county_fips")] <- fips_from_cbg(cbg_pop_sizes$census_block_group)

census_info <- left_join(cbg_pop_sizes, census_data, by = c("state_fips", "county_fips"))
census_info$census_block_group <- as.character(census_info$census_block_group)
census_info$census_block_group <- stringr::str_pad(census_info$census_block_group, width = 12, pad = 0, side = "left")

tract_puma <- read.table("1_Seattle_Mobility_Data/mobility_data/tract_puma_mapping.txt", sep = ",", header = T)[, c(2:3)]
tract_puma$residence_census_tract <- as.character(tract_puma$residence_census_tract)


acs_year <- 2019
seattle <- "033"

# Population of Seattle Area
norm_seattle_pop <-
  getCensus(
    name = "acs/acs1",
    vintage = acs_year,
    region = "county:*",
    regionin = "state:53",
    vars = "B01003_001E"
  ) %>%
  filter(county %in% seattle) %>%
  pull(B01003_001E) # This is an ACS code for population counts.
norm_seattle_pop
norm_seattle_pop <- 2252782

## pop size WA state
norm_wa_pop <- getCensus(
  name = "acs/acs1",
  vintage = acs_year,
  region = "state:53",
  vars = "B01003_001E"
) %>%
  pull(B01003_001E)

#####################################
# Data sets for Scaling Visits
#####################################
## home panel
dir <- "/Users/perofskyamc/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/NIH_Laptop_Updates_Post_Damage/Documents/Seattle_Flu_Study/Seattle_SG_Mobility/"
home_and_visitor_panel <- read_rds(paste0(dir, "SG_data/home_and_visitor_panel_all_years.rds"))

head(home_and_visitor_panel)
range(home_and_visitor_panel$start_date) # "2018-01-01" "2022-09-26"

home_and_visitor_panel <- home_and_visitor_panel %>% as.data.table()

# check for duplicates
home_and_visitor_panel[, .N, by = .(census_block_group, start_date)][N > 1]

home_and_visitor_panel <- as_tibble(home_and_visitor_panel)

## filter to King County, WA
home_and_visitor_panel2 <- home_and_visitor_panel %>%
  dplyr::select(
    start_date, census_block_group, state_fips, county_fips, home_panel_poi_cbg,
    home_panel_poi_county, home_panel_poi_state, home_panel_poi_US,
    visits_per_visitor_state
  ) %>%
  filter(state_fips == "53" & county_fips == "033") %>%
  distinct() %>%
  ungroup()

# check for duplicates
home_and_visitor_panel2 %>%
  group_by(census_block_group, start_date) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel2 %>%
  dplyr::select(-census_block_group, -home_panel_poi_cbg) %>%
  distinct() %>%
  group_by(
    start_date,
    state_fips, county_fips
  ) %>%
  tally() %>%
  filter(n != 1)

# make sure every cbg has an entry for every week
home_and_visitor_panel3 <- home_and_visitor_panel2 %>% pad(group = c("census_block_group"))

# check for duplicates
home_and_visitor_panel3 %>%
  group_by(census_block_group, start_date) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel3 %>% filter(is.na(home_panel_poi_cbg))

## a few CBGs have missing data on 6/25/2018, 10/14/2019, 6/222/2020 = fill in panel size for county, state, and US and visitors per visitor state
home_and_visitor_panel3 %>%
  distinct(start_date, home_panel_poi_county, home_panel_poi_state, home_panel_poi_US, visits_per_visitor_state) %>%
  filter(start_date %in% c(as.Date("2018-06-25"), as.Date("2019-10-14"), as.Date("2020-06-22")))

home_and_visitor_panel3 <-
  home_and_visitor_panel3 %>%
  mutate(
    state_fips = ifelse(is.na(state_fips), "53", state_fips),
    county_fips = ifelse(is.na(county_fips), "033", county_fips)
  ) %>%
  mutate(home_panel_poi_county = ifelse(start_date == as.Date("2018-06-25"), 83601, home_panel_poi_county)) %>%
  mutate(home_panel_poi_county = ifelse(start_date == as.Date("2019-10-14"), 90393, home_panel_poi_county)) %>%
  mutate(home_panel_poi_county = ifelse(start_date == as.Date("2020-06-22"), 74837, home_panel_poi_county)) %>%
  mutate(home_panel_poi_state = ifelse(start_date == as.Date("2018-06-25"), 263225, home_panel_poi_state)) %>%
  mutate(home_panel_poi_state = ifelse(start_date == as.Date("2019-10-14"), 307386, home_panel_poi_state)) %>%
  mutate(home_panel_poi_state = ifelse(start_date == as.Date("2020-06-22"), 287977, home_panel_poi_state)) %>%
  mutate(home_panel_poi_US = ifelse(start_date == as.Date("2018-06-25"), 15614433, home_panel_poi_US)) %>%
  mutate(home_panel_poi_US = ifelse(start_date == as.Date("2019-10-14"), 16145289, home_panel_poi_US)) %>%
  mutate(home_panel_poi_US = ifelse(start_date == as.Date("2020-06-22"), 15897548, home_panel_poi_US)) %>%
  mutate(visits_per_visitor_state = ifelse(start_date == as.Date("2018-06-25"), 7.407728, visits_per_visitor_state)) %>%
  mutate(visits_per_visitor_state = ifelse(start_date == as.Date("2019-10-14"), 9.997938, visits_per_visitor_state)) %>%
  mutate(visits_per_visitor_state = ifelse(start_date == as.Date("2020-06-22"), 7.562118, visits_per_visitor_state)) %>%
  distinct()

# check that all cbgs have data for these dates filled in now
home_and_visitor_panel3 %>%
  distinct(start_date, home_panel_poi_county, home_panel_poi_state, home_panel_poi_US, visits_per_visitor_state) %>%
  filter(start_date %in% c(as.Date("2018-06-25"), as.Date("2019-10-14"), as.Date("2020-06-22")))

home_and_visitor_panel3 <- home_and_visitor_panel3 %>% fill(everything(), .direction = "down")

# check for duplicates
data.table(home_and_visitor_panel3)[, .N, by = .(census_block_group, start_date)][N > 1]

# check for missing data
home_and_visitor_panel3[!complete.cases(home_and_visitor_panel3), ]

# check for duplicates
home_and_visitor_panel3 %>%
  dplyr::select(-census_block_group, -home_panel_poi_cbg) %>%
  distinct() %>%
  group_by(
    start_date,
    state_fips, county_fips
  ) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel3$home_panel_poi_US <- as.integer(home_and_visitor_panel3$home_panel_poi_US)
home_and_visitor_panel3$home_panel_poi_county <- as.integer(home_and_visitor_panel3$home_panel_poi_county)
home_and_visitor_panel3$home_panel_poi_cbg <- as.integer(home_and_visitor_panel3$home_panel_poi_cbg)
home_and_visitor_panel3$home_panel_poi_state <- as.integer(home_and_visitor_panel3$home_panel_poi_state)

## home panel size becomes really volatile in late 2021
## large increases with no concurrent increase in foot traffic to POIs; no explanations from SG
## apply panel size from November 2021 to remaining dates
home_and_visitor_panel4 <-
  home_and_visitor_panel3 %>%
  mutate(across(
    c(
      home_panel_poi_US, home_panel_poi_state,
      home_panel_poi_county, home_panel_poi_cbg
    ),
    ~ if_else(start_date > as.Date("2021-11-15") &
      start_date < as.Date("2022-08-22"), NA_integer_, .)
  )) %>%
  group_by(census_block_group, state_fips, county_fips) %>%
  fill(everything(), .direction = "down") %>%
  ungroup() %>%
  filter(!is.na(state_fips))

# check for duplicates
home_and_visitor_panel4 %>%
  dplyr::select(-census_block_group, -home_panel_poi_cbg) %>%
  distinct() %>%
  group_by(
    start_date,
    state_fips, county_fips
  ) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel4 %>%
  distinct(start_date, home_panel_poi_county) %>%
  group_by(home_panel_poi_county, start_date) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel4 %>%
  group_by(census_block_group, start_date) %>%
  tally() %>%
  filter(n != 1)

## remove big files we don't need anymore
rm(home_and_visitor_panel)
rm(home_and_visitor_panel2)
rm(home_and_visitor_panel3)
write_rds(home_and_visitor_panel4, file = paste0(dir, "SG_data/home_and_visitor_panel_2018_2022.rds"))
# home_and_visitor_panel4 <- read_rds(paste0(dir, "SG_data/home_and_visitor_panel_2018_2022.rds"))

####################################################
# Adjust POI Visitor Counts with Home Panel Sizes
####################################################

####################################################
### visits to POIs based on visitor home location
####################################################
weekly_visitors_cbg_all_years <- read_rds(paste0(dir, "SG_data/SG_patterns_weekly_visitor_home_cbg_all_years_2018_to_2022.rds"))
weekly_visitors_cbg_all_years <- as_tibble(weekly_visitors_cbg_all_years)
head(weekly_visitors_cbg_all_years)

weekly_visitors_cbg_all_years$start_date <- as.Date(weekly_visitors_cbg_all_years$start_date)
range(weekly_visitors_cbg_all_years$start_date) # "2018-10-29" "2022-09-26"

# check for duplicates
data.table(weekly_visitors_cbg_all_years)[, .N, by = .(home_cbg, poi_cbg, start_date)][N > 1]

unique(nchar(weekly_visitors_cbg_all_years$home_cbg))

# some visitors come from Canada CBGs
weekly_visitors_cbg_all_years %>%
  mutate(char = nchar(home_cbg)) %>%
  filter(char == 11)

## all POIs are in King Co
weekly_visitors_cbg_all_years %>%
  mutate(char = nchar(poi_cbg)) %>%
  filter(char == 11)

weekly_visitors_cbg_all_years %>%
  as_tibble() %>%
  filter(is.na(weekly_visits))

weekly_visitors_cbg_all_years$poi_cbg <- as.character(weekly_visitors_cbg_all_years$poi_cbg)
names(weekly_visitors_cbg_all_years)[1] <- "visitor_home_cbg"
weekly_visitors_cbg_all_years$visitor_home_cbg <- as.character(weekly_visitors_cbg_all_years$visitor_home_cbg)

## link county and state fips info to POI CBGs and visitor home CBGs
weekly_visitors_cbg_all_years <-
  left_join(weekly_visitors_cbg_all_years, census_info %>% dplyr::select(-class_code),
    by = c("poi_cbg" = "census_block_group")
  ) %>%
  dplyr::rename(
    poi_cbg_pop = cbg_pop,
    poi_state_fips = state_fips,
    poi_county_fips = county_fips,
    poi_fips = fips,
    poi_state = state,
    poi_county = county
  ) %>%
  left_join(census_info %>% dplyr::select(-class_code),
    by = c("visitor_home_cbg" = "census_block_group")
  ) %>%
  dplyr::rename(
    visitor_cbg_pop = cbg_pop,
    visitor_state_fips = state_fips,
    visitor_county_fips = county_fips,
    visitor_fips = fips,
    visitor_state = state,
    visitor_county = county
  )

# census tract
weekly_visitors_cbg_all_years$poi_tract <- substr(weekly_visitors_cbg_all_years$poi_cbg, 1, 11)
weekly_visitors_cbg_all_years$visitor_tract <- substr(weekly_visitors_cbg_all_years$visitor_home_cbg, 1, 11)

names(weekly_visitors_cbg_all_years)
nchar(unique(weekly_visitors_cbg_all_years$visitor_fips)) %>% unique()
sort(unique(weekly_visitors_cbg_all_years$visitor_state))

## add census tract
weekly_visitors_cbg_all_years <-
  left_join(weekly_visitors_cbg_all_years,
    tract_puma,
    by = c("poi_tract" = "residence_census_tract")
  ) %>%
  dplyr::rename(puma_poi = puma) %>%
  left_join(tract_puma, by = c("visitor_tract" = "residence_census_tract")) %>%
  dplyr::rename(visitor_puma = puma)

tail(weekly_visitors_cbg_all_years)

# check for duplicates
data.table(weekly_visitors_cbg_all_years)[, .N, by = .(visitor_home_cbg, poi_cbg, start_date)][N > 1]

# north and south king county PUMAs
south_king <- c("11612", "11614", "11613", "11615", "11611", "11610")
north_king <- c("11606", "11601", "11602", "11607", "11616", "11608", "11609", "11604", "11603", "11605")

## add whether poi is south or north king co
weekly_visitors_cbg_all_years <-
  weekly_visitors_cbg_all_years %>%
  mutate(
    poi_geo = case_when(
      puma_poi %in% south_king ~ "South King",
      puma_poi %in% north_king ~ "North King"
    ),
    visitor_geo = case_when(
      visitor_puma %in% south_king ~ "South King",
      visitor_puma %in% north_king ~ "North King"
    )
  )

df <- weekly_visitors_cbg_all_years %>%
  group_by(start_date) %>%
  summarize(total = sum(weekly_visits))

ggplot(df) +
  geom_line(aes(x = start_date, y = total)) +
  geom_vline(xintercept = as.Date("2022-06-15"), lty = "dashed")

## limit POIs to King County, WA
seattle_weekly_visitors_cbg_all_years <- weekly_visitors_cbg_all_years %>% filter(poi_fips == "53033")
unique(seattle_weekly_visitors_cbg_all_years$poi_geo)
unique(seattle_weekly_visitors_cbg_all_years$visitor_geo)

# join patterns dataset to home and visit panel data
seattle_weekly_visitors_cbg_all_years <-
  left_join(seattle_weekly_visitors_cbg_all_years,
    home_and_visitor_panel4 %>% ungroup(),
    by = c(
      "poi_state_fips" = "state_fips", "poi_county_fips" = "county_fips",
      "start_date", "poi_cbg" = "census_block_group"
    )
  )
tail(seattle_weekly_visitors_cbg_all_years)

## fill in missing data
seattle_weekly_visitors_cbg_all_years <- seattle_weekly_visitors_cbg_all_years %>%
  group_by(poi_state_fips, poi_county_fips, poi_cbg) %>%
  fill(contains(c("home_panel", "visits_per_visitor")), .direction = "down") %>%
  ungroup()

seattle_weekly_visitors_cbg_all_years$king_co_pop <- norm_seattle_pop # king county pop size
seattle_weekly_visitors_cbg_all_years$wa_state_pop <- norm_wa_pop # wa state pop size

## population size divided by devices captured by SG each month
seattle_weekly_visitors_cbg_all_years <-
  seattle_weekly_visitors_cbg_all_years %>%
  mutate(cbg_multiplier = poi_cbg_pop / home_panel_poi_cbg) %>% # cbg adjustment
  mutate(
    king_co_multiplier = king_co_pop / home_panel_poi_county, # county adjustment
    wa_state_multiplier = wa_state_pop / home_panel_poi_state # state adjustment
  ) %>%
  mutate(
    scaled_visits_cbg = weekly_visits * cbg_multiplier,
    scaled_visits_county = weekly_visits * king_co_multiplier,
    scaled_visits_state = weekly_visits * wa_state_multiplier
  )

## check for missing data
seattle_weekly_visitors_cbg_all_years %>%
  filter(is.na(cbg_multiplier)) %>%
  distinct(poi_cbg, poi_cbg_pop)

seattle_weekly_visitors_cbg_all_years %>%
  filter(is.na(king_co_multiplier)) %>%
  distinct(start_date, poi_cbg)

seattle_weekly_visitors_cbg_all_years %>%
  filter(is.na(wa_state_multiplier)) %>%
  distinct(start_date, poi_cbg)

## scale by visits per visitor (not used in manuscript)
seattle_weekly_visitors_cbg_all_years <-
  left_join(seattle_weekly_visitors_cbg_all_years,
    home_and_visitor_panel4 %>%
      dplyr::select(start_date, census_block_group, visits_per_visitor_state) %>%
      rename(SG_visits_per_visitor_poi_state = visits_per_visitor_state) %>%
      ungroup(),
    by = c("poi_cbg" = "census_block_group", "start_date")
  ) %>%
  group_by(poi_cbg) %>%
  fill(c(SG_visits_per_visitor_poi_state), .direction = "down") %>%
  ungroup() %>%
  mutate(
    scaled_visits_adj = scaled_visits_cbg * SG_visits_per_visitor_poi_state,
    scaled_visits_adj_county = scaled_visits_county * SG_visits_per_visitor_poi_state,
    scaled_visits_adj_state = scaled_visits_state * SG_visits_per_visitor_poi_state
  )

seattle_weekly_visitors_cbg_all_years %>% filter(is.na(scaled_visits_adj_state))
sort(unique(seattle_weekly_visitors_cbg_all_years$visitor_state))
head(seattle_weekly_visitors_cbg_all_years)
write_rds(seattle_weekly_visitors_cbg_all_years, file = paste0(dir, "SG_data/seattle_weekly_visitors_cbg_all_years_2018_2022.rds"))
# seattle_weekly_visitors_cbg_all_years <- read_rds(paste0(dir, "SG_data/seattle_weekly_visitors_cbg_all_years_2018_2022.rds"))

####################################################
## Seattle Mobility Network Maps: Figure 3 in main text
####################################################
# downloaded from: https://data-seattlecitygis.opendata.arcgis.com/datasets/SeattleCityGIS::census-block-groups-2010/about
## Seattle CBGs
cbg_gps <- read.csv("1_Seattle_Mobility_Data/mobility_data/Seattle_Census_Block_Groups_2010.csv", colClasses = c("GEOID10" = "character"))
head(cbg_gps)
unique(nchar(cbg_gps$GEOID10))

## limit POIs and visitor home locations to seattle CBGs
weekly_visitors_cbg_within_seattle <-
  seattle_weekly_visitors_cbg_all_years %>%
  filter(poi_cbg %in% cbg_gps$GEOID10 & visitor_home_cbg %in% cbg_gps$GEOID10)
head(weekly_visitors_cbg_within_seattle)

weekly_visitors_cbg_within_seattle <-
  left_join(weekly_visitors_cbg_within_seattle,
    cbg_gps[, c("GEOID10", "INTPTLAT10", "INTPTLON10")],
    by = c("poi_cbg" = "GEOID10")
  ) %>%
  rename(poi_lat = INTPTLAT10, poi_lon = INTPTLON10) %>%
  left_join(cbg_gps[, c("GEOID10", "INTPTLAT10", "INTPTLON10")],
    by = c("visitor_home_cbg" = "GEOID10")
  ) %>%
  rename(visitor_lat = INTPTLAT10, visitor_lon = INTPTLON10) %>%
  mutate(home_location = if_else(poi_cbg != visitor_home_cbg, "non_resident", "resident"))

## requires API key
king <- get_acs(
  state = "WA", county = "King", geography = "block group",
  variables = "B01003_001", geometry = TRUE, year = 2019
)

seattle_map <- king %>% filter(GEOID %in% cbg_gps$GEOID10)

tail(weekly_visitors_cbg_within_seattle)
range(weekly_visitors_cbg_within_seattle$scaled_visits_adj, na.rm = T)
range(weekly_visitors_cbg_within_seattle$scaled_visits_adj_county, na.rm = T)

## remove CBGs with < 5 visits to another CBG to reduce noise; 4 visits could be 2-4 visits
weekly_visitors_cbg_within_seattle <- weekly_visitors_cbg_within_seattle %>% filter(!is.na(scaled_visits_cbg) & weekly_visits > 4)

range(weekly_visitors_cbg_within_seattle$scaled_visits_cbg, na.rm = T)

write_rds(weekly_visitors_cbg_within_seattle, file = paste0(dir, "SG_data/weekly_visitors_cbg_within_seattle_2018_2022.rds"))
# weekly_visitors_cbg_within_seattle <- read_rds(paste0(dir, "SG_data/weekly_visitors_cbg_within_seattle_2018_2022.rds"))

range(weekly_visitors_cbg_within_seattle$scaled_visits_adj)

seattle_map$estimate <- as.numeric(seattle_map$estimate)

## check for duplicates
weekly_visitors_cbg_within_seattle %>%
  group_by(visitor_home_cbg, poi_cbg, start_date) %>%
  tally() %>%
  filter(n != 1)

## switch from scaled_visits_adj to scaled_visits_county because those counts are more stable
## change in panel in May 2022 causes big drop in total visits
range(weekly_visitors_cbg_within_seattle$start_date)

weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2019-02-11") & home_location == "non_resident") %>%
  summarize(net_movement = sum(scaled_visits_county)) #  984572

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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-02-11") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Feb 2019 Snowstorm") +
  theme(plot.title = element_text(hjust = 0.5)) +
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

mean(cbg_degree$degree) # 16.67
median(cbg_degree$degree) # 11
range(cbg_degree$degree) # 1 184

feb_deg <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median",
  rug = F, binwidth = 2,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 190)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))
feb_deg

# log scale on x-axis
feb_deg_log <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(0, 120))
feb_deg_log

## weighted degree
cbg_w_degree <- strength(mygraph, loops = F) %>% as.data.frame()
colnames(cbg_w_degree) <- "degree"
mean(cbg_w_degree$degree) 
median(cbg_w_degree$degree)
range(cbg_w_degree$degree)
feb_w_deg <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 500,
  fill = "#00AFBB"
) +
  xlab("Weighted Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 66000)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))
feb_w_deg

feb_w_deg_log <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75))
# expand_limits(x = c(0, 120))
feb_w_deg_log # similar shape to integer degree

feb_2019_map_w_deg <- plot_grid(feb_2019_map, feb_deg, nrow = 2, rel_heights = c(6, 2))
feb_2019_map_ww_deg <- plot_grid(feb_2019_map, feb_w_deg, nrow = 2, rel_heights = c(6, 2))
feb_2019_map_w_deg_log <- plot_grid(feb_2019_map, feb_deg_log, nrow = 2, rel_heights = c(6, 2), labels = c("A", "B"))


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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2019-07-22") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("July 2019") +
  theme(plot.title = element_text(hjust = 0.5)) +
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

jul_2019_deg <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 2,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(1, 190)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

# log scale on x-axis
jul_2019_log <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(0, 120))
jul_2019_log

# weighted degree
cbg_w_degree <- strength(mygraph, loops = F) %>% as.data.frame()
colnames(cbg_w_degree) <- "degree"
mean(cbg_w_degree$degree)
median(cbg_w_degree$degree)
range(cbg_w_degree$degree)

jul_2019_w_deg <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 500,
  color = "black", fill = "#00AFBB"
) +
  xlab("Weighted Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 66000)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

jul_2019_w_deg_log <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75))
# expand_limits(x = c(0, 120))
jul_2019_w_deg_log

jul_2019_map
jul_2019_map_w_deg <- plot_grid(jul_2019_map, jul_2019_deg, nrow = 2, rel_heights = c(6, 2))
jul_2019_map_ww_deg <- plot_grid(jul_2019_map, jul_2019_w_deg, nrow = 2, rel_heights = c(6, 2))
jul_2019_map_w_deg_log <- plot_grid(jul_2019_map, jul_2019_log, nrow = 2, rel_heights = c(6, 2))

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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2020-03-16") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("March 2020 Lockdown") +
  theme(plot.title = element_text(hjust = 0.5)) +
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
mar_2020_deg <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = TRUE, binwidth = 2,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(1, 190)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

mar_2020_log <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(0, 120))
mar_2020_log

cbg_w_degree <- strength(mygraph, loops = F) %>% as.data.frame()
colnames(cbg_w_degree) <- "degree"
mean(cbg_w_degree$degree)
median(cbg_w_degree$degree)
range(cbg_w_degree$degree)

# degree histogram
mar_2020_w_deg <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 500,
  color = "black", fill = "#00AFBB"
) +
  xlab("Weighted Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 66000)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))
mar_2020_w_deg

mar_2020_w_deg_log <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75))
# expand_limits(x = c(0, 120))
mar_2020_w_deg_log

mar_2020_map_w_deg <- plot_grid(march_2020_map, mar_2020_deg, nrow = 2, rel_heights = c(6, 2))
mar_2020_map_w_deg

mar_2020_map_ww_deg <- plot_grid(march_2020_map, mar_2020_w_deg, nrow = 2, rel_heights = c(6, 2))
mar_2020_map_ww_deg

mar_2020_map_w_deg_log <- plot_grid(march_2020_map, mar_2020_log, nrow = 2, rel_heights = c(6, 2))
mar_2020_map_w_deg_log

####################################################
## Delta wave: July 2021
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(start_date == as.Date("2021-07-05") & home_location == "non_resident") %>%
  summarize(net_movement = sum(scaled_visits_county)) #721159

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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") &
        home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visits", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Delta Wave (July 2021)") +
  theme(plot.title = element_text(hjust = 0.5)) +
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
jul_2021_deg <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 2,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(1, 190)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

jul_2021_log <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(0, 120))
jul_2021_log

cbg_w_degree <- strength(mygraph, loops = F) %>% as.data.frame()
colnames(cbg_w_degree) <- "degree"
mean(cbg_w_degree$degree)
median(cbg_w_degree$degree)
range(cbg_w_degree$degree)

# weighted degree
jul_2021_w_deg <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 500,
  color = "black", fill = "#00AFBB"
) +
  xlab("Weighted Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 66000)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

jul_2021_w_deg_log <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75))
# expand_limits(x = c(0, 120))
jul_2021_w_deg_log

jul_2021_map_w_deg <- plot_grid(july_2021_map, jul_2021_deg, nrow = 2, rel_heights = c(6, 2))
jul_2021_map_w_deg

jul_2021_map_ww_deg <- plot_grid(july_2021_map, jul_2021_w_deg, nrow = 2, rel_heights = c(6, 2))
jul_2021_map_ww_deg

jul_2021_map_w_deg_log <- plot_grid(july_2021_map, jul_2021_log, nrow = 2, rel_heights = c(6, 2))
jul_2021_map_w_deg_log

####################################################
## Omicron wave Jan 2022 (avoid potential effects of Christmas/NYE on movement)
####################################################

weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "non_resident" &
    start_date == as.Date("2022-01-17")) %>%
  summarize(net_movement = sum(scaled_visits_county))#   703591

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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2022-01-17") & home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visitors", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Omicron Wave (January 2022)") +
  theme(plot.title = element_text(hjust = 0.5)) +
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
jan_2022_deg <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 2,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(1, 190)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

jan_2022_log <- ggpubr::gghistogram(cbg_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75)) +
  expand_limits(x = c(0, 120))
jan_2022_log

## weighted degree centrality
cbg_w_degree <- strength(mygraph, loops = F) %>% as.data.frame()
colnames(cbg_w_degree) <- "degree"
mean(cbg_w_degree$degree) 
median(cbg_w_degree$degree)
range(cbg_w_degree$degree)

jan_2022_w_deg <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F, binwidth = 500,
  color = "black", fill = "#00AFBB"
) +
  xlab("Weighted Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 66000)) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100))

jan_2022_w_deg_log <- ggpubr::gghistogram(cbg_w_degree,
  x = "degree",
  add = "median", rug = F,
  color = "black", fill = "#00AFBB"
) +
  xlab("Degree k") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 75))
# expand_limits(x = c(0, 120))
jan_2022_w_deg_log


jan_2022_map_w_deg <- plot_grid(jan_2022_map, jan_2022_deg, nrow = 2, rel_heights = c(6, 2))
jan_2022_map_w_deg

jan_2022_map_ww_deg <- plot_grid(jan_2022_map, jan_2022_w_deg, nrow = 2, rel_heights = c(6, 2))
jan_2022_map_ww_deg

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
    ), lwd = 0.4
  ) +
  geom_curve(
    data = weekly_visitors_cbg_within_seattle %>%
      filter(start_date == as.Date("2021-07-05") &
        home_location == "non_resident" & scaled_visits_county >= 700),
    aes(
      x = poi_lon, y = poi_lat, xend = visitor_lon, yend = visitor_lat,
      color = scaled_visits_county
    ), lwd = 0.8
  ) +
  scale_color_viridis(
    trans = "log", na.value = "white", option = "magma", name = "Number of Visits", direction = -1,
    breaks = my_breaks, labels = my_breaks, limits = c(first(my_breaks), last(my_breaks))
  ) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  scale_alpha(guide = "none") +
  scale_size(guide = "none") +
  ggtitle("Delta Wave (July 2021)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Net Movement = 721,159") +
  ylab(NULL)
july_2021_map

## maps only
map_combined <- plot_grid(
  feb_2019_map + theme(legend.position = "none"),
  jul_2019_map + theme(legend.position = "none"),
  march_2020_map + theme(legend.position = "none"),
  july_2021_map + theme(legend.position = "none"),
  jan_2022_map + theme(legend.position = "none"),
  nrow = 1
)
leg <- get_legend(july_2021_map + theme(legend.direction = "horizontal"))
map_combined2 <- plot_grid(map_combined, leg, nrow = 2, rel_heights = c(1, 0.1))
map_combined2

## maps with degree histogram
map_combined <- plot_grid(
  feb_2019_map_w_deg + theme(legend.position = "none"),
  jul_2019_map_w_deg + theme(legend.position = "none"),
  mar_2020_map_w_deg + theme(legend.position = "none"),
  jul_2021_map_w_deg + theme(legend.position = "none"),
  jan_2022_map_w_deg + theme(legend.position = "none"),
  nrow = 1
)
leg <- get_legend(july_2021_map + theme(legend.direction = "vertical"))
map_combined2 <- plot_grid(map_combined, leg, nrow = 1, rel_widths = c(1, 0.1))
map_combined2

## maps with weighted degree histogram
map_combined <- plot_grid(
  feb_2019_map_ww_deg + theme(legend.position = "none"),
  jul_2019_map_ww_deg + theme(legend.position = "none"),
  mar_2020_map_ww_deg + theme(legend.position = "none"),
  jul_2021_map_ww_deg + theme(legend.position = "none"),
  jan_2022_map_ww_deg + theme(legend.position = "none"),
  nrow = 1
)
leg <- get_legend(july_2021_map + theme(legend.direction = "vertical"))
map_combined2 <- plot_grid(map_combined, leg, nrow = 1, rel_widths = c(1, 0.1))
map_combined2

## maps with degree histogram on log scale = Figure 3
map_combined <- plot_grid(
  feb_2019_map_w_deg_log + theme(legend.position = "none"),
  jul_2019_map_w_deg_log + theme(legend.position = "none"),
  mar_2020_map_w_deg_log + theme(legend.position = "none"),
  jul_2021_map_w_deg_log + theme(legend.position = "none"),
  jan_2022_map_w_deg_log + theme(legend.position = "none"),
  nrow = 1
)
leg <- get_legend(july_2021_map + theme(legend.direction = "vertical"))
map_combined2 <- plot_grid(map_combined, leg, nrow = 1, rel_widths = c(1, 0.1))
map_combined2

## Figure 3
save_plot(map_combined2,
  filename = "1_Seattle_Mobility_Data/mobility_figures/fig_3_within_seattle_movement_up_to_Jan_2022_with_log_trans_deg_histogram.png",
  base_width = 16, base_height = 8
)

save_plot(map_combined2,
  filename = "figures/fig_3_within_seattle_movement_up_to_Jan_2022_with_log_trans_deg_histogram.png",
  base_width = 16, base_height = 8
)

# for conference abstract
# save_plot(map_combined2,
#           filename = "figures/fig_3_within_seattle_movement_up_to_Jan_2022_with_log_trans_deg_histogram_small.png",
#           base_width = 16, base_height = 5.5
# )

############################################################
## Large scale movements for Seattle and King Co.
############################################################

############################################################
## Within-neighborhood (CBG) movement
############################################################
## only seattle
within_cbg_seattle <-
  weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "resident") %>%
  group_by(start_date) %>%
  summarize(within_neighborhood_movement = sum(scaled_visits_state)) %>%
  ungroup()

## all of king co.
within_cbg_king <-
  seattle_weekly_visitors_cbg_all_years %>%
  filter(poi_fips == "53033") %>%
  filter(weekly_visits > 4) %>%
  # filter(home_location == "resident") %>%
  filter(visitor_home_cbg == poi_cbg) %>%
  group_by(start_date) %>%
  summarize(within_neighborhood_movement = sum(scaled_visits_state)) %>%
  ungroup()

## very similar trends between Seattle and King Co. as a whole
ggplot() +
  geom_line(aes(x = start_date, y = within_neighborhood_movement * 4),
    data = within_cbg_seattle %>% filter(start_date < as.Date("2022-07-01"))
  ) + # seattle only
  geom_line(aes(x = start_date, y = within_neighborhood_movement),
    data = within_cbg_king %>% filter(start_date < as.Date("2022-07-01")), color = "red"
  ) # king county wa

############################################################
## between-neighborhood movement (CBG)
############################################################
## only seattle
within_city_seattle <-
  weekly_visitors_cbg_within_seattle %>%
  filter(home_location == "non_resident") %>%
  group_by(start_date) %>%
  summarize(within_city_movement = sum(scaled_visits_state)) %>%
  ungroup()
range(within_city_seattle$within_city_movement)

# all of king co.
within_city_king <-
  seattle_weekly_visitors_cbg_all_years %>%
  filter(poi_fips == "53033") %>%
  filter(weekly_visits > 4) %>%
  filter(visitor_home_cbg != poi_cbg) %>%
  group_by(start_date) %>%
  summarize(within_city_movement = sum(scaled_visits_state)) %>%
  ungroup()
range(within_city_king$within_city_movement)

# volume of movement in king co as a whole picks up to a greater extent than Seattle after SAH lifted, but similar trends
ggplot() +
  geom_line(aes(x = start_date, y = 8 * within_city_movement),
    data = within_city_seattle %>% filter(start_date < as.Date("2022-07-01"))
  ) + # seattle only
  geom_line(aes(x = start_date, y = within_city_movement),
    data = within_city_king %>% filter(start_date < as.Date("2022-07-01")), color = "red"
  ) # king county wa

############################################################
## Inflow from other WA counties
############################################################
# only seattle
## keep visits >4 so that counties with fewer visits aren't excluded
county_visitors_out_of_county_within_state_seattle <-
  seattle_weekly_visitors_cbg_all_years %>%
  # filter(weekly_visits>4)%>%
  filter(poi_cbg %in% cbg_gps$GEOID10) %>%
  filter(!is.na(scaled_visits_adj)) %>%
  filter(poi_fips == "53033" & visitor_state_fips == "53") %>% # visitors from within WA
  filter(poi_fips != visitor_fips) %>%
  group_by(start_date, poi_fips, visitor_fips, visitor_state_fips, visitor_county_fips, visitor_county) %>%
  summarize(num_visits = sum(scaled_visits_state, na.rm = T)) %>%
  ungroup()

## King Co.
county_visitors_out_of_county_within_state_king <-
  seattle_weekly_visitors_cbg_all_years %>%
  # filter(weekly_visits>4)%>%
  filter(!is.na(scaled_visits_adj)) %>%
  filter(poi_fips == "53033" & visitor_state_fips == "53") %>% # visitors from within WA
  filter(poi_fips != visitor_fips) %>%
  group_by(start_date, poi_fips, visitor_fips, visitor_state_fips, visitor_county_fips, visitor_county) %>%
  summarize(num_visits = sum(scaled_visits_state, na.rm = T)) %>%
  ungroup()

# similar trends between Seattle vs all of King Co.
plot1 <- ggplot() +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_county), lwd = 0.8, data = county_visitors_out_of_county_within_state_king) +
  geom_line(aes(x = start_date, y = 3 * num_visits), color = "black", linetype = "dashed", lwd = 0.5, data = county_visitors_out_of_county_within_state_seattle) +
  facet_wrap(~visitor_county, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from WA Counties (SafeGraph)") +
  ylab("Number of Visitors")
plot1

# which counties have the most visitors to Seattle and King Co.?
county_visitors_out_of_county_within_state_king %>%
  group_by(visitor_county) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 30)

county_visitors_out_of_county_within_state_king %>%
  group_by(visitor_county) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 40)

county_visitors_out_of_county_within_state_seattle %>%
  group_by(visitor_county) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 40)

## get rid of counties with low numbers of visits (more volatile time series)
keep <- county_visitors_out_of_county_within_state_seattle %>%
  group_by(visitor_county) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  filter(mean >= 1000) %>%
  pull(visitor_county)
keep

## visitors to King Co. vs Seattle
plot2 <- ggplot() +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_county),
    lwd = 0.8,
    data = county_visitors_out_of_county_within_state_king %>% filter(visitor_county %in% keep)
  ) +
  geom_line(aes(x = start_date, y = 3 * num_visits),
    color = "black", linetype = "dashed", lwd = 0.5,
    data = county_visitors_out_of_county_within_state_seattle %>% filter(visitor_county %in% keep)
  ) +
  facet_wrap(~visitor_county, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from WA Counties (SafeGraph)") +
  ylab("Number of Visitors")
plot2

## visitors to King Co. vs Seattle = Pierce County only
plot2 <- ggplot() +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_county),
    lwd = 0.8,
    data = county_visitors_out_of_county_within_state_king %>% filter(visitor_county == "Pierce County")
  ) +
  geom_line(aes(x = start_date, y = 5 * num_visits),
    color = "black", linetype = "dashed", lwd = 0.5,
    data = county_visitors_out_of_county_within_state_seattle %>% filter(visitor_county == "Pierce County")
  ) +
  facet_wrap(~visitor_county, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from WA Counties (SafeGraph)") +
  ylab("Number of Visitors")
plot2

length(keep)
length(unique(county_visitors_out_of_county_within_state_seattle$visitor_county))

## aggregate visitors from all other WA counties
within_state_seattle <-
  county_visitors_out_of_county_within_state_seattle %>%
  filter(visitor_county %in% keep) %>%
  group_by(start_date) %>%
  summarize(within_state_movement = sum(num_visits, na.rm = T)) %>%
  ungroup()
tail(within_state_seattle)

within_state_king <-
  county_visitors_out_of_county_within_state_king %>%
  filter(visitor_county %in% keep) %>%
  group_by(start_date) %>%
  summarize(within_state_movement = sum(num_visits, na.rm = T)) %>%
  ungroup()
tail(within_state_king)
############################################################
## inflow from other US states
############################################################
# only seattle
county_visitors_outside_state_seattle <-
  seattle_weekly_visitors_cbg_all_years %>%
  filter(poi_cbg %in% cbg_gps$GEOID10) %>%
  filter(poi_fips == "53033" & visitor_state_fips != "53") %>%
  group_by(start_date, visitor_state_fips, visitor_state) %>%
  summarise(
    num_visits = sum(scaled_visits_state, na.rm = T)
  ) %>%
  ungroup()

# all of King Co.
county_visitors_outside_state_king <-
  seattle_weekly_visitors_cbg_all_years %>%
  filter(poi_fips == "53033" & visitor_state_fips != "53") %>%
  group_by(start_date, visitor_state_fips, visitor_state) %>%
  summarise(
    num_visits = sum(scaled_visits_state, na.rm = T)
  ) %>%
  ungroup()

county_visitors_outside_state_king %>%
  filter(visitor_state == "KS" & start_date > as.Date("2020-01-01")) %>%
  print(n = 30)

county_visitors_outside_state_king %>%
  filter(visitor_state == "FL" & start_date > as.Date("2022-01-01")) %>%
  print(n = 30)

## spot fix abberations
county_visitors_outside_state_king <-
  county_visitors_outside_state_king %>%
  mutate(num_visits = ifelse(visitor_state == "KS" &
    start_date >= as.Date("2020-02-24") &
    start_date <= as.Date("2020-03-30"), NA_integer_, num_visits)) %>%
  mutate(num_visits = ifelse(visitor_state == "FL" &
    start_date >= as.Date("2022-06-06") &
    start_date <= as.Date("2022-06-13"), NA_integer_, num_visits)) %>%
  group_by(visitor_state) %>%
  mutate(num_visits = forecast::na.interp(num_visits)) %>%
  ungroup()

## seattle has same issues with KS and FL on these dates
county_visitors_outside_state_seattle %>%
  filter(visitor_state == "KS" & start_date > as.Date("2020-01-01")) %>%
  print(n = 30)

county_visitors_outside_state_seattle %>%
  filter(visitor_state == "FL" & start_date > as.Date("2022-01-01")) %>%
  print(n = 30)

county_visitors_outside_state_seattle <-
  county_visitors_outside_state_seattle %>%
  mutate(num_visits = ifelse(visitor_state == "KS" &
    start_date >= as.Date("2020-02-24") &
    start_date <= as.Date("2020-03-30"), NA_integer_, num_visits)) %>%
  mutate(num_visits = ifelse(visitor_state == "FL" &
    start_date >= as.Date("2022-06-06") &
    start_date <= as.Date("2022-06-13"), NA_integer_, num_visits)) %>%
  group_by(visitor_state) %>%
  mutate(num_visits = forecast::na.interp(num_visits)) %>%
  ungroup()

# almost identical trends between Seattle and all of King Co.
ggplot() +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_state), lwd = 0.8, data = county_visitors_outside_state_king) +
  geom_line(aes(x = start_date, y = 2 * num_visits), color = "black", linetype = "dashed", lwd = 0.5, data = county_visitors_outside_state_seattle) +
  facet_wrap(~visitor_state, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from other states (SafeGraph)") +
  ylab("Number of Visitors")

county_visitors_outside_state_king %>%
  group_by(visitor_state) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 50)

county_visitors_outside_state_king %>%
  filter(start_date > as.Date("2021-07-01")) %>%
  group_by(visitor_state) %>%
  slice_max(num_visits) %>%
  arrange(start_date) %>%
  print(n = 40)

## filter out states with low number of visits (more volatile time series)
keep <- county_visitors_outside_state_king %>%
  group_by(visitor_state) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  filter(!(visitor_state %in% c("NE", "DC"))) %>%
  filter(mean > 900) %>%
  pull(visitor_state)
keep

# spot fix date with abnormally large spike for many states
county_visitors_outside_state_king <-
  county_visitors_outside_state_king %>%
  mutate(num_visits = ifelse(start_date == as.Date("2021-08-30"), NA_integer_, num_visits)) %>%
  group_by(visitor_state) %>%
  mutate(num_visits = forecast::na.interp(num_visits)) %>%
  ungroup()

county_visitors_outside_state_seattle <-
  county_visitors_outside_state_seattle %>%
  mutate(num_visits = ifelse(start_date == as.Date("2021-08-30"), NA_integer_, num_visits)) %>%
  group_by(visitor_state) %>%
  mutate(num_visits = forecast::na.interp(num_visits)) %>%
  ungroup()

ggplot(county_visitors_outside_state_king %>% filter(visitor_state %in% keep & start_date < as.Date("2022-07-01"))) +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_state), lwd = 0.8) +
  facet_wrap(~visitor_state, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from other states (SafeGraph)") +
  ylab("Number of Visitors")

# trends are almost identical
ggplot() +
  geom_line(aes(x = start_date, y = num_visits, color = visitor_state),
    lwd = 0.8,
    data = county_visitors_outside_state_king %>% filter(visitor_state %in% keep & start_date < as.Date("2022-07-01"))
  ) +
  geom_line(aes(x = start_date, y = 2 * num_visits),
    lwd = 0.8, color = "black", lty = "dashed",
    data = county_visitors_outside_state_seattle %>% filter(visitor_state %in% keep & start_date < as.Date("2022-07-01"))
  ) +
  facet_wrap(~visitor_state, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", legend.justification = 0.5) +
  ggtitle("Mobility In Flows from other states (SafeGraph)") +
  ylab("Number of Visitors")

county_visitors_outside_state_king %>%
  group_by(visitor_state) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 50)

county_visitors_outside_state_seattle %>%
  group_by(visitor_state) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  print(n = 50)

county_visitors_outside_state_king %>%
  filter(start_date > as.Date("2021-07-01")) %>%
  group_by(visitor_state) %>%
  slice_max(num_visits) %>%
  arrange(start_date) %>%
  print(n = 40)

keep <- county_visitors_outside_state_king %>%
  group_by(visitor_state) %>%
  summarize(mean = mean(num_visits, na.rm = T)) %>%
  arrange(-mean) %>%
  filter(!(visitor_state %in% c("NE", "DC"))) %>%
  filter(mean > 900) %>%
  pull(visitor_state)
keep

## aggregate visitors from all other states
out_of_state_king <-
  county_visitors_outside_state_king %>%
  filter(visitor_state %in% keep) %>%
  group_by(start_date) %>%
  summarize(out_of_state_movement = sum(num_visits, na.rm = T)) %>%
  ungroup()

out_of_state_seattle <-
  county_visitors_outside_state_seattle %>%
  filter(visitor_state %in% keep) %>%
  group_by(start_date) %>%
  summarize(out_of_state_movement = sum(num_visits, na.rm = T)) %>%
  ungroup()

names(out_of_state_king)[2] <- "out_of_state_movement_king"
names(out_of_state_seattle)[2] <- "out_of_state_movement_seattle"
names(within_state_king)[2] <- "within_state_movement_king"
names(within_state_seattle)[2] <- "within_state_movement_seattle"
names(within_city_king)[2] <- "within_city_movement_king"
names(within_city_seattle)[2] <- "within_city_movement_seattle"
names(within_cbg_king)[2] <- "within_neighborhood_movement_king"
names(within_cbg_seattle)[2] <- "within_neighborhood_movement_seattle"

## combine all metrics
all_between_movement <- plyr::join_all(list(
  out_of_state_king, out_of_state_seattle,
  within_state_king, within_state_seattle,
  within_city_king, within_city_seattle,
  within_cbg_king, within_cbg_seattle
), by = "start_date")

head(all_between_movement)
range(all_between_movement$start_date)

# check for duplicates
all_between_movement %>%
  group_by(start_date) %>%
  tally() %>%
  filter(n != 1)

range(all_between_movement$start_date)

write_rds(all_between_movement, file = paste0(dir, "SG_data/within_and_between_cbg_movement_indicators_up_to_2022_09_26.rds"))
# all_between_movement <- read_rds(paste0(dir, "SG_data/within_and_between_cbg_movement_indicators_up_to_2022_09_26.rds"))

all_between_movement_lim <- all_between_movement %>%
  dplyr::select(
    start_date, within_neighborhood_movement_seattle, within_city_movement_seattle,
    within_state_movement_king, out_of_state_movement_king
  ) %>%
  rename(
    within_neighborhood_movement = within_neighborhood_movement_seattle,
    within_city_movement = within_city_movement_seattle,
    within_state_movement = within_state_movement_king,
    out_of_state_movement = out_of_state_movement_king
  )
write_rds(all_between_movement_lim, file = paste0(dir, "SG_data/within_and_between_cbg_movement_indicators_up_to_2022_09_26_red.rds"))

############################################################
## visits to different categories of POIs
############################################################

daily_visits_naics_all_years <- read_rds(paste0(dir, "SG_data/SG_patterns_naics_visits_by_day_2018_to_2022.rds"))
head(daily_visits_naics_all_years)

## documented issue: parks in Seattle are "sinks" (artifically inflated number of pings)
daily_visits_naics_all_years %>% filter(industry %in% c("Nature Parks and Other Similar Institutions"))

# check for duplicates
daily_visits_naics_all_years %>%
  group_by(industry, date) %>%
  tally() %>%
  filter(n != 1)

home_and_visitor_panel4 %>%
  dplyr::select(-census_block_group, -home_panel_poi_cbg) %>%
  distinct() %>%
  group_by(
    start_date,
    state_fips, county_fips
  ) %>%
  tally() %>%
  filter(n != 1)

## combine weekly patterns dataset with home and visitor panel data
seattle_daily_visits_naics_all_years <-
  left_join(daily_visits_naics_all_years,
    home_and_visitor_panel4 %>%
      dplyr::select(-census_block_group, -home_panel_poi_cbg) %>%
      distinct(),
    by = "start_date"
  )

seattle_daily_visits_naics_all_years$king_co_pop <- norm_seattle_pop # king co. pop size
seattle_daily_visits_naics_all_years$wa_state_pop <- norm_wa_pop # wa state pop size

seattle_daily_visits_naics_all_years <-
  seattle_daily_visits_naics_all_years %>%
  mutate(
    king_co_multiplier = king_co_pop / home_panel_poi_county,
    wa_state_multiplier = wa_state_pop / home_panel_poi_state
  ) %>%
  # filter(weekly_visits>4)%>%
  mutate(
    scaled_visits_county = visits_by_day * king_co_multiplier,
    scaled_visits_state = visits_by_day * wa_state_multiplier
  ) %>%
  rename(SG_visits_per_visitor_poi_state = visits_per_visitor_state) %>%
  mutate(
    scaled_visits_adj_county = scaled_visits_county * SG_visits_per_visitor_poi_state, ## scale by visits per visitor (not used in manuscript)
    scaled_visits_adj_state = scaled_visits_state * SG_visits_per_visitor_poi_state
  )

seattle_daily_visits_naics_all_years$industry <- trimws(seattle_daily_visits_naics_all_years$industry, which = "right")

## fill in missing data
seattle_daily_visits_naics_all_years <- seattle_daily_visits_naics_all_years %>%
  group_by(industry, naics_code) %>%
  fill(everything(), .direction = "down") %>%
  ungroup()

seattle_daily_visits_naics_all_years %>%
  filter(date > as.Date("2021-01-01")) %>%
  slice_max(scaled_visits_adj_county, n = 10) # issue with nonresidential buildings

write_rds(seattle_daily_visits_naics_all_years, file = paste0(dir, "SG_data/seattle_daily_visits_naics.rds"))
# seattle_daily_visits_naics_all_years = read_rds(paste0(dir, "SG_data/seattle_daily_visits_naics.rds"))
