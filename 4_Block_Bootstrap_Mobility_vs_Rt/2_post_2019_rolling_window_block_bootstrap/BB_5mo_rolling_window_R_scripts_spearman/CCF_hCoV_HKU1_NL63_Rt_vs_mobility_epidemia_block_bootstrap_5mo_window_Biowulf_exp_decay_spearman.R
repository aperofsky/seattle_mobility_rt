setwd("//home//perofskyamc//SFS_Rt_Block_Bootstrap")
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(tseries)
source("utils.R")
source("ts_block_bootstrap_sp.R")

########################################################################################
## import data
########################################################################################

combined <- read_rds("combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()

###########################################
## Actual CCF Loop
###########################################

metrics <- c(
  "within_neighborhood_movement",
  "within_city_movement",
  "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  "groceries_and_pharmacies",
  "transit",
  "religious_orgs",
  "child_day_care",
  "elementary_and_secondary_schools",
  "colleges"
)

combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined[2][is.na(combined[2])] <- 0

comb_weekly <- combined %>%
  mutate_at(vars(all_of(metrics)), ~ scale(.x) %>% as.vector()) %>%
  filter(epi_date >= as.Date("2019-06-01") & epi_date < as.Date("2020-08-01")) %>%
  group_by(epi_date) %>%
  summarize_at(c("scov_HKU1_NL63", metrics), ~ mean(.x, na.rm = T))
comb_weekly[2][is.na(comb_weekly[2])] <- 0

actual_data_and_perm <- ts_block_bootstrap_sp(df1=comb_weekly,pathogen="scov_HKU1_NL63",window=10,l=0.077,perms=1000)

save(actual_data_and_perm, file = "//data//perofskyamc//hCoV_HKU1_NL63_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData")
