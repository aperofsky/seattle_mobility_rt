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

combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined$not_masking <- 100 - combined$mask_wearing_final

metrics <- c(
  "within_neighborhood_movement",
  "within_city_movement",
  "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  "groceries_and_pharmacies",
  "full_service_restaurants",
  "transit",
  "religious_orgs",
  "child_day_care",
  "oxford_stringency_index",
  "not_masking",
  "elementary_and_secondary_schools",
  "colleges"
)

comb_weekly <- combined %>%
  mutate_at(vars(all_of(metrics)), ~ scale(.x) %>% as.vector()) %>%
  filter(epi_date >= as.Date("2020-06-01") & epi_date < as.Date("2022-05-22")) %>%
  group_by(epi_date) %>%
  summarize_at(c("rhino", metrics), ~ mean(.x, na.rm = T))

actual_data_and_perm <- ts_block_bootstrap_sp(df1=comb_weekly,pathogen="rhino",window=10,l=0.077,perms=1000)

save(actual_data_and_perm, file = "//data//perofskyamc//rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData")
