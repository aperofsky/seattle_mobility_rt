library(readr)
library(dplyr)
library(zoo)
library(tidyr)

#######################################################################################
## Combine mobility and Rt data into one data frame for statistical analyses
########################################################################################

# mobility data
combined_mob <- read_rds("1_Seattle_Mobility_Data/mobility_data/mobility_metrics_for_epidemia.rds") %>% as_tibble()
names(combined_mob)
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

## use data that have been smoothed with 2 week moving averages
combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

# Rt estimates
rt_df <- read_csv("2_Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

rt_df_wide <- rt_df %>%
  dplyr::select(date, median, organism) %>%
  distinct() %>%
  pivot_wider(names_from = "organism", values_from = "median") %>%
  dplyr::filter(date > as.Date("2018-11-20"))

names(rt_df_wide)[2:length(names(rt_df_wide))]

## rename pathogen columns
names(rt_df_wide)[names(rt_df_wide) %in% c(
  "RSV A", "RSV B", "IAV/H1N1", "IAV/H3N2", "IBV", "EV", "RV", "HPIV-3 + HPIV-4",
  "HPIV-1 + HPIV-2", "hMPV", "Adenovirus", "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63",
  "SARS-CoV-2"
)] <- c(
  "rsv_a",
  "rsv_b",
  "h1n1",
  "h3n2",
  "ivb",
  "entero",
  "rhino",
  "hpiv_3_4",
  "hpiv_1_2",
  "hmpv",
  "adeno",
  "scov_229E_OC43",
  "scov_HKU1_NL63",
  "covid"
)

names(rt_df_wide)

# combine 2 datasets
combined <- left_join(rt_df_wide,
  combined_mob_avg,
  by = "date"
) %>%
  filter(date < as.Date("2022-06-01")) #limit to before June 2022 because % staying home ends in late May 2022


## df with raw Rt values
# combined = left_join(rt_df_wide,
#                      combined_mob,by="date")

write_rds(combined,file="4_Block_Bootstrap_Mobility_vs_Rt/1_snowstorm_rolling_window_block_bootstrap/combined_rt_mobility_15day_mv_avg.rds")
write_rds(combined,file="4_Block_Bootstrap_Mobility_vs_Rt/2_post_2019_rolling_window_block_bootstrap/combined_rt_mobility_15day_mv_avg.rds")
write_rds(combined,file="3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds")
