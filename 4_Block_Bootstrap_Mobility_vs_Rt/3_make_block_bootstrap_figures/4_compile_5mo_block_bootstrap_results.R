library(dplyr)
library(readr)
library(cdcfluview)

########################################################################################
## Compile 5mo block bootstrap results
########################################################################################

dir <- "4_Block_Bootstrap_Mobility_vs_Rt/2_post_2019_rolling_window_block_bootstrap/BB_5mo_rolling_window_output_spearman/"

## adenovirus
load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
head(actual_data_and_perm)
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))

load(paste0(dir, "adeno_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi_spearman.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))

all_adeno <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
range(all_adeno$start_week)
unique(all_adeno$mobility_metric)
unique(all_adeno$pathogen)
########################################################################################
### covid
load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))

load(paste0(dir, "covid_mobility_CCF_5mo_sliding_window_actual_and_null_output_osi_spearman.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))

all_covid <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
unique(all_covid$mobility_metric)
unique(all_covid$pathogen)
########################################################################################
## seasonal CoV
load(paste0(dir, "hCoV_229E_OC43_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "scov_229E_OC43_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm

all_scov_229E_OC43 <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_scov_229E_OC43$mobility_metric)
unique(all_scov_229E_OC43$pathogen)

load(paste0(dir, "hCoV_HKU1_NL63_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "scov_HKU1_NL63_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm

all_scov_HKU1_NL63 <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_scov_HKU1_NL63$mobility_metric)
unique(all_scov_HKU1_NL63$pathogen)
########################################################################################
## hmpv
load(paste0(dir, "hmpv_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "hmpv_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm

all_hmpv <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_hmpv$mobility_metric)
unique(all_hmpv$pathogen)
########################################################################################
## HPIV
load(paste0(dir, "hpiv_3_4_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "hpiv_3_4_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm

all_hpiv <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
range(all_hpiv$start_week)
unique(all_hpiv$mobility_metric)
unique(all_hpiv$pathogen)

load(paste0(dir, "hpiv_1_2_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
all_hpiv_1_2 <- actual_data_and_perm
unique(all_hpiv_1_2$mobility_metric)
unique(all_hpiv_1_2$pathogen)
########################################################################################
## rhino
load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))

load(paste0(dir, "rhino_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))

all_rhino <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
sort(unique(all_rhino$mobility_metric))
unique(all_rhino$pathogen)
########################################################################################
## entero
load(paste0(dir, "entero_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "entero_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm %>% filter(mobility_metric %in% c("not_masking"))

load(paste0(dir, "entero_mobility_CCF_5mo_sliding_window_actual_and_null_output_masking_spearman.RData"))
actual_data_and_perm3 <- actual_data_and_perm %>% filter(mobility_metric %in% c("oxford_stringency_index"))

all_entero <- bind_rows(actual_data_and_perm1, actual_data_and_perm2, actual_data_and_perm3)
sort(unique(all_entero$mobility_metric))
unique(all_entero$pathogen)
########################################################################################
## rsv a
load(paste0(dir, "rsv_a_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
all_rsv_a <- actual_data_and_perm
unique(all_rsv_a$pathogen)
########################################################################################
## rsv b
load(paste0(dir, "rsv_b_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
actual_data_and_perm1 <- actual_data_and_perm

load(paste0(dir, "rsv_b_mobility_CCF_5mo_sliding_window_actual_and_null_output_rebound_2021_spearman.RData"))
actual_data_and_perm2 <- actual_data_and_perm

all_rsv_b <- bind_rows(actual_data_and_perm1, actual_data_and_perm2)
unique(all_rsv_b$mobility_metric)
unique(all_rsv_b$pathogen)
########################################################################################
### h1n1
load(paste0(dir, "h1n1_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
all_h1n1 <- actual_data_and_perm
unique(all_h1n1$pathogen)
########################################################################################
# ivb
load(paste0(dir, "ivb_mobility_CCF_5mo_sliding_window_actual_and_null_output_spearman.RData"))
all_ivb <- actual_data_and_perm
unique(all_ivb$pathogen)

########################################################################################
## combine data for all pathogens
########################################################################################
all_path_ccf <- bind_rows(
  all_entero,
  all_adeno, all_covid, all_hmpv, all_h1n1, all_ivb,
  all_hpiv, all_hpiv_1_2, all_rhino, all_rsv_a,
  all_rsv_b, all_scov_229E_OC43, all_scov_HKU1_NL63
)

all_path_ccf$pathogen <- as.factor(all_path_ccf$pathogen)
levels(all_path_ccf$pathogen)
levels(all_path_ccf$mobility_metric)
unique(all_path_ccf$pathogen)
all_path_ccf$pathogen <- factor(all_path_ccf$pathogen, levels = c(
  "h1n1", "ivb", "rsv_a",
  "rsv_b", "hmpv",
  "scov_229E_OC43",
  "scov_HKU1_NL63",
  "hpiv_1_2",
  "hpiv_3_4",
  "rhino", "entero", "adeno", "covid"
))

levels(all_path_ccf$pathogen)

levels(all_path_ccf$pathogen) <- c(
  "Influenza A/H1N1",
  "Influenza B",
  "RSV A",
  "RSV B",
  "hMPV",
  "hCoV 229E + OC43",
  "hCoV HKU1 + NL63",
  "hPIV 1 + 2",
  "hPIV 3 + 4",
  "Rhinovirus",
  "Enterovirus",
  "Adenovirus",
  "SARS-CoV-2"
)

levels(all_path_ccf$mobility_metric)

all_path_ccf <- all_path_ccf %>%
  filter(!(mobility_metric %in% c("performing_arts_or_sports_events"))) %>%
  droplevels()
unique(all_path_ccf$mobility_metric)

all_path_ccf$mobility_metric <- factor(all_path_ccf$mobility_metric,
                                       levels = c(
                                         "within_neighborhood_movement",
                                         "within_city_movement",
                                         "within_state_movement",
                                         "out_of_state_movement",
                                         "fb_leaving_home_custom",
                                         "full_service_restaurants",
                                         "groceries_and_pharmacies",
                                         "transit",
                                         "religious_orgs", "child_day_care",
                                         "elementary_and_secondary_schools",
                                         "colleges", "not_masking",
                                         "oxford_stringency_index"
                                       )
)
levels(all_path_ccf$mobility_metric)

levels(all_path_ccf$mobility_metric) <- c(
  "Within-neighborhood\nmovement",
  "Between-neighborhood\nmovement",
  "Infux visitors\nother WA counties",
  "Influx out-of-state\nvisitors",
  "% Devices\nleaving home",
  "Restaurants",
  "Groceries and\npharmacies",
  "Transit",
  "Religious\norganizations",
  "Child daycare",
  "Elementary and\nhigh schools",
  "Colleges",
  "% not wearing masks",
  "Oxford Stringency\nIndex"
)

write_rds(all_path_ccf,file="4_Block_Bootstrap_Mobility_vs_Rt/3_make_block_bootstrap_figures/all_5mo_block_bootstrap_results.rds")
