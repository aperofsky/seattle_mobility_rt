library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(cdcfluview)
library(lubridate)
library(tidymodels)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

########################################################################################
## Import Nextstrain sequence data
########################################################################################

nextstrain_lineage <- fread("6_LASSO_Predictive_Models/3_VOC_data/gisaid_metadata.tsv.gz")
head(nextstrain_lineage)
sort(names(nextstrain_lineage))
sort(unique(nextstrain_lineage$country))
sort(unique(nextstrain_lineage$Nextstrain_clade))

## Clean data
## filter to US sequences with good QC scores
nextstrain_lineage_us <- nextstrain_lineage %>%
  filter(country == "USA") %>%
  filter(QC_overall_status == "good")
rm(nextstrain_lineage)

## filter to King County WA sequences
KC_variants <- nextstrain_lineage_us %>%
  filter(location == "King County WA") %>%
  dplyr::select(
    strain, date, region, country, division, location,
    region_exposure, country_exposure, division_exposure,
    host, Nextstrain_clade, pango_lineage, clade_nextstrain, clade_who,
    Nextclade_pango, originating_lab, submitting_lab
  )
rm(nextstrain_lineage_us)

KC_variants$date_char <- nchar(KC_variants$date)
KC_variants <- KC_variants %>%
  filter(date_char == 10) %>%
  dplyr::select(-date_char)
KC_variants$date <- as.Date(KC_variants$date)
range(KC_variants$date) # "2020-02-21" "2024-02-23"

KC_variants_2020_2022 <- KC_variants %>% filter(date < as.Date("2022-07-01"))
names(KC_variants_2020_2022)
sort(unique(KC_variants_2020_2022$Nextstrain_clade))
write_csv(KC_variants_2020_2022, file = "6_LASSO_Predictive_Models/3_VOC_data/king_county_SC2_clades_2020_2022.csv")

range(KC_variants_2020_2022$date) #"2020-02-21" "2022-06-30"
nrow(KC_variants_2020_2022)#9626
length(unique(KC_variants_2020_2022$strain))#9626

## rename clades
KC_variants_2020_2022 <- KC_variants_2020_2022 %>% mutate(clade_new = case_when(
  Nextstrain_clade == "20G" ~ "20G (B.1.2)",
  Nextstrain_clade == "20A" ~ "20A (B.1)",
  Nextstrain_clade == "20C" ~ "20C",
  Nextstrain_clade == "20I (Alpha, V1)" ~ "Alpha",
  Nextstrain_clade == "20H (Beta, V2)" ~ "Beta",
  Nextstrain_clade == "21J (Delta)" ~ "Delta",
  Nextstrain_clade == "21I (Delta)" ~ "Delta",
  Nextstrain_clade == "21A (Delta)" ~ "Delta",
  Nextstrain_clade == "21C (Epsilon)" ~ "Epsilon",
  Nextstrain_clade == "19B" ~ "Ancestral",
  Nextstrain_clade == "19A" ~ "Ancestral",
  Nextstrain_clade == "20B" ~ "20B (B.1.1)",
  Nextstrain_clade == "20D" ~ "20B (B.1.1)",
  Nextstrain_clade == "20G" ~ "20G (B.1.2)",
  Nextstrain_clade == "20J (Gamma, V3)" ~ "Gamma",
  Nextstrain_clade == "21K (Omicron)" ~ "Omicron BA.1",
  Nextstrain_clade == "21L (Omicron)" ~ "Omicron BA.2",
  Nextstrain_clade == "22C (Omicron)" ~ "Omicron BA.2",
  Nextstrain_clade == "22A (Omicron)" ~ "Omicron BA.4",
  Nextstrain_clade == "22B (Omicron)" ~ "Omicron BA.5",
  TRUE ~ as.character(Nextstrain_clade)
))

variant_table <- KC_variants_2020_2022 %>% distinct(Nextstrain_clade, clade_nextstrain, clade_who, clade_new)
variant_table %>% arrange(clade_new)

range(KC_variants_2020_2022$date) # "2020-02-21" "2022-06-30"

## fill in dataframe so each clade has an entry for each date
filled_df <- KC_variants_2020_2022 %>%
  group_by(date, clade_new) %>%
  tally() %>%
  ungroup() %>%
  group_by(clade_new) %>%
  complete(date = seq.Date(as.Date("2020-02-21"), as.Date("2022-06-30"), by = "day")) %>%
  replace_na(list(n = 0)) %>%
  group_by(clade_new) %>%
  mutate(cum_count = cumsum(n)) %>%
  ungroup()

filled_df %>%
  group_by(clade_new) %>%
  summarize(sum = sum(n)) %>%
  arrange(-sum) %>%
  print(n = 50)

## assign low count clades to "other"
low_count_clades <- filled_df %>%
  group_by(clade_new) %>%
  summarize(sum = sum(n)) %>%
  arrange(-sum) %>%
  filter(sum < 50) %>%
  pull(clade_new)
low_count_clades

ggplot(filled_df) +
  geom_col(aes(x = date, y = n, fill = clade_new))

ggplot(filled_df) +
  geom_line(aes(x = date, y = cum_count, color = clade_new))

## when do cum counts exceed 20?
filled_df %>%
  group_by(clade_new) %>%
  filter(cum_count >= 20) %>%
  slice_head(n = 1) %>%
  arrange(date)

filled_df <- bind_cols(filled_df, mmwr_week(filled_df$date)[, 1:2])
filled_df$epi_week <- mmwr_week_to_date(filled_df$mmwr_year, filled_df$mmwr_week)

filled_df <- filled_df %>%
  group_by(clade_new, epi_week) %>%
  summarize(n = sum(n)) %>%
  ungroup()

filled_df %>%
  filter(clade_new == "Alpha") %>%
  filter(n > 0) %>%
  print(n = 20)

filled_df %>%
  filter(clade_new == "Delta") %>%
  filter(n > 0) %>%
  print(n = 20)

filled_df %>%
  filter(clade_new == "Omicron BA.1") %>%
  filter(n > 0) %>%
  print(n = 20)

filled_df %>%
  filter(clade_new == "Omicron BA.2") %>%
  filter(n > 0) %>%
  print(n = 20)

filled_df %>%
  filter(clade_new == "Omicron BA.4") %>%
  filter(n > 0) %>%
  print(n = 20)

KC_variants_wide <- KC_variants_2020_2022 %>%
  group_by(date, location, clade_new) %>%
  tally() %>%
  pivot_wider(names_from = clade_new, values_from = n) %>%
  ungroup()

KC_variants_wide[is.na(KC_variants_wide)] <- 0
data.table(KC_variants_wide)

KC_variants_wide <-
  KC_variants_wide %>%
  rowwise() %>%
  mutate(total = sum(c_across(Ancestral:`Omicron BA.4`)))

KC_variants_wide <-
  KC_variants_wide %>%
  mutate_at(vars(Ancestral:`Omicron BA.4`), ~ . / total)

names(KC_variants_wide)

data.table(KC_variants_wide)

KC_variants_long <- KC_variants_wide %>%
  pivot_longer(cols = c(Ancestral:`Omicron BA.4`), names_to = "clade_new", values_to = "prop")

KC_variants_long$date_char <- nchar(KC_variants_long$date)
KC_variants_long$date <- as.Date(KC_variants_long$date)
range(KC_variants_long$date)

########################################################################################
## Multinomial logistic regression model
########################################################################################

KC_variants_lm <-
  KC_variants_2020_2022 %>%
  mutate(
    times = 1 + (date - min(date)) / ddays(), # Adding numeric times
    clade_new = as.factor(clade_new)
  )
KC_variants_lm %>% arrange(clade_new, date)

# Setting workflow for tidymodels
multinomial_recipe <- recipe(clade_new ~ times, data = KC_variants_lm) %>%
  step_normalize(all_predictors()) # normalize standard dev and mean

multinomial_prep <- prep(multinomial_recipe)

multinomial_spec <- multinom_reg() %>%
  set_engine("nnet") # setting nnet package for multinom reg fit

multinomial_wf <- workflow() %>%
  add_recipe(multinomial_recipe) %>%
  add_model(multinomial_spec)

model_results <- fit(multinomial_wf, KC_variants_lm)

times <- seq(min(KC_variants_lm$times), max(KC_variants_lm$times))
date <- times + min(KC_variants_lm$date) - 1

times_to_predict <- data.frame(times = times, date = date)

pred_df <- cbind(times_to_predict, predict(model_results, times_to_predict, type = "prob")) %>%
  distinct(date, .keep_all = TRUE) %>%
  pivot_longer(starts_with(".pred"), names_to = "clade_new", values_to = "pred_freq") %>%
  mutate(clade_new = str_replace(clade_new, "^.pred_", ""))
sort(unique(pred_df$clade_new))

pred_df <- pred_df %>%
  mutate(clade_new = if_else(clade_new %in% low_count_clades, "Other", clade_new)) %>%
  group_by(times, date, clade_new) %>%
  summarize(pred_freq = sum(pred_freq)) %>%
  droplevels()
length(unique(pred_df$clade_new)) # 14
unique(pred_df$clade_new)

pred_df$clade_new <- factor(pred_df$clade_new, levels = c(
  "Ancestral", "20A (B.1)", "20B (B.1.1)", "20G (B.1.2)", "20C", "Epsilon",
  "Alpha", "Other", "Gamma", "Delta", "Omicron BA.1", "Omicron BA.2",
  "Omicron BA.4", "Omicron BA.5"
))

write_rds(pred_df, file = "6_LASSO_Predictive_Models/3_VOC_data/KC_WA_covid_cases_by_variant_2020_2022.rds")

nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols)

voc_freq <- ggplot(pred_df, aes(x = date, y = pred_freq, fill = clade_new)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(expand = c(0, 0), date_labels = "%b %Y", date_breaks = "4 months") +
  ylab("Predicted Frequency") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = mycolors, name = "Clade") +
  xlab("Date")
voc_freq
save_plot(voc_freq, filename = "figures/fig_s22_sc2_variant_frequency.png", base_width = 14, base_height = 7)
