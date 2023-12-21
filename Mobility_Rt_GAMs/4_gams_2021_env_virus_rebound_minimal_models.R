########################################################################################
## GAMS: enveloped virus rebound, 2021
########################################################################################
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(lubridate)
library(mgcv)
library(MuMIn)
library(gratia)
library(patchwork)

########################################################################################
## import data
########################################################################################

## mobility data
combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

## daily rt
rt_df <- read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

rt_df_wide <- rt_df %>%
  dplyr::select(date, median, organism) %>%
  distinct() %>%
  pivot_wider(names_from = "organism", values_from = "median") %>%
  dplyr::filter(date > as.Date("2018-11-20"))

names(rt_df_wide)[2:14] <- c(
  "rsv_a", "rsv_b",
  "h1n1", "h3n2", "ivb", "rhino",
  "hpiv_3_4", "hpiv_1_2", "hmpv", "adeno", "scov_229E_OC43",
  "scov_HKU1_NL63", "covid"
)
combined <- left_join(rt_df_wide,
  combined_mob_avg,
  by = "date"
) %>% filter(date < as.Date("2022-06-01"))

combined[is.na(combined)] <- 0
combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined$not_masking <- 100 - combined$mask_wearing_final

cols2 <- c(
  # "within_neighborhood_movement",
  "within_city_movement",
  # "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  # "full_service_restaurants",
  "religious_orgs",
  # "child_day_care",
  # "not_masking",
  # "oxford_stringency_index",
  "elementary_and_secondary_schools"
  # "colleges"
)

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  filter(date >= as.Date("2019-01-08") & date <= as.Date("2022-05-01")) %>%
  complete(date = seq.Date(as.Date("2019-01-08"), as.Date("2022-05-01"), by = "day"))

###############################################################
######### Jan 2021 to August 2021 ##################
###############################################################

com_2021 <- com_df %>% filter(date > as.Date("2021-01-15") & date < as.Date("2021-08-15"))

com_2021_df_lm <-
  com_2021 %>%
  mutate(time = 1 + (date - min(date)) / ddays()) %>% # Adding numeric times
  mutate(week = isoweek(date)) %>%
  mutate(stay_at_home = if_else(date >= as.Date("2020-02-28") & date <= as.Date("2020-06-06"), 1, 0))

names(com_2021_df_lm)

com_2021_df_lm_scaled <-
  com_2021_df_lm %>%
  mutate_at(vars(within_city_movement:colleges), ~ scale(.x) %>% as.vector())

############################
### RSV B, Summer 2021
############################
rsv_df_lim <- com_2021_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), rsv_b) %>%
  filter(rsv_b > 0)

week_df <- data.frame(
  week = unique(rsv_df_lim$week),
  week_index = seq(1:length(unique(rsv_df_lim$week)))
)
week_df

rsv_df_lim <- left_join(rsv_df_lim, week_df, by = "week")

fm <- paste("s(", names(rsv_df_lim[!(names(rsv_df_lim) %in% c("date", "time", "week", "rsv_b", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rsv_b ~", fm))
fm
b1 <- mgcv::gam(formula = fm, data = rsv_df_lim, method = "ML", select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rsv_b_b1 <- update(best_mod, method = "REML")

p <- draw(rsv_b_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(rsv_b_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p1
p2 <- draw(rsv_b_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rsv_b_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rsv_b_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "RSV B", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rsv_b_p_all

############################
### HMPV
############################
hmpv_df_lim <- com_2021_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), hmpv) %>%
  filter(hmpv > 0)
range(hmpv_df_lim$date) # "2021-04-06" "2021-07-05"

week_df <- data.frame(
  week = unique(hmpv_df_lim$week),
  week_index = seq(1:length(unique(hmpv_df_lim$week)))
)
week_df

hmpv_df_lim <- left_join(hmpv_df_lim, week_df, by = "week") %>% filter(date < as.Date("2021-07-06"))

fm <- paste("s(", names(hmpv_df_lim[!(names(hmpv_df_lim) %in% c("date", "time", "week", "hmpv", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hmpv ~", fm))
fm
b1 <- mgcv::gam(
  formula = fm, data = hmpv_df_lim, method = "ML", gamma = 2, select = TRUE,
  na.action = "na.fail", family = Gamma(link = "log")
)
gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod) # colleges and out of state movement
hmpv_b1 <- update(best_mod, method = "REML")

p <- draw(hmpv_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hmpv_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p1
p2 <- draw(hmpv_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(hmpv_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")

hmpv_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hMPV", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hmpv_p_all

############################
### scov_229E_OC43
############################

scov_229E_OC43_df_lim <- com_2021_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), scov_229E_OC43) %>%
  filter(scov_229E_OC43 > 0) %>%
  filter(date < as.Date("2021-06-01"))
range(scov_229E_OC43_df_lim$date)

week_df <- data.frame(
  week = unique(scov_229E_OC43_df_lim$week),
  week_index = seq(1:length(unique(scov_229E_OC43_df_lim$week)))
)
week_df

scov_229E_OC43_df_lim <- left_join(scov_229E_OC43_df_lim, week_df, by = "week")

fm <- paste("s(", names(scov_229E_OC43_df_lim[!(names(scov_229E_OC43_df_lim) %in% c("date", "time", "week", "scov_229E_OC43", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_229E_OC43 ~", fm))
fm
b1 <- mgcv::gam(formula = fm, data = scov_229E_OC43_df_lim, gamma = 2, method = "ML", select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))
gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_229E_OC43_b1 <- update(best_mod, method = "REML")

p <- draw(scov_229E_OC43_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p1

p2 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2

p3 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

scov_229E_OC43_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hCoV 229E + OC43", theme = theme(plot.title = element_text(size = 18, face = "bold")))
scov_229E_OC43_p_all

############################
### scov HKU1 NL63
############################
scov_HKU1_NL63_df_lim <- com_2021_df_lm %>%
  dplyr::select(date, time, week, all_of(cols2), scov_HKU1_NL63) %>%
  filter(scov_HKU1_NL63 > 0) %>%
  filter(date < as.Date("2021-04-21"))
range(scov_HKU1_NL63_df_lim$date)

week_df <- data.frame(
  week = unique(scov_HKU1_NL63_df_lim$week),
  week_index = seq(1:length(unique(scov_HKU1_NL63_df_lim$week)))
)
week_df

scov_HKU1_NL63_df_lim <- left_join(scov_HKU1_NL63_df_lim, week_df, by = "week")

fm <- paste("s(", names(scov_HKU1_NL63_df_lim[!(names(scov_HKU1_NL63_df_lim) %in% c("date", "time", "week", "scov_HKU1_NL63", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_HKU1_NL63 ~", fm))
fm
b1 <- mgcv::gam(formula = fm, data = scov_HKU1_NL63_df_lim, method = "ML", gamma = 2, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))
gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_HKU1_NL63_b1 <- update(best_mod, method = "REML")

p <- draw(scov_HKU1_NL63_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
p1
p2 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

scov_HKU1_NL63_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hCoV HKU1 + NL63", theme = theme(plot.title = element_text(size = 18, face = "bold")))
scov_HKU1_NL63_p_all

############################
### hpiv 3+4
############################
hpiv_3_4_df_lim <- com_2021_df_lm %>%
  dplyr::select(date, time, week, all_of(cols2), hpiv_3_4) %>%
  filter(hpiv_3_4 > 0) %>%
  filter(date < as.Date("2021-06-10"))
range(hpiv_3_4_df_lim$date)

week_df <- data.frame(
  week = unique(hpiv_3_4_df_lim$week),
  week_index = seq(1:length(unique(hpiv_3_4_df_lim$week)))
)
week_df

hpiv_3_4_df_lim <- left_join(hpiv_3_4_df_lim, week_df, by = "week")
range(hpiv_3_4_df_lim$date)

fm <- paste("s(", names(hpiv_3_4_df_lim[!(names(hpiv_3_4_df_lim) %in% c("date", "time", "week", "hpiv_3_4", "colleges", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hpiv_3_4 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = hpiv_3_4_df_lim, gamma = 2, select = TRUE, method = "ML", na.action = "na.fail", family = Gamma(link = "log"))
gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hpiv_3_4_b1 <- update(best_mod, method = "REML")

p <- draw(hpiv_3_4_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hpiv_3_4_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p1
p2 <- draw(hpiv_3_4_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(hpiv_3_4_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hpiv_3_4_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hPIV 3 + 4", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hpiv_3_4_p_all

########################################################################################
## combine all plots into one figure
########################################################################################
all_path_2021 <- plot_grid(rsv_b_p_all, scov_229E_OC43_p_all, scov_HKU1_NL63_p_all, hmpv_p_all,
  hpiv_3_4_p_all,
  ncol = 2, byrow = F
)
save_plot(all_path_2021, base_height = 10, base_width = 20, filename = "figures/fig_s15_partial_effects_env_virus_rebound_2021.png")
