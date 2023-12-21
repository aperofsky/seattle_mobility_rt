########################################################################################
## GAMs: SARS-CoV-2 waves
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
names(combined_mob)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

## daily rt
rt_df <- read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

rt_df_wide <- rt_df %>%
  dplyr::select(date, median, organism) %>%
  distinct() %>%
  pivot_wider(names_from = "organism", values_from = "median") %>%
  dplyr::filter(date >= as.Date("2020-06-01"))

names(rt_df_wide)[2:14] <- c(
  "rsv_a", "rsv_b",
  "h1n1", "h3n2", "ivb", "rhino",
  "hpiv_3_4", "hpiv_1_2", "hmpv", "adeno", "scov_229E_OC43",
  "scov_HKU1_NL63", "covid"
)

## join rt and mobility data
combined <- left_join(rt_df_wide, combined_mob_avg, by = "date") %>%
  filter(date < as.Date("2022-06-01"))

combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

## join rt and mobility data
combined <- left_join(rt_df_wide,
  combined_mob_avg,
  by = "date"
) %>% filter(date <= as.Date("2022-05-15"))

combined[is.na(combined)] <- 0
combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined$not_masking <- 100 - combined$mask_wearing_final

## model covariates
## keep variables with univariate correlations with Rt; remove some that are redundant (e.g., remove childcare; keep schools)
cols2 <- c(
  # "within_neighborhood_movement",
  "within_city_movement",
  # "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  # "transit",
  "oxford_stringency_index",
  "not_masking"
  # "religious_orgs",
  # "child_day_care",
  # "elementary_and_secondary_schools",
  # "colleges"
)
cols2

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  complete(date = seq.Date(as.Date("2020-06-01"), as.Date("2022-05-15"), by = "day"))

###############################################################
## Filter to post SAH period
###############################################################

com_20_22 <- com_df %>% filter(date > as.Date("2020-06-01"))

com_20_22_df_lm <-
  com_20_22 %>%
  mutate(time = 1 + (date - min(date)) / ddays()) %>% # Adding numeric times
  mutate(week = isoweek(date)) %>%
  mutate(stay_at_home = if_else(date >= as.Date("2020-02-28") & date <= as.Date("2020-06-06"), 1, 0))

data.table::data.table(com_20_22_df_lm)

com_20_22_df_lm_scaled <-
  com_20_22_df_lm %>%
  mutate_at(vars(within_city_movement:colleges), ~ scale(.x) %>% as.vector())

########################################################################################
## GAMs: COVID-19 waves
########################################################################################

############################
### COVID winter 2020 - 2021
############################
covid_df_lim <- com_20_22_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), covid) %>%
  filter(covid > 0) %>%
  filter(date > as.Date("2020-08-20") & date < as.Date("2020-11-20"))

week_df <- data.frame(
  week = unique(covid_df_lim$week),
  week_index = seq(1:length(unique(covid_df_lim$week)))
)
week_df

covid_df_lim <- left_join(covid_df_lim, week_df, by = "week")
head(covid_df_lim)
names(covid_df_lim)

fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% c("date", "time", "week", "covid", "within_city_movement", "colleges"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm
b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 2, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)

gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(full_service_restaurants)") & theme_bw() & ggtitle("s(restaurants)") & xlab("restaurants")
p2
p3 <- draw(covid_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

covid_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Winter 2020-2021 Wave", theme = theme(plot.title = element_text(size = 18, face = "bold")))
covid_p_all

############################
### COVID alpha wave
############################
ggplot(com_df) +
  geom_line(aes(x = date, y = covid)) +
  geom_vline(xintercept = as.Date("2021-02-01"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-04-15"), lty = "dashed")

## approx first 3 months of wave
covid_df_lim <- com_20_22_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), covid) %>%
  filter(covid > 0) %>%
  filter(date > as.Date("2021-02-01") & date < as.Date("2021-04-15"))

week_df <- data.frame(
  week = unique(covid_df_lim$week),
  week_index = seq(1:length(unique(covid_df_lim$week)))
)
week_df

covid_df_lim <- left_join(covid_df_lim, week_df, by = "week")
head(covid_df_lim)
names(covid_df_lim)

## no statistically significant univariate correlations; remove as many covariates as possible to reduce concurvity
fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% c(
  "date", "time", "week", "covid",
  # "full_service_restaurants",
  "not_masking",
  "colleges", "oxford_stringency_index"
))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm
b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 2, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
gam.check(b1)
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(covid_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

covid_p_alpha <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Alpha Wave Spring 2021", theme = theme(plot.title = element_text(size = 18, face = "bold")))
covid_p_alpha
############################
### COVID Delta wave
############################
ggplot(com_df) +
  geom_line(aes(x = date, y = covid)) +
  geom_vline(xintercept = as.Date("2021-05-20"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-07-25"), lty = "dashed")

covid_df_lim <- com_20_22_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), covid) %>%
  filter(covid > 0) %>%
  filter(date > as.Date("2021-05-20") & date < as.Date("2021-07-25"))

week_df <- data.frame(
  week = unique(covid_df_lim$week),
  week_index = seq(1:length(unique(covid_df_lim$week)))
)
week_df

covid_df_lim <- left_join(covid_df_lim, week_df, by = "week")

fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% c("date", "time", "week", "covid", "colleges"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm
## remove gamma = 2 because causes convergence issues
b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
gam.check(b1)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
draw(best_mod, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)

covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(not_masking)") & theme_bw() & ggtitle("s(% not masking)") & xlab("% not masking")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(covid_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

covid_p_delta <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Delta Wave Summer 2021", theme = theme(plot.title = element_text(size = 18, face = "bold")))
covid_p_delta

############################
### Omicron 2021-2022
############################
ggplot(com_df) +
  geom_line(aes(x = date, y = covid)) +
  geom_vline(xintercept = as.Date("2021-10-25"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2022-01-10"), lty = "dashed")

## approx first 3 months of wave
covid_df_lim <- com_20_22_df_lm_scaled %>%
  dplyr::select(date, time, week, all_of(cols2), covid) %>%
  filter(covid > 0) %>%
  filter(date > as.Date("2021-10-25") & date < as.Date("2022-01-10"))

week_df <- data.frame(
  week = unique(covid_df_lim$week),
  week_index = seq(1:length(unique(covid_df_lim$week)))
)
week_df

covid_df_lim <- left_join(covid_df_lim, week_df, by = "week")

fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% c("date", "time", "week", "covid", "colleges", "oxford_stringency_index"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm
b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 2, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
gam.check(b1)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(not_masking)") & theme_bw() & ggtitle("s(% not masking)") & xlab("% not masking")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(full_service_restaurants)") & theme_bw() & ggtitle("s(restaurants)") & xlab("restaurants")
p2
p3 <- draw(covid_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

covid_p_omicron <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Omicron Wave 2021-2022", theme = theme(plot.title = element_text(size = 18, face = "bold")))

########################################################################################
## combine all plots into one figure
########################################################################################

covid_2020_2022 <- plot_grid(covid_p_all, covid_p_alpha, covid_p_delta, covid_p_omicron, ncol = 2, byrow = F)
covid_2020_2022
save_plot(covid_2020_2022, base_height = 8, base_width = 18, filename = "figures/fig_s11_partial_effects_waves_covid_2020_2022.png")
