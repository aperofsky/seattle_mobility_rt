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

source("utils.R")
########################################################################################
## import data
########################################################################################
combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()
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

com_2021 <- com_df %>% filter(date > as.Date("2021-01-01") & date < as.Date("2021-09-01"))

com_2021_df_lm <- com_2021 %>% mutate(week = isoweek(date))

make_gam_df <- function(df, pathogen, date1) {
  df2 <- df %>%
    dplyr::select(date, week, all_of(cols2), all_of(pathogen)) %>%
    filter(.data[[pathogen]] > 0) %>%
    filter(date < as.Date(date1))

  week_df <- data.frame(
    week = unique(df2$week),
    week_index = seq(1:length(unique(df2$week)))
  )

  df3 <- left_join(df2, week_df, by = "week")

  df3_scaled <-
    df3 %>%
    mutate_at(vars(within_city_movement:elementary_and_secondary_schools), ~ scale(.x) %>% as.vector())

  return(df3_scaled)
}

############################
### RSV B, Summer 2021
############################
rsv_df_lim <- make_gam_df(com_2021_df_lm, "rsv_b", "2021-08-15") # first bump in transmission
range(rsv_df_lim$date)

ggplot(rsv_df_lim) +
  geom_line(aes(x = date, y = rsv_b)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_bw()

fm <- paste("s(", names(rsv_df_lim[!(names(rsv_df_lim) %in% c("date", "week", "rsv_b"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rsv_b ~", fm))
fm

## remove gamma = 1.4 because causes convergence issues
b1 <- mgcv::gam(
  formula = fm, data = rsv_df_lim, method = "ML", select = TRUE, gamma = 1.3,
  na.action = "na.fail", family = Gamma(link = "log")
)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(rsv_df_lim %>% dplyr::select(-date, -week), "rsv_b")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
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

hmpv_df_lim <- make_gam_df(com_2021_df_lm %>% filter(date > as.Date("2021-04-15")), "hmpv", "2021-08-10")

ggplot(hmpv_df_lim) +
  geom_line(aes(x = date, y = hmpv))

fm <- paste("s(", names(hmpv_df_lim[!(names(hmpv_df_lim) %in% c("date", "week", "hmpv"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hmpv ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = hmpv_df_lim, method = "ML", gamma = 1.4, select = TRUE,
  na.action = "na.fail", family = Gamma(link = "log")
)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(hmpv_df_lim %>% dplyr::select(-date, -week), "hmpv")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod) # colleges and out of state movement
hmpv_b1 <- update(best_mod, method = "REML")

p <- draw(hmpv_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hmpv_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p1
p2 <- draw(hmpv_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(hmpv_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hmpv_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hMPV", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hmpv_p_all

############################
### scov_229E_OC43
############################

scov_229E_OC43_df_lim <- make_gam_df(com_2021_df_lm %>% filter(date > as.Date("2021-03-15")), "scov_229E_OC43", "2021-06-15")
range(scov_229E_OC43_df_lim$date)

ggplot(scov_229E_OC43_df_lim) +
  geom_line(aes(x = date, y = scov_229E_OC43)) +
  geom_hline(yintercept = 1, lty = "dashed")

fm <- paste("s(", names(scov_229E_OC43_df_lim[!(names(scov_229E_OC43_df_lim) %in% c("date", "week", "scov_229E_OC43"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_229E_OC43 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_229E_OC43_df_lim, gamma = 1.4, method = "ML", select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(scov_229E_OC43_df_lim %>% dplyr::select(-date, -week), "scov_229E_OC43")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
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

scov_HKU1_NL63_df_lim <- make_gam_df(com_2021_df_lm %>% filter(date > as.Date("2021-02-15")), "scov_HKU1_NL63", "2021-06-05")

ggplot(scov_HKU1_NL63_df_lim) +
  geom_line(aes(x = date, y = scov_HKU1_NL63)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_bw()

fm <- paste("s(", names(scov_HKU1_NL63_df_lim[!(names(scov_HKU1_NL63_df_lim) %in% c("date", "week", "scov_HKU1_NL63"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_HKU1_NL63 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_HKU1_NL63_df_lim, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(scov_HKU1_NL63_df_lim %>% dplyr::select(-date, -week), "scov_HKU1_NL63")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_HKU1_NL63_b1 <- update(best_mod, method = "REML")

p <- draw(scov_HKU1_NL63_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
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

hpiv_3_4_df_lim <- make_gam_df(com_2021_df_lm %>% filter(date > as.Date("2021-04-15")), "hpiv_3_4", "2021-07-01")

ggplot(hpiv_3_4_df_lim) +
  geom_line(aes(x = date, y = hpiv_3_4)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_bw()

fm <- paste("s(", names(hpiv_3_4_df_lim[!(names(hpiv_3_4_df_lim) %in% c("date", "week", "hpiv_3_4", "colleges"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hpiv_3_4 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = hpiv_3_4_df_lim, gamma = 1.4, select = TRUE, method = "ML", na.action = "na.fail", family = Gamma(link = "log"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(hpiv_3_4_df_lim %>% dplyr::select(-date, -week), "hpiv_3_4")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
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
all_path_2021 <- plot_grid(rsv_b_p_all, scov_229E_OC43_p_all, scov_HKU1_NL63_p_all, hmpv_p_all, hpiv_3_4_p_all,
  ncol = 2, byrow = F
)
all_path_2021

save_plot(all_path_2021, base_height = 10, base_width = 20, filename = "figures/fig_s16_partial_effects_env_virus_rebound_2021_updated.png")
