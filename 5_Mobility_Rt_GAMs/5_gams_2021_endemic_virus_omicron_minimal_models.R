########################################################################################
## GAMS: Endemic virus decline during Omicron BA.1 wave, 2021
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
  "elementary_and_secondary_schools",
  "colleges"
)

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  filter(date >= as.Date("2019-01-08") & date <= as.Date("2022-05-01")) %>%
  complete(date = seq.Date(as.Date("2019-01-08"), as.Date("2022-05-01"), by = "day"))

###############################################################
######### Decline in late 2021 ##################
###############################################################

com_2021 <- com_df %>% filter(date >= as.Date("2021-11-01") & date < as.Date("2022-01-15"))

com_2021_df_lm <- com_2021 %>% mutate(week = isoweek(date))

make_gam_df <- function(df, pathogen) {
  df2 <- df %>%
    dplyr::select(date, week, all_of(cols2), all_of(pathogen)) %>%
    filter(.data[[pathogen]] > 0)

  week_df <- data.frame(
    week = unique(df2$week),
    week_index = seq(1:length(unique(df2$week)))
  )

  df3 <- left_join(df2, week_df, by = "week")

  df3_scaled <-
    df3 %>%
    mutate_at(vars(within_city_movement:colleges), ~ scale(.x) %>% as.vector())

  return(df3_scaled)
}

############################
### RSV B
############################

rsv_b_df_lim_scaled <- make_gam_df(com_2021_df_lm, "rsv_b")

fm <- paste("s(", names(rsv_b_df_lim_scaled[!(names(rsv_b_df_lim_scaled) %in% c("date", "week", "rsv_b"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rsv_b ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = rsv_b_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
b1
dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)

subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rsv_b_b1 <- update(best_mod, method = "REML")

p <- draw(rsv_b_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(rsv_b_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
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
hmpv_df_lim_scaled <- make_gam_df(com_2021_df_lm, "hmpv")

fm <- paste("s(", names(hmpv_df_lim_scaled[!(names(hmpv_df_lim_scaled) %in% c("date", "week", "hmpv"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hmpv ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = hmpv_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE,
  na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hmpv_b1 <- update(best_mod, method = "REML")

p <- draw(hmpv_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hmpv_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
p1
p2 <- draw(hmpv_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p2
p3 <- draw(hmpv_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hmpv_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hMPV", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hmpv_p_all

############################
### scov_229E_OC43 19_20
############################

scov_229E_OC43_df_lim_scaled <- make_gam_df(com_2021_df_lm, "scov_229E_OC43")

fm <- paste("s(", names(scov_229E_OC43_df_lim_scaled[!(names(scov_229E_OC43_df_lim_scaled) %in% c("date", "week", "scov_229E_OC43"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_229E_OC43 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_229E_OC43_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_229E_OC43_b1 <- update(best_mod, method = "REML")

p <- draw(scov_229E_OC43_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
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

scov_HKU1_NL63_df_lim_scaled <- make_gam_df(com_2021_df_lm, "scov_HKU1_NL63")

fm <- paste("s(", names(scov_HKU1_NL63_df_lim_scaled[!(names(scov_HKU1_NL63_df_lim_scaled) %in% c("date", "week", "scov_HKU1_NL63"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_HKU1_NL63 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_HKU1_NL63_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_HKU1_NL63_b1 <- update(best_mod, method = "REML")

p <- draw(scov_HKU1_NL63_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1

p2 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p2

p3 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

scov_HKU1_NL63_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hCoV HKU1 + NL63", theme = theme(plot.title = element_text(size = 18, face = "bold")))
scov_HKU1_NL63_p_all

############################
### hpiv 3+4
############################
hpiv_3_4_df_lim_scaled <- make_gam_df(com_2021_df_lm, "hpiv_3_4")

fm <- paste("s(", names(hpiv_3_4_df_lim_scaled[!(names(hpiv_3_4_df_lim_scaled) %in% c("date", "week", "hpiv_3_4"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hpiv_3_4 ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = hpiv_3_4_df_lim_scaled, method = "ML",
  gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hpiv_3_4_b1 <- update(best_mod, method = "REML")

p <- draw(hpiv_3_4_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hpiv_3_4_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
p2 <- draw(hpiv_3_4_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p2
p3 <- draw(hpiv_3_4_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hpiv_3_4_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hPIV 3 + 4", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hpiv_3_4_p_all

############################
### rhino
############################
rhino_df_lim_scaled <- make_gam_df(com_2021_df_lm, "rhino")

fm <- paste("s(", names(rhino_df_lim_scaled[!(names(rhino_df_lim_scaled) %in% c("date", "week", "rhino"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rhino ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = rhino_df_lim_scaled, method = "ML",
  gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rhino_b1 <- update(best_mod, method = "REML")

p <- draw(rhino_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(rhino_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
p1
p2 <- draw(rhino_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rhino_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rhino_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hRV", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rhino_p_all

############################
### adenovirus
############################
adeno_df_lim_scaled <- make_gam_df(com_2021_df_lm, "adeno")

fm <- paste("s(", names(adeno_df_lim_scaled[!(names(adeno_df_lim_scaled) %in% c("date", "week", "adeno"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("adeno ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = adeno_df_lim_scaled, method = "ML",
  gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
adeno_b1 <- update(best_mod, method = "REML")

p <- draw(adeno_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(adeno_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
p1
p2 <- draw(adeno_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(adeno_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

adeno_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "AdV", theme = theme(plot.title = element_text(size = 18, face = "bold")))

############################
### enterovirus
############################
entero_df_lim_scaled <- make_gam_df(com_2021_df_lm, "entero")

fm <- paste("s(", names(entero_df_lim_scaled[!(names(entero_df_lim_scaled) %in% c("date", "week", "entero"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("entero ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = entero_df_lim_scaled, method = "ML",
  gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)


dd <- dredge(b1,
  m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc",
  subset = !(`s(colleges)` && `s(elementary_and_secondary_schools)`)
)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
entero_b1 <- update(best_mod, method = "REML")

p <- draw(entero_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(entero_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)") & xlab("colleges")
p1
p2 <- draw(entero_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
p2
p3 <- draw(entero_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

entero_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "EV", theme = theme(plot.title = element_text(size = 18, face = "bold")))

########################################################################################
## combine all plots into one figure
########################################################################################
all_path_2021 <- plot_grid(
  rsv_b_p_all,
  hmpv_p_all,
  scov_229E_OC43_p_all,
  scov_HKU1_NL63_p_all,
  hpiv_3_4_p_all,
  rhino_p_all,
  entero_p_all,
  adeno_p_all,
  ncol = 2, byrow = F
)
all_path_2021
save_plot(all_path_2021, base_height = 12, base_width = 18, filename = "figures/fig_s18_partial_effects_endemic_pathogens_omicron_2021_2022_updated.png")
