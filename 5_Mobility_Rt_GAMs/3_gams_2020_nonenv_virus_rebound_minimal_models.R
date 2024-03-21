########################################################################################
## GAMS: rhinovirus and adenovirus rebound, June 2020 - Dec 2020
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
  "full_service_restaurants",
  "religious_orgs"
  # "child_day_care",
  # "not_masking",
  # "oxford_stringency_index",
  # "elementary_and_secondary_schools"
  # "colleges"
)

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  filter(date >= as.Date("2019-01-08") & date <= as.Date("2022-05-01")) %>%
  complete(date = seq.Date(as.Date("2019-01-08"), as.Date("2022-05-01"), by = "day"))

###############################################################
######### June to December 2020 ##################
###############################################################
## first 6 months of rebound
com_2020 <- com_df %>% filter(date >= as.Date("2020-06-01") & date < as.Date("2020-12-01"))

com_2020_df_lm <- com_2020 %>% mutate(week = isoweek(date))

make_gam_df <- function(df, pathogen, date1, date2) {
  df2 <- df %>%
    dplyr::select(date, week, all_of(cols2), all_of(pathogen)) %>%
    filter(.data[[pathogen]] > 0) %>%
    filter(date > as.Date(date1) & date < as.Date(date2))

  week_df <- data.frame(
    week = unique(df2$week),
    week_index = seq(1:length(unique(df2$week)))
  )

  df3 <- left_join(df2, week_df, by = "week")

  df3_scaled <-
    df3 %>%
    mutate_at(vars(within_city_movement:religious_orgs), ~ scale(.x) %>% as.vector())

  return(df3_scaled)
}

############################
### rhinovirus
############################

rhino_df_lim <- make_gam_df(com_2020_df_lm, "rhino", date1 = min(com_2020_df_lm$date), date2 = "2020-11-01")
names(rhino_df_lim)

ggplot(rhino_df_lim) +
  geom_line(aes(x = date, y = rhino))

fm <- paste("s(", names(rhino_df_lim[!(names(rhino_df_lim) %in% c("date", "rhino", "week"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rhino ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = rhino_df_lim, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

smat <- subset_cor(rhino_df_lim %>% dplyr::select(-date, -week, -rhino), "rhino")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rhino_b1 <- update(best_mod, method = "REML")
summary(rhino_b1)

p <- draw(rhino_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(rhino_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p1
p2 <- draw(rhino_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rhino_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rhino_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Rhinovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rhino_p_all

############################
### adenovirus
############################
adeno_df_lim <- make_gam_df(com_2020_df_lm, "adeno", date1 = min(com_2020_df_lm$date), date2 = "2020-11-15")

ggplot(adeno_df_lim) +
  geom_line(aes(x = date, y = adeno))

fm <- paste("s(", names(adeno_df_lim[!(names(adeno_df_lim) %in% c("date", "week", "adeno"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("adeno ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = adeno_df_lim, method = "ML", gamma = 1.4, select = TRUE, family = Gamma(link = "log"), na.action = "na.fail")

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
adeno_b1 <- update(best_mod, method = "REML")
summary(adeno_b1)

p <- draw(adeno_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(adeno_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p1
p2 <- draw(adeno_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(adeno_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")

adeno_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Adenovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
adeno_p_all

############################
### enterovirus
############################
entero_df_lim <- make_gam_df(com_2020_df_lm, "entero", date1 = "2020-06-15", date2 = "2020-11-15")

ggplot(entero_df_lim) +
  geom_line(aes(x = date, y = entero)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  theme_bw()

fm <- paste("s(", names(entero_df_lim[!(names(entero_df_lim) %in% c("date", "week", "entero"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("entero ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = entero_df_lim, method = "ML", gamma = 1.4, select = TRUE, family = Gamma(link = "log"), na.action = "na.fail")

smat <- subset_cor(entero_df_lim, "entero")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
entero_b1 <- update(best_mod, method = "REML")
summary(entero_b1)

p <- draw(entero_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(entero_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p1
p2 <- draw(entero_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(entero_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")

entero_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Enterovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
entero_p_all

########################################################################################
## combine all plots into one figure
########################################################################################
non_env_p_2020 <- plot_grid(rhino_p_all, entero_p_all, adeno_p_all, nrow = 3)
non_env_p_2020
save_plot(non_env_p_2020, base_height = 12, base_width = 16, filename = "figures/fig_s14_partial_effects_non_env_viruses_june_to_dec_2020_updated.png")
