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
library(scales)

source("utils.R")

########################################################################################
## import data
########################################################################################
combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()
combined[is.na(combined)] <- 0
combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined$not_masking <- 100 - combined$mask_wearing_final

## model covariates
## keep variables with univariate correlations with Rt; remove some that are redundant (e.g., remove childcare; keep schools)
cols2 <- c(
  "within_neighborhood_movement",
  "within_city_movement",
  # "within_state_movement", #similar to out of state
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  # "groceries_and_pharmacies",
  "transit",
  "oxford_stringency_index",
  "not_masking"
  # "religious_orgs"
  # "child_day_care",
  # "elementary_and_secondary_schools",
  # "colleges"
)
cols2

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  complete(date = seq.Date(as.Date("2020-06-01"), as.Date("2022-05-15"), by = "day"))

voc_cases = read_rds("6_LASSO_Predictive_Models/3_VOC_data/KC_WA_covid_cases_by_variant_2020_2022.rds")

voc_wide = voc_cases %>% dplyr::select(date,clade_new,pred_freq) %>% pivot_wider(names_from = clade_new, values_from = pred_freq)

ggplot(voc_wide %>%  
         dplyr::select(date,Ancestral,contains("B.1"),"20C") %>% 
         mutate(sum = rowSums(across(where(is.numeric)))))+
  geom_line(aes(x=date,y=sum))

###############################################################
## Filter to post SAH period
###############################################################

com_20_22 <- com_df %>% filter(date > as.Date("2020-06-01"))

com_20_22_df_lm <- com_20_22 %>% mutate(week = isoweek(date))

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
    mutate_at(vars(within_neighborhood_movement:not_masking), ~ scale(.x) %>% as.vector())

  return(df3_scaled)
}

########################################################################################
## GAMs: COVID-19 waves
########################################################################################

############################
### COVID winter 2020 - 2021
############################
## no vaccination yet
covid_df_lim <- make_gam_df(com_20_22_df_lm, "covid", "2020-08-20", "2020-11-20")
head(covid_df_lim)

ggplot(covid_df_lim)+
  geom_line(aes(x=date,y=covid))

ggplot(covid_df_lim)+
  geom_line(aes(x=date,y=not_masking))

voc_winter_2020_2021 <- voc_wide %>% filter(date>as.Date("2020-08-20") & date < as.Date("2020-11-20"))

## still mostly ancestral strain (98%)
ggplot(voc_winter_2020_2021 %>%  
         dplyr::select(date,Ancestral,contains("B.1"),"20C") %>% 
         mutate(sum = rowSums(across(where(is.numeric)))))+
  geom_line(aes(x=date,y=sum))

rem = c("date","week", "covid",
        "oxford_stringency_index",
        "within_neighborhood_movement",
        "within_city_movement",
        # "groceries_and_pharmacies",
        "transit")
fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% rem)]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)

#limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = covid_df_lim %>% dplyr::select(-all_of(rem)), pathogen = "covid")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset=smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(full_service_restaurants)") & theme_bw() & ggtitle("s(restaurants)") & xlab("restaurants")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
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
  geom_vline(xintercept = as.Date("2021-02-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-04-20"), lty = "dashed")

ggplot(com_df%>% filter(date>as.Date("2021-02-15") & date < as.Date("2021-04-20"))) +
  geom_line(aes(x = date, y = covid)) 

## approx first 3 months of wave
covid_df_lim <- make_gam_df(com_20_22_df_lm, "covid", "2021-01-20", "2021-04-15")
range(covid_df_lim$covid)
unique(covid_df_lim$week_index)
names(covid_df_lim)

voc_alpha <- voc_wide %>% filter(date>as.Date("2021-02-15") & date < as.Date("2021-04-20"))

## no vaccination yet
ggplot(voc_alpha %>%  
         dplyr::select(date,Ancestral,contains("B.1"),"20C") %>% 
         mutate(sum = rowSums(across(where(is.numeric)))))+
  geom_line(aes(x=date,y=sum))

voc_alpha = voc_alpha %>%
  dplyr::select(date,Ancestral,contains("B.1"),"20C","Alpha") %>% 
  mutate(ancestralsum = rowSums(across(c(Ancestral,contains("B.1"),"20C"))))%>%
  mutate(per_voc = 1 - ancestralsum)%>%
  mutate(per_alpha = Alpha)

ggplot(voc_alpha)+
  geom_line(aes(x=date,y=Alpha),color="blue")+
  geom_line(aes(x=date,y=per_voc))

rem = c("date","week", "covid",
        # "within_city_movement",
        # "fb_leaving_home_custom",
        "not_masking",
        "transit",
        "within_neighborhood_movement",
        # "religious_orgs",
        "oxford_stringency_index")

## no statistically significant univariate correlations; remove as many covariates as possible to reduce concurvity
fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% rem)]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)

#limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = covid_df_lim %>% dplyr::select(-all_of(rem)), pathogen = "covid")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")

p <- draw(covid_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(full_service_restaurants)") & theme_bw() & ggtitle("s(restaurants)") & xlab("restaurants")
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
  geom_vline(xintercept = as.Date("2021-05-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2021-08-01"), lty = "dashed")

covid_df_lim <- make_gam_df(com_20_22_df_lm, "covid", "2021-05-15", "2021-07-25")

voc_delta <- voc_wide %>% filter(date>as.Date("2021-05-15") & date < as.Date("2021-07-25"))

voc_delta = voc_delta %>%
  dplyr::select(date,Ancestral,contains("B.1"),"20C",Delta) %>% 
  mutate(ancestralsum = rowSums(across(c(Ancestral,contains("B.1"),"20C"))))%>%
  mutate(per_voc = 1 - ancestralsum)%>%
  mutate(per_delta = Delta)

ggplot(voc_delta)+
  geom_line(aes(x=date,y=per_delta),color="blue")+
  geom_line(aes(x=date,y=per_voc))

rem = c("date","week", "covid",
        "within_neighborhood_movement",
        "within_city_movement",
        "transit",
        # "religious_orgs",
        "oxford_stringency_index")

fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% rem)]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm

## remove gamma = 1.4 because causes convergence issues
b1 <- mgcv::gam(
  formula = fm, 
  data = covid_df_lim, gamma = 1.1,
  method = "ML", select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)

#limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = covid_df_lim %>% dplyr::select(-all_of(rem)), pathogen = "covid")
smat

## everything is collinear so remove subset 
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
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
ggplot(com_df %>% filter(date>as.Date("2021-09-01"))) +
  geom_line(aes(x = date, y = covid)) +
  geom_vline(xintercept = as.Date("2021-10-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2022-01-05"), lty = "dashed")

## approx first 3 months of wave
covid_df_lim <- make_gam_df(com_20_22_df_lm, "covid", "2021-10-15", "2022-01-05")

voc_omicron <- voc_wide %>% filter(date>as.Date("2021-10-15") & date < as.Date("2022-01-05"))

voc_omicron = voc_omicron %>% 
  dplyr::select(date,`Omicron BA.1`,`Omicron BA.2`) %>% 
  mutate(per_omicron = `Omicron BA.1`)

ggplot(voc_omicron)+
  geom_line(aes(x=date,y=per_omicron),color="blue")+
  geom_vline(aes(xintercept=as.Date("2021-12-01")),lty="dashed")+
  theme_bw()

rem = c("date","week",
        "within_city_movement",
        "transit",
        "not_masking",
        "oxford_stringency_index")

fm <- paste("s(", names(covid_df_lim[!(names(covid_df_lim) %in% c("covid",rem))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("covid ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = covid_df_lim, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
summary(b1)
gam.check(b1)

#limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = covid_df_lim %>% dplyr::select(-all_of(rem)), pathogen = "covid")
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank="AICc",subset=smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
covid_b1 <- update(best_mod, method = "REML")
gam.check(covid_b1)

p <- draw(covid_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(covid_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1
p2 <- draw(covid_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(covid_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

covid_p_omicron <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Omicron Wave 2021-2022", theme = theme(plot.title = element_text(size = 18, face = "bold")))
covid_p_omicron

########################################################################################
## combine all plots into one figure
########################################################################################
covid_2020_2022 <- plot_grid(covid_p_all, covid_p_alpha, covid_p_delta, covid_p_omicron, ncol = 2, byrow = F)
covid_2020_2022
save_plot(covid_2020_2022, base_height = 8, base_width = 18, filename = "figures/fig_s12_partial_effects_waves_covid_2020_2022_updated.png")
