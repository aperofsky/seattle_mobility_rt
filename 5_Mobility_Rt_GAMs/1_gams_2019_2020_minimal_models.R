########################################################################################
## GAMS: 2019-2020 season (pre-pandemic)
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
library(performance)
library(scales)

source("utils.R")

########################################################################################
## import data
########################################################################################
combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()
combined[is.na(combined)] <- 0
combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom

## model covariates
cols2 <- c(
  # "within_neighborhood_movement",
  "within_city_movement",
  # "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  # "full_service_restaurants",
  # "transit",
  "religious_orgs",
  # "child_day_care",
  "elementary_and_secondary_schools",
  "colleges"
)

com_df <- combined %>%
  dplyr::select(date, all_of(cols2), rsv_a:covid) %>%
  filter(date >= as.Date("2019-01-08") & date <= as.Date("2022-05-01")) %>%
  complete(date = seq.Date(as.Date("2019-01-08"), as.Date("2022-05-01"), by = "day"))

########################################################################################
# make input dataframe for GAMs
########################################################################################

make_gam_df <- function(df, pathogen) {
  
  first_date <- df %>%
    dplyr::select(date, week, all_of(cols2), all_of(pathogen)) %>%
    filter(.data[[pathogen]] > 0) %>%
    slice_min(date) %>%
    pull(date)
  
  df2 <- df %>%
    dplyr::select(date, week, all_of(cols2), all_of(pathogen)) %>%
    filter(date >= as.Date(first_date))%>%
    filter(.data[[pathogen]] > 0)
  
  week_df <- data.frame(
    week = unique(df2$week),
    week_index = seq(1:length(unique(df2$week)))
  )
  
  df3 <- left_join(df2, week_df, by = "week")
  
  df4_scaled <-
    df3 %>%
    mutate_at(vars(within_city_movement:colleges), ~ scale(.x) %>% as.vector())
  
  return(df4_scaled)
}

###############################################################
## Filter to pre-pandemic period
###############################################################

com_19_20 <- com_df %>% filter(date > as.Date("2019-09-01") & date <= as.Date("2020-03-01"))

com_19_20_df_lm <- com_19_20 %>% mutate(week = isoweek(date))
head(com_19_20_df_lm)

###############################################################
## view incidence data
###############################################################

## incidence data
inc_df <- read_rds("2_Epidemia_Models/Epidemia_Models_Biowulf/extended_daily_ILI_per_pos_scaled_for_select_pathogens_2023_10_22.rds")
head(inc_df)

# 3-week moving average to reduce noise but not oversmooth
inc_df_sm <- inc_df %>%
  rename(inc = all_sites_per_pos_ili_scaled_sum) %>%
  group_by(organism) %>%
  mutate(inc = round(zoo::rollmean(x = inc, k = 21, fill = F, align = "center") * 1000)) %>%
  mutate(inc = zoo::na.approx(inc, maxgap = 3, na.rm = F)) %>%
  mutate(inc = tidyr::replace_na(inc, 0)) %>%
  dplyr::select(organism, date, inc)

head(inc_df_sm)
unique(inc_df_sm$organism)

inc_df_sm <- inc_df_sm %>%
  filter(organism %in% c(
    "RSV A", "RSV B", "IVA/H1N1",
    # "IVA/H3N2",
    "IVB", "Rhinovirus",
    "HPIV-3 + HPIV-4", "HPIV-1 + HPIV-2", "hMPV", "Adenovirus", "EV",
    "CoV-229E + CoV-OC43", "CoV-HKU1 + CoV-NL63"
  )) %>%
  filter(date > as.Date("2019-09-01") & date <= as.Date("2020-03-01")) %>%
  group_by(organism) %>%
  mutate(cum_inc = cumsum(inc)) %>%
  dplyr::select(organism, date, inc, cum_inc) %>%
  droplevels()

head(inc_df_sm)
names(com_19_20_df_lm)
levels(inc_df_sm$organism)

levels(inc_df_sm$organism) <- c(
  "adeno", "scov_229E_OC43", "scov_HKU1_NL63", "entero",
  "hmpv", "hpiv_1_2", "hpiv_3_4", "h1n1",
  "ivb", "rsv_a", "rsv_b", "rhino"
)
head(inc_df_sm)

ggplot(inc_df_sm) +
  geom_line(aes(x = date, y = cum_inc)) +
  scale_y_continuous(trans = log_trans()) +
  facet_wrap(~organism, scales = "free")

ggplot(inc_df_sm) +
  geom_line(aes(x = date, y = inc)) +
  facet_wrap(~organism, scales = "free") +
  geom_vline(xintercept = as.Date("2019-10-01"), lty = "dashed") +
  # geom_vline(xintercept = as.Date("2020-01-15"), lty = "dashed") +
  theme_bw()

first_peak_df <- inc_df_sm %>%
  group_by(organism) %>%
  slice_max(inc, n = 1) %>%
  rename(peak_date = date) %>%
  arrange(peak_date) %>%
  ungroup()
first_peak_df

first_peak_df %>%
  group_by(organism) %>%
  slice_max(peak_date, n = 1) %>%
  arrange(organism, peak_date) %>%
  print(n = 30)
# organism       peak_date    inc cum_inc
# <fct>          <date>     <dbl>   <dbl>
# 1 adeno          2020-03-01   263   17792
# 2 scov_229E_OC43 2020-02-07    73    3325
# 3 scov_HKU1_NL63 2020-01-01   559   16564
# 4 entero         2019-11-29   167    2734
# 5 hmpv           2020-02-18   299   19531
# 6 hpiv_1_2       2019-11-16   356   11422
# 7 hpiv_3_4       2019-11-29    66    1218
# 8 h1n1           2020-01-28   343   10448
# 9 ivb            2019-12-25  1130   28520
# 10 rsv_a          2019-12-29   540   15503
# 11 rsv_b          2019-12-08   211    3176
# 12 rhino          2020-02-29   393   23476

inc_df_sm2 <- left_join(inc_df_sm, first_peak_df %>% dplyr::select(-contains("inc")), by = "organism")

########################################################################################
## Pathogen-specific GAMs
########################################################################################

############################
### RSV B 19_20
############################

rsv_b_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "rsv_b")

ggplot(inc_df_sm2 %>% filter(organism == "rsv_b")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(rsv_b_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date("2020-01-14"), lty = "dashed") +
  theme_bw()

rsv_b_df_lim_scaled %>%
  filter(date <= as.Date("2020-01-14")) %>%
  pull(week_index) %>%
  max()

ggplot(rsv_b_df_lim_scaled) +
  geom_line(aes(x = date, y = rsv_b)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-01-14"), lty = "dashed") +
  theme_bw()

rsv_b_df_lim_scaled <- rsv_b_df_lim_scaled %>% filter(date < as.Date("2020-01-15"))

range(rsv_b_df_lim_scaled$date)
head(rsv_b_df_lim_scaled)
names(rsv_b_df_lim_scaled)
unique(rsv_b_df_lim_scaled$week_index)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = rsv_b_df_lim_scaled, pathogen = "rsv_b")
smat[7, ] <- c(rep(TRUE, 6), NA)
smat

fm <- paste("s(", names(rsv_b_df_lim_scaled[!(names(rsv_b_df_lim_scaled) %in% c("date", "week", "rsv_b"))]), ")", sep = "", collapse = " + ")
fm

fm <- as.formula(paste("rsv_b ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = rsv_b_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail",
  family = Gamma(link = "log")
)
summary(b1)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
rsv_b_b1 <- update(best_mod, method = "REML")

rsvb_p <- draw(rsv_b_b1, residuals = TRUE, title = "RSV B", scales = "free", nrow = 1)
rsvb_p <- rsvb_p & theme_bw()
rsvb_p

p1 <- draw(rsv_b_b1, residuals = T, select = "s(colleges)") & theme_bw()
p1
p2 <- draw(rsv_b_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rsv_b_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rsv_b_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "RSV B", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rsv_b_p_all

############################
### H1N1 19_20
############################
h1n1_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "h1n1")

peak <- inc_df_sm2 %>%
  filter(organism == "h1n1") %>%
  pull(peak_date) %>%
  unique()

ggplot(inc_df_sm2 %>% filter(organism == "h1n1")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(h1n1_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-14"), lty = "dashed") +
  theme_bw()

ggplot(h1n1_df_lim_scaled) +
  geom_line(aes(x = date, y = h1n1)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-20"), lty = "dashed") +
  theme_bw()

h1n1_df_lim_scaled <- h1n1_df_lim_scaled %>% filter(date < as.Date("2020-02-20"))

range(h1n1_df_lim_scaled$date)
h1n1_df_lim_scaled %>% slice_max(h1n1, n = 1)
range(h1n1_df_lim_scaled$week_index)
head(h1n1_df_lim_scaled)

ggplot(h1n1_df_lim_scaled) +
  geom_line(aes(x = date, y = h1n1))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = h1n1_df_lim_scaled, pathogen = "h1n1")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(h1n1_df_lim_scaled[!(names(h1n1_df_lim_scaled) %in% c("date", "week", "h1n1"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("h1n1 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = h1n1_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
h1n1_b1 <- update(best_mod, method = "REML")

h1n1_p <- draw(h1n1_b1, residuals = TRUE, scales = "free", nrow = 1)
h1n1_p <- h1n1_p & theme_bw()
h1n1_p

p1 <- draw(h1n1_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
p2 <- draw(h1n1_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(h1n1_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

h1n1_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Influenza A/H1N1", theme = theme(plot.title = element_text(size = 18, face = "bold")))
h1n1_p_all

############################
### RSV A 19_20
############################
rsv_a_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "rsv_a")

peak <- inc_df_sm2 %>%
  filter(organism == "rsv_a") %>%
  pull(peak_date) %>%
  unique()

ggplot(inc_df_sm2 %>% filter(organism == "rsv_a")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(rsv_a_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-14"), lty = "dashed") +
  theme_bw()

ggplot(rsv_a_df_lim_scaled) +
  geom_line(aes(x = date, y = rsv_a)) +
  geom_vline(xintercept = as.Date("2020-02-01"), lty = "dashed")

rsv_a_df_lim_scaled <- rsv_a_df_lim_scaled %>% filter(date < as.Date("2020-02-01"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(rsv_a_df_lim_scaled, "rsv_a")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(rsv_a_df_lim_scaled[!(names(rsv_a_df_lim_scaled) %in% c("date", "week", "rsv_a"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rsv_a ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = rsv_a_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE,
  na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rsv_a_b1 <- update(best_mod, method = "REML")
summary(rsv_a_b1)

rsva_p <- draw(rsv_a_b1, residuals = TRUE, scales = "free", nrow = 1)
rsva_p <- rsva_p & theme_bw()
rsva_p

p1 <- draw(rsv_a_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
p2 <- draw(rsv_a_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rsv_a_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rsv_a_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "RSV A", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rsv_a_p_all

############################
### IBV 19_20
############################

ivb_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "ivb")

peak <- inc_df_sm2 %>%
  filter(organism == "ivb") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "ivb")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(ivb_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-14"), lty = "dashed") +
  theme_bw()

ggplot(ivb_df_lim_scaled) +
  geom_line(aes(x = date, y = ivb)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date("2019-12-23"), lty = "dashed") +
  theme_bw()

ivb_df_lim_scaled <- ivb_df_lim_scaled %>% filter(date < as.Date("2019-12-24"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(ivb_df_lim_scaled, "ivb")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(ivb_df_lim_scaled[!(names(ivb_df_lim_scaled) %in% c("date", "week", "ivb"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("ivb ~", fm))
fm

# remove  gamma = 1.4 because of convergence issues
b1 <- mgcv::gam(
  formula = fm, data = ivb_df_lim_scaled, method = "ML", gamma = 1.3,
  select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
ivb_b1 <- update(best_mod, method = "REML")
summary(ivb_b1)

ivb_p <- draw(ivb_b1, residuals = TRUE, scales = "free", nrow = 1)
ivb_p <- ivb_p & theme_bw()
ivb_p

p1 <- draw(ivb_b1, residuals = T, select = "s(colleges)") & theme_bw()
p1

p2 <- draw(ivb_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2

p3 <- draw(ivb_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

ivb_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Influenza B", theme = theme(plot.title = element_text(size = 18, face = "bold")))
ivb_p_all

############################
### HMPV 19_20
############################
hmpv_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "hmpv") %>% filter(hmpv > 0)

peak <- inc_df_sm2 %>%
  filter(organism == "hmpv") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "hmpv")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(hmpv_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-14"), lty = "dashed") +
  theme_bw()

ggplot(hmpv_df_lim_scaled) +
  geom_line(aes(x = date, y = hmpv)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(peak), lty = "dashed") +
  theme_bw()

hmpv_df_lim_scaled <- hmpv_df_lim_scaled %>% filter(date < as.Date(peak))

range(hmpv_df_lim_scaled$date)
length(unique(hmpv_df_lim_scaled$week))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(df = hmpv_df_lim_scaled, pathogen = "hmpv")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(hmpv_df_lim_scaled[!(names(hmpv_df_lim_scaled) %in% c("date", "week", "hmpv"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hmpv ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = hmpv_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE,
  na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hmpv_b1 <- update(best_mod, method = "REML")

hmpv_p <- draw(hmpv_b1, residuals = TRUE, scales = "free", nrow = 1)
hmpv_p <- hmpv_p & theme_bw()
hmpv_p

p1 <- draw(hmpv_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
p2 <- draw(hmpv_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(hmpv_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hmpv_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hMPV", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hmpv_p_all

############################
### scov_229E_OC43 19_20
############################
scov_229E_OC43_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "scov_229E_OC43")

peak <- inc_df_sm2 %>%
  filter(organism == "scov_229E_OC43") %>%
  pull(peak_date) %>%
  unique() %>%
  max()
peak

ggplot(inc_df_sm2 %>% filter(organism == "scov_229E_OC43")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(scov_229E_OC43_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-14"), lty = "dashed") +
  theme_bw()

ggplot(scov_229E_OC43_df_lim_scaled) +
  geom_line(aes(x = date, y = scov_229E_OC43)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  theme_bw()

scov_229E_OC43_df_lim_scaled <- scov_229E_OC43_df_lim_scaled %>% filter(date < as.Date(max(peak)))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(scov_229E_OC43_df_lim_scaled, "scov_229E_OC43")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(scov_229E_OC43_df_lim_scaled[!(names(scov_229E_OC43_df_lim_scaled) %in% c("date", "week", "scov_229E_OC43"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_229E_OC43 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_229E_OC43_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_229E_OC43_b1 <- update(best_mod, method = "REML")

p <- draw(scov_229E_OC43_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
# p2 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(within_city_movement)") & theme_bw() & ggtitle("s(between-CBG movement)") & xlab("between-CBG movement")
# p2
p3 <- draw(scov_229E_OC43_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

scov_229E_OC43_p_all <- p1 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hCoV 229E + OC43", theme = theme(plot.title = element_text(size = 18, face = "bold")))

############################
### scov HKU1 NL63 19_20
############################

scov_HKU1_NL63_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "scov_HKU1_NL63")
range(scov_HKU1_NL63_df_lim_scaled$week_index)

peak <- inc_df_sm2 %>%
  filter(organism == "scov_HKU1_NL63") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "scov_HKU1_NL63")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(scov_HKU1_NL63_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-01"), lty = "dashed") +
  theme_bw()

ggplot(scov_HKU1_NL63_df_lim_scaled) +
  geom_line(aes(x = date, y = scov_HKU1_NL63)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-15"), lty = "dashed") +
  theme_bw()

scov_HKU1_NL63_df_lim_scaled <- scov_HKU1_NL63_df_lim_scaled %>% filter(date < as.Date("2020-01-01"))
range(scov_HKU1_NL63_df_lim_scaled$date)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(scov_HKU1_NL63_df_lim_scaled, "scov_HKU1_NL63")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(scov_HKU1_NL63_df_lim_scaled[!(names(scov_HKU1_NL63_df_lim_scaled) %in% c("date", "week", "time", "scov_HKU1_NL63", "child_day_care"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("scov_HKU1_NL63 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = scov_HKU1_NL63_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
scov_HKU1_NL63_b1 <- update(best_mod, method = "REML")

p <- draw(scov_HKU1_NL63_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(colleges)") & theme_bw()
p1

p2 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p2

p3 <- draw(scov_HKU1_NL63_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

scov_HKU1_NL63_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hCoV HKU1 + NL63", theme = theme(plot.title = element_text(size = 18, face = "bold")))
scov_HKU1_NL63_p_all

############################
### hpiv 1+2
############################
hpiv_1_2_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "hpiv_1_2")
range(hpiv_1_2_df_lim_scaled$week_index)

peak <- inc_df_sm2 %>%
  filter(organism == "hpiv_1_2") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "hpiv_1_2")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(hpiv_1_2_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2019-12-15"), lty = "dashed") +
  theme_bw()

ggplot(hpiv_1_2_df_lim_scaled) +
  geom_line(aes(x = date, y = hpiv_1_2)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2019-12-15"), lty = "dashed") +
  theme_bw()

hpiv_1_2_df_lim_scaled <- hpiv_1_2_df_lim_scaled %>% filter(date < as.Date("2019-12-15"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(hpiv_1_2_df_lim_scaled, "hpiv_1_2")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(hpiv_1_2_df_lim_scaled[!(names(hpiv_1_2_df_lim_scaled) %in% c("date", "week", "hpiv_1_2"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hpiv_1_2 ~", fm))
fm

b1 <- mgcv::gam(formula = fm, data = hpiv_1_2_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log"))
dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hpiv_1_2_b1 <- update(best_mod, method = "REML")

p <- draw(hpiv_1_2_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hpiv_1_2_b1, residuals = T, select = "s(colleges)") & theme_bw()
p1

p2 <- draw(hpiv_1_2_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2

p3 <- draw(hpiv_1_2_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hpiv_1_2_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hPIV 1 + 2", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hpiv_1_2_p_all

############################
### hpiv 3+4
############################
hpiv_3_4_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "hpiv_3_4")
range(hpiv_3_4_df_lim_scaled$week_index)

peak <- inc_df_sm2 %>%
  filter(organism == "hpiv_3_4") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "hpiv_3_4")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(hpiv_3_4_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2019-12-15"), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-01"), lty = "dashed") +
  theme_bw()

ggplot(hpiv_3_4_df_lim_scaled) +
  geom_line(aes(x = date, y = hpiv_3_4)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-02-01"), lty = "dashed") +
  theme_bw()

hpiv_3_4_df_lim_scaled <- hpiv_3_4_df_lim_scaled %>% filter(date < as.Date("2020-02-01"))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(hpiv_3_4_df_lim_scaled, "hpiv_3_4")
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(hpiv_3_4_df_lim_scaled[!(names(hpiv_3_4_df_lim_scaled) %in% c("date", "week", "hpiv_3_4"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("hpiv_3_4 ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = hpiv_3_4_df_lim_scaled, method = "ML",
  gamma = 1.4, select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
hpiv_3_4_b1 <- update(best_mod, method = "REML")

p <- draw(hpiv_3_4_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(hpiv_3_4_b1, residuals = T, select = "s(colleges)") & theme_bw() & ggtitle("s(colleges)")
p1
p2 <- draw(hpiv_3_4_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(hpiv_3_4_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

hpiv_3_4_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "hPIV 3 + 4", theme = theme(plot.title = element_text(size = 18, face = "bold")))
hpiv_3_4_p_all

############################
### rhinovirus
############################
rhino_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "rhino")

peak <- inc_df_sm2 %>%
  filter(organism == "rhino") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "rhino")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(rhino_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  # geom_vline(xintercept = as.Date("2019-12-15"),lty="dashed")+
  # geom_vline(xintercept = as.Date("2020-02-01"),lty="dashed")+
  theme_bw()

ggplot(rhino_df_lim_scaled) +
  geom_line(aes(x = date, y = rhino)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-01-30"), lty = "dashed") +
  theme_bw()

rhino_df_lim_scaled <- rhino_df_lim_scaled %>% filter(date < as.Date("2020-01-30"))
range(rhino_df_lim_scaled$date)
range(rhino_df_lim_scaled$week_index)

ggplot(rhino_df_lim_scaled) +
  geom_line(aes(x = date, y = rhino))

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(rhino_df_lim_scaled, "rhino")
smat
smat[7, ] <- c(rep(TRUE, 6), NA)

fm <- paste("s(", names(rhino_df_lim_scaled[!(names(rhino_df_lim_scaled) %in% c("date", "week", "rhino"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("rhino ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = rhino_df_lim_scaled, method = "ML", gamma = 1.4,
  select = TRUE, na.action = "na.fail", family = Gamma(link = "log")
)

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc")
subset(dd, delta < 4)

best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
rhino_b1 <- update(best_mod, method = "REML")

p <- draw(rhino_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(rhino_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1

p2 <- draw(rhino_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p2

p3 <- draw(rhino_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

rhino_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Rhinovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
rhino_p_all

############################
### adenovirus
############################
adeno_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "adeno")
range(adeno_df_lim_scaled$date)

peak <- inc_df_sm2 %>%
  filter(organism == "adeno") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "adeno")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(adeno_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-01-01"), lty = "dashed") +
  # geom_vline(xintercept = as.Date("2019-12-15"),lty="dashed")+
  # geom_vline(xintercept = as.Date("2020-02-01"),lty="dashed")+
  theme_bw()

ggplot(adeno_df_lim_scaled) +
  geom_line(aes(x = date, y = adeno)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date("2019-10-20"), lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-01-01"), lty = "dashed") +
  theme_bw()

adeno_df_lim_scaled <- adeno_df_lim_scaled %>% filter(date < as.Date("2020-02-01"))

fm <- paste("s(", names(adeno_df_lim_scaled[!(names(adeno_df_lim_scaled) %in% c("date", "week", "time", "adeno"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("adeno ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = adeno_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE,
  family = Gamma(link = "log"), na.action = "na.fail"
)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(adeno_df_lim_scaled, "adeno")
smat[7, ] <- c(rep(TRUE, 6), NA)
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
adeno_b1 <- update(best_mod, method = "REML")

p <- draw(adeno_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(adeno_b1, residuals = T, select = "s(elementary_and_secondary_schools)") & theme_bw() & ggtitle("s(elem. and high schools)") & xlab("elem. and high schools")
p1
p2 <- draw(adeno_b1, residuals = T, select = "s(religious_orgs)") & theme_bw() & ggtitle("s(religious services)") & xlab("religious services")
p2
p3 <- draw(adeno_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

adeno_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Adenovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
adeno_p_all

############################
### enterovirus
############################
entero_df_lim_scaled <- make_gam_df(com_19_20_df_lm, "entero")
range(entero_df_lim_scaled$date)

peak <- inc_df_sm2 %>%
  filter(organism == "entero") %>%
  pull(peak_date) %>%
  unique()
peak

ggplot(inc_df_sm2 %>% filter(organism == "entero")) +
  geom_line(aes(x = date, y = inc)) +
  geom_vline(xintercept = as.Date(min(entero_df_lim_scaled$date))) +
  geom_vline(xintercept = as.Date(peak)) +
  geom_vline(xintercept = as.Date("2020-02-01"), lty = "dashed") +
  # geom_vline(xintercept = as.Date("2019-12-15"),lty="dashed")+
  # geom_vline(xintercept = as.Date("2020-02-01"),lty="dashed")+
  theme_bw()

ggplot(entero_df_lim_scaled) +
  geom_line(aes(x = date, y = entero)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_vline(xintercept = as.Date("2019-10-20"), lty = "dashed") +
  geom_vline(xintercept = as.Date(max(peak)), lty = "dashed") +
  geom_vline(xintercept = as.Date("2020-01-28"), lty = "dashed") +
  theme_bw()

entero_df_lim_scaled <- entero_df_lim_scaled %>% filter(date < as.Date("2020-01-30"))

fm <- paste("s(", names(entero_df_lim_scaled[!(names(entero_df_lim_scaled) %in% c("date", "week", "time", "entero"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste("entero ~", fm))
fm

b1 <- mgcv::gam(
  formula = fm, data = entero_df_lim_scaled, method = "ML", gamma = 1.4, select = TRUE,
  family = Gamma(link = "log"), na.action = "na.fail"
)

# limit combinations of variables to pairs that are not highly collinear
smat <- subset_cor(entero_df_lim_scaled, "entero")
smat[7, ] <- c(rep(TRUE, 6), NA)
smat

dd <- dredge(b1, m.lim = c(1, 3), fixed = "s(week_index)", rank = "AICc", subset = smat)
subset(dd, delta < 4)
best_mod <- get.models(dd, subset = 1)[[1]]
summary(best_mod)
entero_b1 <- update(best_mod, method = "REML")

p <- draw(entero_b1, residuals = TRUE, scales = "free", nrow = 1)
p <- p & theme_bw()
p

p1 <- draw(entero_b1, residuals = T, select = "s(fb_leaving_home_custom)") & theme_bw() & ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1
p2 <- draw(entero_b1, residuals = T, select = "s(out_of_state_movement)") & theme_bw() & ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(entero_b1, residuals = T, select = "s(week_index)") & theme_bw() & ggtitle("s(week)") & xlab("week")
p3

entero_p_all <- p1 + p2 + p3 + plot_layout(ncol = 3) +
  plot_annotation(title = "Enterovirus", theme = theme(plot.title = element_text(size = 18, face = "bold")))
entero_p_all

########################################################################################
## combine all plots into one figure
########################################################################################

all_path_2019_2020 <- plot_grid(h1n1_p_all,
  ivb_p_all,
  rsv_a_p_all,
  rsv_b_p_all,
  hmpv_p_all,
  scov_229E_OC43_p_all,
  scov_HKU1_NL63_p_all,
  hpiv_1_2_p_all,
  hpiv_3_4_p_all,
  rhino_p_all,
  entero_p_all,
  adeno_p_all,
  ncol = 2, byrow = F
)
all_path_2019_2020

save_plot(all_path_2019_2020, base_height = 18, base_width = 20, filename = "figures/fig_s9_partial_effects_all_pathogens_2019_2020_updated.png")
