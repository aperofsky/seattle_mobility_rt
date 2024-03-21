########################################################################################
## Predictive model of Rhinovirus Rt
########################################################################################

library(readr)
library(dplyr)
library(tidyverse)
library(glue)
library(cdcfluview)
library(forecast)
library(tibbletime)
library(caret)
library(lubridate)
library(zoo)
library(glmnet)
library(Metrics)
library(MLmetrics)
library(ggplot2)
library(cowplot)
library(ggExtra)
library(tidyquant)
library(ggpubr)

########################################################################################
## Import data
########################################################################################
rt_df <- read_csv("2_Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

climate_df <- read_csv("6_LASSO_Predictive_Models/1_Climate_data/seattle_daily_weather_variables_2023_01_30.csv")

## interpolate missing Precipitation data
climate_df <- climate_df %>%
  complete(date = seq.Date(as.Date("2018-01-01"), as.Date("2022-09-30"), by = "day")) %>%
  mutate(
    DailyPrecipitation = zoo::na.approx(DailyPrecipitation, maxgap = 4),
    DailyAverageRelativeHumidity = zoo::na.approx(DailyAverageRelativeHumidity, maxgap = 4),
    DailyAverageWetBulbTemperature = zoo::na.approx(DailyAverageWetBulbTemperature, maxgap = 4)
  )
head(climate_df)
range(climate_df$date)

combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()
combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined[!complete.cases(combined), ]

combined <- left_join(combined, climate_df, by = "date")
head(combined)

combined <- combined %>%
  mutate(across(c(within_neighborhood_movement:DailyPrecipitation), ~ lag(.x, n = 1), .names = paste0("{.col}_lag", 1))) %>%
  mutate(across(c(within_neighborhood_movement:DailyPrecipitation), ~ lag(.x, n = 7), .names = paste0("{.col}_lag", 7))) %>%
  mutate(across(c(within_neighborhood_movement:DailyPrecipitation), ~ lag(.x, n = 14), .names = paste0("{.col}_lag", 14)))

max_date <- combined %>%
  filter(adeno > 0) %>%
  slice_max(date) %>%
  pull(date)
max_date

min_date <- combined %>%
  filter(adeno > 0) %>%
  slice_min(date) %>%
  pull(date)
min_date # "2018-11-21"

combined %>% filter(!is.na(fb_leaving_home_custom_lag14))

cols2 <- c(
  # "within_neighborhood_movement",
  "within_city_movement",
  "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  # "groceries_and_pharmacies",
  # "transit",
  "religious_orgs",
  # "oxford_stringency_index",
  # "mask_wearing_final",
  "child_day_care",
  "elementary_and_secondary_schools",
  "colleges",
  # "Perc_Full_Vax",
  # "Ancestral",
  # "Alpha",
  # "Delta",
  # "omicron_ba1",
  # "omicron_ba2",
  "DailyAverageRelativeHumidity",
  "DailyPrecipitation",
  "DailyAverageWetBulbTemperature"
)
cols2

combined %>%
  select(date, matches(cols2), adeno, covid) %>%
  dplyr::select(date, adeno, covid, contains("lag"))

adeno_df <- combined %>%
  select(date, matches(cols2), adeno, covid) %>%
  dplyr::select(date, adeno, covid, contains("lag")) %>%
  filter(date >= as.Date("2019-01-22") & date <= as.Date("2022-05-22")) %>%
  complete(date = seq.Date(as.Date("2019-01-22"), as.Date("2022-05-22"), by = "day")) %>%
  filter(!is.na(fb_leaving_home_custom_lag1))

data.table::data.table(adeno_df)
adeno_df <- adeno_df %>% replace_na(list(covid = 0))
data.table::data.table(adeno_df[!complete.cases(adeno_df), ])
########################################################################################
## Moving window LASSO predictive model
########################################################################################
## limit to after vaccination and Delta emergence
adeno_df_scaled <- adeno_df
adeno_df_scaled[!complete.cases(adeno_df_scaled), ]
head(adeno_df_scaled) # 2019-01-22
data.table::data.table(adeno_df_scaled[!complete.cases(adeno_df_scaled), ])

window_train <- 30 # number of days of the training window
w <- c(rep(1, window_train - 28), rep(2, 28)) # observation weights (last month has double the weight)
lag_AR <- 1 # minimum lag (in days) of the autoregressive term (sometimes surveillance data may have more than 1 day lag)
window_AR <- 14 # time window of the autoregressive terms (in days)
t0_train <- 15 # starting day of the training
dt_predict <- 7 # prediction ahead relative to the last data point of the training set
alpha <- 1 # weight between L1 and L2 penalty alpha = 1 is Lasso regression
nfolds <- 10 # nfolds crossvalidation to evaluate hyperparameter lambda for elastic net regression
s <- "lambda.min" # choose lambda either min cross validation error (lambda.min) or most regularized (lambda.1se)

names(adeno_df_scaled)
var_target <- "adeno" # dependent variable
var_regressors_GT <- select(adeno_df_scaled, -date, -adeno, -covid, -c(DailyAverageRelativeHumidity_lag1:DailyAverageWetBulbTemperature_lag14)) %>% names()
var_regressors_GT
var_regressors_GT_COVID <- select(adeno_df_scaled, -date, -adeno, -c(DailyAverageRelativeHumidity_lag1:DailyAverageWetBulbTemperature_lag14)) %>% names()
var_regressors_GT_COVID
var_regressors_climate <- select(adeno_df_scaled, -date, -adeno, -covid) %>% names()
var_regressors_climate
var_regressors_AR <- sprintf("AR%02d", lag_AR:(lag_AR + window_AR - 1)) # regressors, autoregressive
var_regressors_AR
var_regressors_ARGO <- unique(c(var_regressors_AR, var_regressors_GT_COVID, var_regressors_climate)) # all regressors, ARGO
var_regressors_ARGO
var_regressors_GT_AR <- unique(c(var_regressors_AR, var_regressors_GT))
var_regressors_GT_AR
climate_var <- var_regressors_climate[grepl("Daily", var_regressors_climate)]
var_regressors_climate_AR <- unique(c(var_regressors_AR, climate_var))
var_regressors_climate_AR
var_regressors_mob_climate_AR <- unique(c(var_regressors_AR, var_regressors_GT, climate_var))
var_regressors_mob_climate_AR
var_predicts <- sprintf("predict_%dd_ahead", 1:dt_predict)
var_data <- sprintf("data_%dd_ahead", 1:dt_predict)

min_date <- min(adeno_df_scaled$date)

# make reference dates
reference_dates <- data.frame(num_days = seq(from = -1500, to = 1500, by = 1))
reference_dates[["date"]] <- reference_dates[["num_days"]] + ymd(min_date)
reference_dates[["epi_year"]] <- epiyear(reference_dates[["date"]])
reference_dates[["epi_week"]] <- epiweek(reference_dates[["date"]])

rdf <- left_join(adeno_df_scaled, reference_dates)
day_max <- max(rdf$num_days)

t0_predict_start <- t0_train + window_train
tmax_predict_start <- day_max - dt_predict + 1

# Set up lags
lags <- c(1:14)

# Name the columns that will contain the lags, with appropriate index
lag_names <- glue::glue('AR{str_pad(lags, nchar(max(lags)), pad = "0")}')

# Create list of lag functions, as eval-ed/parse-d labmda functions
lag_functions <-
  map(lags, ~ eval(parse(text = glue::glue("~ dplyr::lag(.x, {.x})")))) %>%
  set_names(lag_names)

adeno_agg <- rdf %>%
  dplyr::select(date, num_days, epi_year, epi_week, adeno) %>%
  mutate_at(vars(adeno), .funs = lag_functions)

adeno_agg <- gather(adeno_agg, key = "metric", value = "adeno", adeno:AR14)

temp_adeno <- adeno_agg %>%
  rename(val = adeno) %>%
  filter(num_days >= 15 & num_days <= day_max)

sgtrends_adeno <- adeno_df_scaled %>%
  left_join(reference_dates) %>%
  dplyr::select(
    date, all_of(unique(c(var_regressors_GT, var_regressors_GT_COVID, var_regressors_climate))),
    num_days, epi_year, epi_week
  ) %>%
  pivot_longer(
    cols = all_of(unique(c(var_regressors_GT, var_regressors_GT_COVID, var_regressors_climate))),
    names_to = "metric", values_to = "val"
  ) %>%
  filter(num_days >= 15 & num_days <= day_max)
unique(sgtrends_adeno$metric)

sg <- unique(sgtrends_adeno$metric)
sg

col_order <- c("date", "num_days", "epi_year", "epi_week", "adeno", sprintf("AR%02d", 1:14), sg)
glm_adeno <- rbind(temp_adeno, sgtrends_adeno)
glm_adeno <- spread(glm_adeno, key = "metric", value = "val")
glm_adeno <- glm_adeno[col_order]
names(glm_adeno)

res_GT <- data.frame()
res_GT_AR <- data.frame()
res_GT_COVID <- data.frame()
res_climate_AR <- data.frame()
res_climate <- data.frame()
res_mob_climate_AR <- data.frame()
res_AR <- data.frame()
res_ARGO <- data.frame()
# t_predict_start=45
for (t_predict_start in t0_predict_start:tmax_predict_start) {
  print(t_predict_start)
  t_predict_end <- t_predict_start + dt_predict - 1
  t_train_start <- t_predict_start - window_train
  t_train_end <- t_predict_start - 1
  date_train_end <- reference_dates[reference_dates$num_days == t_train_end, c("date")]

  glm_train <- filter(glm_adeno, num_days >= t_train_start & num_days <= t_train_end)
  range(glm_train$date)
  glm_predict <- filter(glm_adeno, num_days >= t_predict_start & num_days <= t_predict_end) # predict up to 7 days ahead
  range(glm_predict$date)
  # names(glm_predict)

  x_train_GT <- data.matrix(glm_train[, var_regressors_GT])
  x_train_GT_AR <- data.matrix(glm_train[, var_regressors_GT_AR])
  x_train_GT_COVID <- data.matrix(glm_train[, var_regressors_GT_COVID])
  x_train_climate <- data.matrix(glm_train[, var_regressors_climate])
  x_train_climate_AR <- data.matrix(glm_train[, var_regressors_climate_AR])
  x_train_mob_climate_AR <- data.matrix(glm_train[, var_regressors_mob_climate_AR])
  x_train_AR <- data.matrix(glm_train[, var_regressors_AR])
  x_train_ARGO <- data.matrix(glm_train[, var_regressors_ARGO])

  x_predict_GT <- data.matrix(glm_predict[, var_regressors_GT])
  x_predict_GT_AR <- data.matrix(glm_predict[, var_regressors_GT_AR])
  x_predict_GT_COVID <- data.matrix(glm_predict[, var_regressors_GT_COVID])
  x_predict_climate <- data.matrix(glm_predict[, var_regressors_climate])
  x_predict_climate_AR <- data.matrix(glm_predict[, var_regressors_climate_AR])
  x_predict_mob_climate_AR <- data.matrix(glm_predict[, var_regressors_mob_climate_AR])
  x_predict_AR <- data.matrix(glm_predict[, var_regressors_AR])
  x_predict_ARGO <- data.matrix(glm_predict[, var_regressors_ARGO])

  y_train <- glm_train[, var_target] %>% pull(adeno)
  # range(glm_train$date)
  y_train <- as.numeric(y_train)

  cvfit_GT <- cv.glmnet(x_train_GT, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_GT_AR <- cv.glmnet(x_train_GT_AR, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_GT_COVID <- cv.glmnet(x_train_GT_COVID, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_climate <- cv.glmnet(x_train_climate, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_climate_AR <- cv.glmnet(x_train_climate_AR, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_mob_climate_AR <- cv.glmnet(x_train_mob_climate_AR, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_AR <- cv.glmnet(x_train_AR, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_ARGO <- cv.glmnet(x_train_ARGO, y_train, weights = w, type.measure = "mse", nfolds = nfolds, alpha = alpha)

  coef_GT <- cvfit_GT %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_GT_AR <- cvfit_GT_AR %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_GT_COVID <- cvfit_GT_COVID %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_climate <- cvfit_climate %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_climate_AR <- cvfit_climate_AR %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_mob_climate_AR <- cvfit_mob_climate_AR %>%
    coef(s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_AR <- coef(cvfit_AR, s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  coef_ARGO <- coef(cvfit_ARGO, s = s) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rename(intercept = "(Intercept)")

  predict_GT <- cvfit_GT %>%
    predict(newx = x_predict_GT, s = s) %>%
    t() %>%
    data.frame()

  predict_GT_AR <- cvfit_GT_AR %>%
    predict(newx = x_predict_GT_AR, s = s) %>%
    t() %>%
    data.frame()

  predict_GT_COVID <- cvfit_GT_COVID %>%
    predict(newx = x_predict_GT_COVID, s = s) %>%
    t() %>%
    data.frame()

  predict_climate <- cvfit_climate %>%
    predict(newx = x_predict_climate, s = s) %>%
    t() %>%
    data.frame()

  predict_climate_AR <- cvfit_climate_AR %>%
    predict(newx = x_predict_climate_AR, s = s) %>%
    t() %>%
    data.frame()

  predict_mob_climate_AR <- cvfit_mob_climate_AR %>%
    predict(newx = x_predict_mob_climate_AR, s = s) %>%
    t() %>%
    data.frame()

  predict_AR <- cvfit_AR %>%
    predict(newx = x_predict_AR, s = s) %>%
    t() %>%
    data.frame()

  predict_ARGO <- cvfit_ARGO %>%
    predict(newx = x_predict_ARGO, s = s) %>%
    t() %>%
    data.frame()

  colnames(predict_GT) <- var_predicts
  colnames(predict_GT_AR) <- var_predicts
  colnames(predict_GT_COVID) <- var_predicts
  colnames(predict_climate) <- var_predicts
  colnames(predict_climate_AR) <- var_predicts
  colnames(predict_mob_climate_AR) <- var_predicts
  colnames(predict_AR) <- var_predicts
  colnames(predict_ARGO) <- var_predicts

  data_epi <- glm_predict[, "adeno"] %>%
    t() %>%
    data.frame()

  colnames(data_epi) <- var_data

  res_GT <- bind_rows(res_GT, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_GT, coef_GT))
  head(res_GT)

  res_GT_AR <- bind_rows(res_GT_AR, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_GT_AR, coef_GT_AR))
  head(res_GT_AR)

  res_GT_COVID <- bind_rows(res_GT_COVID, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_GT_COVID, coef_GT_COVID))
  head(res_GT_COVID)

  res_climate <- bind_rows(res_climate, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_climate, coef_climate))
  head(res_climate)

  res_climate_AR <- bind_rows(res_climate_AR, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_climate_AR, coef_climate_AR))
  head(res_climate_AR)

  res_mob_climate_AR <- bind_rows(res_mob_climate_AR, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_mob_climate_AR, coef_mob_climate_AR))
  head(res_mob_climate_AR)

  res_AR <- bind_rows(res_AR, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_AR, coef_AR))

  res_ARGO <- bind_rows(res_ARGO, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_ARGO, coef_ARGO))
}

########################################################################################
## Summarize model results
########################################################################################

head(res_ARGO)
tail(res_ARGO)

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
coef_ARGO$regressor <- gsub("_", ".", coef_ARGO$regressor)
orders <- rev(c("intercept", gsub("_", ".", var_regressors_ARGO)))

# look at coefficients by correlation
summary_coef <- coef_ARGO %>%
  group_by(regressor) %>%
  summarise(avg = mean(abs(coef), na.rm = TRUE)) %>%
  ungroup()

summary_coef %>%
  arrange(-abs(avg)) %>%
  filter(!grepl("AR", regressor)) %>%
  print(n = 40)

summary_coef[summary_coef$avg == 0, ]

summary_coef <- coef_ARGO %>%
  filter(coef != 0) %>%
  group_by(regressor) %>%
  tally() %>%
  arrange(-n)
summary_coef %>%
  print(n = 50)

########################################################################################
## Plots
########################################################################################

########################################################################################
## Plot coefficient values at each time point
########################################################################################

patterns <- c("AR", "intercept")
keep <- unique(coef_ARGO$regressor)[!grepl(paste(patterns, collapse = "|"), unique(coef_ARGO$regressor))]
keep
ar_terms <- unique(coef_ARGO$regressor)[c(1:15)]

coef_ARGO %>%
  filter(regressor %in% ar_terms) %>%
  pull(coef) %>%
  range()

label_df <- data.frame(
  var = unique(coef_ARGO$regressor),
  label = c(
    unique(coef_ARGO$regressor)[1:15],
    "SARS-CoV-2 Rt",
    # "within-neighborhood movement",
    # "between-neighborhood movement",
    "between-neighborhood movement, Lag 1",
    "between-neighborhood movement, Lag 7",
    "between-neighborhood movement, Lag 14",
    # "influx visitors other WA counties",
    "influx visitors other WA counties, Lag 1",
    "influx visitors other WA counties, Lag 7",
    "influx visitors other WA counties, Lag 14",
    # "influx out-of-state visitors",
    "influx out-of-state visitors, Lag 1",
    "influx out-of-state visitors, Lag 7",
    "influx out-of-state visitors, Lag 14",
    # "% devices leaving home",
    "% devices leaving home, Lag 1",
    "% devices leaving home, Lag 7",
    "% devices leaving home, Lag 14",
    # "restaurants",
    "restaurants, Lag 1",
    "restaurants, Lag 7",
    "restaurants, Lag 14",
    # "groceries and pharmacies",
    # "transit",
    # "religious organizations",
    "religious organizations, Lag 1",
    "religious organizations, Lag 7",
    "religious organizations, Lag 14",
    # "Oxford Stringency Index",
    # "Oxford Stringency Index, Lag 1",
    # "Oxford Stringency Index, Lag 7",
    # "Oxford Stringency Index, Lag 14",
    # "% not wearing masks",
    # "% not wearing masks, Lag 1",
    # "% not wearing masks, Lag 7",
    # "% not wearing masks, Lag 14",
    # "child daycare",
    "child daycare, Lag 1",
    "child daycare, Lag 7",
    "child daycare, Lag 14",
    # "elementary and high schools",
    "elementary and high schools, Lag 1",
    "elementary and high schools, Lag 7",
    "elementary and high schools, Lag 14",
    # "colleges",
    "colleges, Lag 1",
    "colleges, Lag 7",
    "colleges, Lag 14",
    # "Relative Humidity",
    "Relative Humidity, Lag 1",
    "Relative Humidity, Lag 7",
    "Relative Humidity, Lag 14",
    # "Precipitation",
    "Precipitation, Lag 1",
    "Precipitation, Lag 7",
    "Precipitation, Lag 14",
    # "Temperature",
    "Temperature, Lag 1",
    "Temperature, Lag 7",
    "Temperature, Lag 14"
  )
)

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
range(coef_ARGO$coef)

coef_ARGO[coef_ARGO$coef < -0.006, "coef"] <- -0.006
coef_ARGO[coef_ARGO$coef > 0.006, "coef"] <- 0.006
coef_ARGO$regressor <- gsub("_", ".", coef_ARGO$regressor)
plot_df <- left_join(coef_ARGO, label_df, by = c("regressor" = "var"))

orders <- c(
  unique(plot_df$label)[1:15],
  "SARS-CoV-2 Rt",
  # "within-neighborhood movement",
  # "between-neighborhood movement",
  "between-neighborhood movement, Lag 1",
  "between-neighborhood movement, Lag 7",
  "between-neighborhood movement, Lag 14",
  # "influx visitors other WA counties",
  "influx visitors other WA counties, Lag 1",
  "influx visitors other WA counties, Lag 7",
  "influx visitors other WA counties, Lag 14",
  # "influx out-of-state visitors",
  "influx out-of-state visitors, Lag 1",
  "influx out-of-state visitors, Lag 7",
  "influx out-of-state visitors, Lag 14",
  # "% devices leaving home",
  "% devices leaving home, Lag 1",
  "% devices leaving home, Lag 7",
  "% devices leaving home, Lag 14",
  # "restaurants",
  "restaurants, Lag 1",
  "restaurants, Lag 7",
  "restaurants, Lag 14",
  # "groceries and pharmacies",
  # "transit",
  # "religious organizations",
  "religious organizations, Lag 1",
  "religious organizations, Lag 7",
  "religious organizations, Lag 14",
  # "Oxford Stringency Index",
  # "Oxford Stringency Index, Lag 1",
  # "Oxford Stringency Index, Lag 7",
  # "Oxford Stringency Index, Lag 14",
  # "% not wearing masks",
  # "% not wearing masks, Lag 1",
  # "% not wearing masks, Lag 7",
  # "% not wearing masks, Lag 14",
  # "child daycare",
  "child daycare, Lag 1",
  "child daycare, Lag 7",
  "child daycare, Lag 14",
  # "elementary and high schools",
  "elementary and high schools, Lag 1",
  "elementary and high schools, Lag 7",
  "elementary and high schools, Lag 14",
  # "colleges",
  "colleges, Lag 1",
  "colleges, Lag 7",
  "colleges, Lag 14",
  # "Relative Humidity",
  "Relative Humidity, Lag 1",
  "Relative Humidity, Lag 7",
  "Relative Humidity, Lag 14",
  # "Precipitation",
  "Precipitation, Lag 1",
  "Precipitation, Lag 7",
  "Precipitation, Lag 14",
  # "Temperature",
  "Temperature, Lag 1",
  "Temperature, Lag 7",
  "Temperature, Lag 14"
)

p <- ggplot(
  plot_df %>% filter(regressor %in% unique(coef_ARGO$regressor)[16:length(unique(coef_ARGO$regressor))]),
  aes(x = date, y = factor(label, level = orders), fill = coef)
) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(
    name = "Coefficient", high = scales::muted("red"), mid = "white",
    low = scales::muted("blue")
  ) +
  ylab("") +
  xlab("Date") +
  scale_x_date(expand = c(0, 0), date_labels = "%b %y", date_breaks = "4 months") +
  scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 18)
p

p <- p + theme(legend.position = "right") +
  theme(strip.background = element_rect(colour = "white")) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 10)) +
  labs(x = "", y = "") +
  removeGrid() # ggExtra
p

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
range(coef_ARGO$coef)
coef_ARGO %>%
  group_by(regressor) %>%
  summarize(
    min = min(coef),
    max = max(coef)
  ) %>%
  arrange(max) %>%
  print(n = 60)

coef_ARGO[coef_ARGO$coef < -0.08, "coef"] <- -0.08
coef_ARGO[coef_ARGO$coef > 0.08, "coef"] <- 0.08
coef_ARGO$regressor <- gsub("_", ".", coef_ARGO$regressor)
plot_df <- left_join(coef_ARGO, label_df, by = c("regressor" = "var"))

q <- ggplot(
  plot_df %>% filter(!regressor %in% keep),
  aes(x = date, y = factor(label, level = rev(orders)), fill = coef)
) +
  geom_tile(color = "white", lwd = 0.1) +
  scale_fill_gradient2(
    name = "Coefficient", high = scales::muted("red"), mid = "white",
    low = scales::muted("blue")
  ) +
  ylab("") +
  xlab("Date") +
  scale_x_date(expand = c(0, 0), date_labels = "%b %y", date_breaks = "4 months") +
  theme_minimal(base_size = 18)
q

q <- q + theme(legend.position = "right") +
  theme(strip.background = element_rect(colour = "white")) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 10)) +
  labs(x = "", y = "") +
  removeGrid() # ggExtra
q
plot_grid(q, p, nrow = 2, align = "v")

########################################################################################
## Figure S18: Overlay model predictions on observed SC2 Rt and incidence
########################################################################################

res_ARGO2 <- res_ARGO %>%
  mutate(data_7d_ahead = as.numeric(as.character(data_7d_ahead)))

head(rt_df)
unique(rt_df$organism)


ci_df <- rt_df %>%
  filter(organism == "Adenovirus") %>%
  dplyr::select(date, median, lower, upper, level) %>%
  filter(level == 90) %>%
  filter(date >= as.Date("2019-01-01"))
ci_df$date <- as.Date(ci_df$date)


full_df <- readRDS("2_Epidemia_Models/Epidemia_Models_Biowulf/extended_daily_ILI_per_pos_scaled_for_select_pathogens_2023_10_22.rds")

adeno_results <- full_df %>%
  filter(organism == "Adenovirus") %>%
  droplevels()

adeno_results <- adeno_results %>%
  arrange(organism, date) %>%
  group_by(organism) %>%
  mutate(
    all_sites_per_pos_ili_scaled_sum_mv_avg =
      zoo::rollmean(all_sites_per_pos_ili_scaled_sum,
        k = 15, align = "center",
        fill = NA
      )
  ) %>%
  ungroup()

coeff <- 20
r <- ggplot() +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  geom_rect(aes(xmin = as.Date("2019-01-22"), xmax = as.Date("2019-03-13"), ymin = -Inf, ymax = Inf),
    alpha = 0.1,
    fill = "yellow"
  ) +
  geom_rect(
    aes(
      xmin = as.Date("2020-03-23"),
      xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf
    ),
    alpha = 0.1,
    fill = "orange"
  ) +
  geom_line(
    data = adeno_results %>% filter(date >= as.Date("2019-01-22") & date <= as.Date("2022-05-16")),
    aes(x = date, y = (100 * all_sites_per_pos_ili_scaled_sum_mv_avg) / coeff),
    alpha = 0.8, lwd = 1, lty = "solid", color = "#117733"
  ) +
  geom_ribbon(
    data = adeno_results %>% filter(date >= as.Date("2019-01-22") & date <= as.Date("2022-05-16")),
    aes(x = date, ymin = 0, ymax = (100 * all_sites_per_pos_ili_scaled_sum_mv_avg) / coeff), fill = "#117733", alpha = 0.4
  ) +
  geom_line(
    data = ci_df %>% filter(date >= as.Date("2019-01-22") & date <= as.Date("2022-05-16")),
    aes(x = date, y = median), color = "#ABDDA4", alpha = 0.5
  ) +
  geom_ribbon(
    data = ci_df %>% filter(date >= as.Date("2019-01-22") & date <= as.Date("2022-05-16")),
    aes(x = date, ymin = lower, ymax = upper), fill = "#ABDDA4", alpha = 0.3
  ) +
  geom_line(
    data = res_GT_AR %>% mutate(new_date = date + 7),
    aes(x = new_date, y = predict_7d_ahead, colour = "SG-AR"), linewidth = 1
  ) +
  geom_line(
    data = res_climate_AR %>% mutate(new_date = date + 7),
    aes(x = new_date, y = predict_7d_ahead, colour = "Climate-AR"), linewidth = 1
  ) +
  geom_line(
    data = res_mob_climate_AR %>% mutate(new_date = date + 7),
    aes(x = new_date, y = predict_7d_ahead, colour = "Mobility-Climate-AR"), linewidth = 1
  ) +
  geom_line(
    data = res_ARGO %>% mutate(new_date = date + 7),
    aes(x = new_date, y = predict_7d_ahead, colour = "AR-SG-COVID-Climate"), linewidth = 1
  ) +
  geom_line(
    data = res_AR %>% mutate(new_date = date + 7),
    aes(x = new_date, y = predict_7d_ahead, colour = "AR"), linewidth = 1
  ) +
  geom_line(
    data = res_ARGO2 %>% mutate(new_date = date + 7),
    aes(x = new_date, y = data_7d_ahead, colour = "Adenovirus Rt"), linewidth = 1.7
  ) +
  labs(x = "Date", y = "Rt", color = "Model") +
  theme_classic(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = as.Date("2019-03-13"), lty = "dashed", color = "yellow") +
  scale_color_manual(
    name = NULL,
    values = c("#ABDDA4", "#FDAE61", "#2B83BA", "navy", "purple", "#D7191C"),
    labels = c("Adenovirus Rt", "AR", "AR + Mobility", "AR + Climate", "AR + Mobility + Climate", "AR + Mobility + Climate + SARS-CoV-2 Rt"),
    breaks = c("Adenovirus Rt", "AR", "SG-AR", "Climate-AR", "Mobility-Climate-AR", "AR-SG-COVID-Climate")
  ) +
  background_grid() +
  scale_x_date(expand = c(0, 0), date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * coeff, name = "Scaled Incidence"),
    expand = c(0, 0)
  )
r

########################################################################################
## Predictive accuracy: Table S3
########################################################################################
model_predict_df <- bind_rows(
  res_ARGO %>% mutate(model = "AR + Mobility + Climate + SARS-CoV-2 Rt"),
  res_AR %>% mutate(model = "AR"),
  res_GT_AR %>% mutate(model = "AR + Mobility"),
  res_climate_AR %>% mutate(model = "AR + Climate"),
  res_mob_climate_AR %>% mutate(model = "AR + Mobility + Climate"),
)

accuracy_table <-
  model_predict_df %>%
  group_by(model) %>%
  summarize(
    rmse = rmse(predicted = predict_7d_ahead, actual = data_7d_ahead),
    mae = mae(predicted = predict_7d_ahead, actual = data_7d_ahead)
  )

baseline_rmse <- accuracy_table$rmse[1]
baseline_mae <- accuracy_table$mae[1]
accuracy_table$rmse[3] # 0.04018487

accuracy_table %>%
  mutate(
    rmse_perc_diff_from_AR = 100 * (rmse - baseline_rmse) / baseline_rmse,
    mae_perc_diff_from_AR = 100 * (mae - baseline_mae) / baseline_mae
  ) %>%
  arrange(mae_perc_diff_from_AR)
# model                                     rmse    mae rmse_perc_diff_from_AR mae_perc_diff_from_AR
# <chr>                                    <dbl>  <dbl>                  <dbl>                 <dbl>
# 1 AR                                      0.0290 0.0205                   0                     0
# 2 AR + Climate                            0.0313 0.0221                   8.00                  7.94
# 3 AR + Mobility                           0.0402 0.0301                  38.7                  46.7
# 4 AR + Mobility + Climate                 0.0414 0.0310                  43.0                  51.4
# 5 AR + Mobility + Climate + SARS-CoV-2 Rt 0.0424 0.0315                  46.4                  53.8

all_residuals <- model_predict_df %>%
  mutate(
    error = data_7d_ahead - predict_7d_ahead, # Yt - Ft
    percent_error = (predict_7d_ahead - data_7d_ahead) / data_7d_ahead,
    year_mon = as.yearmon(date)
  )
unique(all_residuals$model)

all_residuals_lim <- all_residuals %>%
  filter(grepl("AR", model)) %>%
  droplevels()

unique(all_residuals_lim$model)
range(all_residuals_lim$error)
range(all_residuals_lim$percent_error) #-0.1614007  0.4278382

p <- ggplot() +
  geom_rect(aes(xmin = as.Date("2019-01-22"), xmax = as.Date("2019-03-13"), ymin = -Inf, ymax = Inf),
    alpha = 0.1,
    fill = "yellow"
  ) +
  geom_rect(
    aes(
      xmin = as.Date("2020-03-23"),
      xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf
    ),
    alpha = 0.1,
    fill = "orange"
  ) +
  geom_vline(aes(xintercept = as.Date("2019-03-13")), lty = "dashed", color = "yellow") +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1, color = "gray") +
  geom_line(aes(x = new_date, y = percent_error, color = model, group = model),
    data = all_residuals_lim %>% mutate(new_date = date + 7), lwd = 1, alpha = 0.8
  ) +
  scale_x_date(
    expand = c(0, 0), date_breaks = "3 months", date_labels = "%b %y",
    limits = c(as.Date("2019-01-22"), as.Date("2022-05-16"))
  ) +
  theme_bw(base_size = 16) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.5, 0.5)) +
  theme(legend.position = "bottom") +
  scale_color_manual(
    name = "Model",
    values = c("#ABDDA4", "#FDAE61", "#2B83BA", "navy", "purple", "#D7191C"),
    breaks = c("Adenovirus Rt", "AR", "AR + Mobility", "AR + Climate", "AR + Mobility + Climate", "AR + Mobility + Climate + SARS-CoV-2 Rt")
  ) +
  xlab("Date") +
  ylab("Percentage Error")
p

com <- plot_grid(r, p, nrow = 2, labels = "AUTO", align = "v")
com
save_plot(com, filename = "figures/figure_s21_adeno_predictive_model_fit_and_error_updated.png", base_width = 12, base_height = 10)

# ##### stay at home and initial rebound
accuracy_table <-
  model_predict_df %>%
  filter(date > as.Date("2020-02-28") & date < as.Date("2020-06-15")) %>%
  group_by(model) %>%
  summarize(
    rmse = rmse(predicted = predict_7d_ahead, actual = data_7d_ahead),
    mae = mae(predicted = predict_7d_ahead, actual = data_7d_ahead)
  )

baseline_rmse <- accuracy_table$rmse[1]
baseline_mae <- accuracy_table$mae[1]
accuracy_table$rmse[3] # 0.021626

accuracy_table %>%
  mutate(
    rmse_perc_diff_from_AR = 100 * (rmse - baseline_rmse) / baseline_rmse,
    mae_perc_diff_from_AR = 100 * (mae - baseline_mae) / baseline_mae
  ) %>%
  arrange(model)
# model                                     rmse    mae rmse_perc_diff_from_AR mae_perc_diff_from_AR
# 1 AR                                      0.0193 0.0141                  0                      0
# 2 AR + Climate                            0.0192 0.0138                 -0.535                 -2.20
# 3 AR + Mobility                           0.0216 0.0168                 12.2                   19.5
# 4 AR + Mobility + Climate                 0.0236 0.0181                 22.7                   28.5
# 5 AR + Mobility + Climate + SARS-CoV-2 Rt 0.0252 0.0194                 30.7                   37.5

# change in rmse between whole study period vs SAH
((0.04018487 - 0.021626) / 0.04018487) * 100 # 46.18%
