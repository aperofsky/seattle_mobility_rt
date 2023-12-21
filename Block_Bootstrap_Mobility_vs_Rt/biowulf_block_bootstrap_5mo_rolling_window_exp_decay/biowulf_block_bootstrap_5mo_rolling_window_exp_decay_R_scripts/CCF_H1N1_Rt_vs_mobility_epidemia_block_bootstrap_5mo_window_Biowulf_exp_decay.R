setwd("//home//perofskyamc//SFS_Rt_Block_Bootstrap")
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(tseries)

########################################################################################
## exponential decay function
########################################################################################

exp_decay_func <- function(n, l) {
  # n <- window * 2

  # Define the midpoint index
  midpoint <- (n + 1) / 2

  # Define the decay parameter
  lambda <- l

  # Generate the exponential decay vector
  exponential_decay <- exp(-lambda * abs(1:n - midpoint))

  # Scale the vector to have maximum value at the midpoint
  exponential_decay <- exponential_decay / max(exponential_decay)

  # Print the resulting vector
  return(exponential_decay)
}

# pr = exp_decay_func(n=20,l=0.07)
# ggplot()+
#   geom_line(aes(x=seq(1:length(pr)),y=pr))

########################################################################################
# weighted ccf
########################################################################################

manual_ccf <- function(a, b, w, lag.max) {
  n <- length(a)

  # checking if the sum of the weights are 1, if not normalize them
  if (sum(w) != 1) {
    w <- w / sum(w)
  }
  # weighted mean of x and y
  ma <- weighted.mean(a, w)
  mb <- weighted.mean(b, w)

  # std deviation of weighted x and y
  c_0 <- sqrt(sum(w * (a - ma)^2 / n) * sum(w * (b - mb)^2 / n))
  y <- NULL
  for (t in -lag.max:lag.max) {
    if (t <= 0) {
      w_t <- w[1:(n + t)]
      c_t <- 1 / n * sum(w_t * (a[1:(n + t)] - ma) * (b[(1 - t):n] - mb))
    } else {
      w_t <- w[(1 + t):n]
      c_t <- 1 / n * sum(w_t * (a[(1 + t):n] - ma) * (b[1:(n - t)] - mb))
    }
    r_t <- c_t / c_0
    y <- rbind(y, r_t)
  }
  ccf_df <- y
  ccf_df <- data.frame(cor = y, lag = -lag.max:lag.max)
  return(ccf_df)
}

########################################################################################
## import data
########################################################################################

combined_mob <- read_rds("mobility_metrics_for_epidemia.rds") %>% as_tibble()

combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

rt_df <- read_csv("rt_all_pathogens_15day_mv_avg.csv")
# rt_df = read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

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
# combined = left_join(rt_df_wide,
#                      combined_mob,by="date")
###########################################
## Actual CCF Loop
###########################################

metrics <- c(
  "within_neighborhood_movement",
  "within_city_movement",
  "within_state_movement",
  "out_of_state_movement",
  "fb_leaving_home_custom",
  "full_service_restaurants",
  "groceries_and_pharmacies",
  "full_service_restaurants",
  "transit",
  "religious_orgs",
  "child_day_care",
  "elementary_and_secondary_schools",
  "colleges"
)

combined$fb_leaving_home_custom <- 100 - combined$fb_stay_put_custom
combined["h1n1"][is.na(combined["h1n1"])] <- 0

comb_weekly <- combined %>%
  mutate_at(vars(all_of(metrics)), ~ scale(.x) %>% as.vector()) %>%
  filter(epi_date >= as.Date("2019-06-01") & epi_date < as.Date("2020-08-01")) %>%
  group_by(epi_date) %>%
  summarize_at(c("h1n1", metrics), ~ mean(.x, na.rm = T))
comb_weekly["h1n1"][is.na(comb_weekly["h1n1"])] <- 0

df1 <- comb_weekly
window <- 10
weeks <- unique(df1$epi_date)[window:(length(unique(df1$epi_date)) - window)]

df <- expand.grid(X = metrics, Y = weeks) %>% arrange(X)

y <- function(X, Y, l = 0.07) {
  df <- df1[, c("epi_date", as.character(X), "h1n1")] %>%
    filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window))
  start <- min(df$epi_date)
  end <- max(df$epi_date)

  t1 <- df[, 2] %>% pull(1)
  t2 <- df[, 3] %>% pull(1)

  n <- length(t1)
  exp_w <- exp_decay_func(n, l)

  exp_plot <- ggplot() +
    geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
  exp_plot

  res <- manual_ccf(a = t1, b = t2, lag.max = 4, w = exp_w)
  # cv <- ccf(x = df[,2], y =df[,3], type = c("correlation"),plot = F,lag.max = 4)
  # cor = cv$acf[,,1]
  # lag = cv$lag[,,1]
  # res = data.frame(cor,lag)
  max_ccf_lag <- res[which.max(abs(res$cor)), ]$lag
  max_ccf <- res[which.max(abs(res$cor)), ]$cor
  actual_data <- data.frame(
    mobility_metric = as.character(X),
    start_week = as.character(Y),
    obs_max_ccf_lag = as.numeric(max_ccf_lag),
    obs_max_ccf = as.numeric(max_ccf)
  )
  return(actual_data)
}
# y("within_neighborhood_movement","2021-04-25")
tmp <- as.data.frame(t(df))
result <- mapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = 0.07, SIMPLIFY = F)

actual_data <- do.call(rbind.data.frame, result)
actual_data$start_week <- as.Date(actual_data$start_week)
###########################################
## Null Distribution CCFs
###########################################

extract_ccf <- function(cv) {
  # cor = cv$acf[,,1]
  # lag = cv$lag[,,1]
  # res = data.frame(cor,lag)
  cor <- cv$cor
  lag <- cv$lag
  res <- data.frame(cor, lag)
  max_ccf_lag <- res[which.max(abs(res$cor)), ]$lag
  max_ccf <- res[which.max(abs(res$cor)), ]$cor
  df <- data.frame(max_ccf_lag, max_ccf)
  return(df)
}


filter_epi <- function(x, Y, window) {
  df2 <- x %>% dplyr::filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window))
  return(df2)
}

df <- expand.grid(X = metrics, Y = weeks) %>% arrange(X)

y <- function(X = metrics, Y = weeks, perms = 1000, l = 0.07) {
  actual_hosp <- df1[, c("epi_date", "h1n1")] %>%
    filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window)) %>%
    as.data.frame()

  mob_df <- df1[, c("epi_date", X)] %>%
    filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window)) %>%
    as.data.frame()

  t1 <- mob_df[, 2]
  t2 <- actual_hosp[, 2]

  n <- length(t1)
  exp_w <- exp_decay_func(n, l)

  set.seed(123)
  ts_boot <- tseries::tsbootstrap(x = mob_df[, 2], nb = perms, statistic = NULL, type = "block", b = 2)
  l <- lapply(seq_len(ncol(ts_boot)), function(i) ts_boot[, i])
  list_df <- lapply(l, function(x) data.frame(epi_date = mob_df[, 1], mob = x))

  filtered_df <- lapply(list_df, function(x) x %>% filter_epi(Y = Y, window = window))
  models <- lapply(filtered_df, function(x) manual_ccf(a = x[, 2], b = actual_hosp[, 2], lag.max = 4, w = exp_w))
  output_list <- lapply(models, FUN = extract_ccf)
  output_df <- do.call(rbind.data.frame, output_list)
  return_df <- data.frame(mobility_metric = X, start_week = as.character(Y), output_df)
  return(return_df)
}
tmp <- as.data.frame(t(df))

mc.cores <- parallelly::availableCores()
system.time(result_list <- mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = 0.07, perms = 10, SIMPLIFY = F, mc.cores = mc.cores)) # analysis takes 2.5 hours
# analysis takes 1.7 hrs
#     user    system   elapsed
# 12085.438    74.857  6285.194
result_df <- do.call(rbind.data.frame, result_list)

result_df$start_week <- as.Date(result_df$start_week)
null_and_obs_ccf <- left_join(result_df, actual_data, by = c("mobility_metric", "start_week"))
null_and_obs_ccf$max_ccf <- as.numeric(null_and_obs_ccf$max_ccf)
null_and_obs_ccf$obs_max_ccf <- as.numeric(null_and_obs_ccf$obs_max_ccf)

# d2 <- null_and_obs_ccf %>%
#   group_by(mobility_metric, start_week) %>%
#   summarize(lower = quantile(abs(max_ccf), probs = .025),
#             upper = quantile(abs(max_ccf), probs = .975))%>%
#   ungroup

d2 <- null_and_obs_ccf %>%
  group_by(mobility_metric, start_week) %>%
  summarize(
    lower = quantile(abs(max_ccf), probs = .05),
    upper = quantile(abs(max_ccf), probs = .95)
  ) %>%
  ungroup()

output_df <- left_join(d2, actual_data, by = c("mobility_metric", "start_week")) %>%
  mutate(sig = ifelse(abs(obs_max_ccf) > upper | abs(obs_max_ccf) < lower, "yes", "no"))

null_and_obs_ccf$mobility_metric <- as.factor(null_and_obs_ccf$mobility_metric)

output_df$mobility_metric <- as.factor(output_df$mobility_metric)
d2$mobility_metric <- as.factor(d2$mobility_metric)

actual_data$mobility_metric <- as.factor(actual_data$mobility_metric)

actual_data_and_perm <- left_join(actual_data, output_df %>% dplyr::select(mobility_metric, start_week, sig), by = c("mobility_metric", "start_week"))
save(actual_data_and_perm, file = "//data//perofskyamc//h1n1_mobility_CCF_5mo_sliding_window_actual_and_null_output.RData")
