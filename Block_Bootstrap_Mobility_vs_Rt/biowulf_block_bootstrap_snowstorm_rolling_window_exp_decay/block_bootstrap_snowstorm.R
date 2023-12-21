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
# ed = exp_decay_func(150,0.01)
# ggplot()+geom_line(aes(x=seq(1:length(ed)),y=ed))

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
# import data
########################################################################################

combined_mob <- read_rds("mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)
combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

rt_df <- read_csv("rt_all_pathogens_15day_mv_avg.csv")

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

####################
### block bootstrap

path_ccf_mv_window_snowstorm <- function(pathogen = "h1n1", l = 0.001) {
  metrics <- c(
    "within_neighborhood_movement",
    "within_city_movement",
    "within_state_movement",
    "out_of_state_movement",
    "performing_arts_or_sports_events",
    # "fb_leaving_home_custom",
    "full_service_restaurants",
    "groceries_and_pharmacies",
    "full_service_restaurants",
    "transit",
    "religious_orgs",
    "child_day_care",
    "elementary_and_secondary_schools",
    "colleges"
  )

  path_df <- combined %>%
    dplyr::select(date, all_of(pathogen), all_of(metrics)) %>%
    complete(date = seq.Date(as.Date(min(date)), as.Date(max(date)), by = "day"))
  path_df[2][is.na(path_df[2])] <- 0

  path_df_scaled <- path_df %>%
    mutate_at(vars(all_of(metrics)), ~ scale(.x) %>% as.vector())

  min_date <- rt_df_wide %>%
    filter(date < as.Date("2019-06-01")) %>%
    dplyr::select(date, pathogen) %>%
    na.omit() %>%
    pull(date) %>%
    first()

  last_date <- combined %>%
    filter(date < as.Date("2019-06-01")) %>%
    filter(eval(parse(text = pathogen)) > 0) %>%
    slice_max(date) %>%
    pull(date)

  max_date <- last_date
  max_date

  path_df_part1 <- path_df_scaled %>%
    filter(date >= as.Date(min_date) & date <= as.Date(max_date)) %>%
    as_tibble()

  window <- 15

  days <- unique(path_df_part1$date)[window:(length(unique(path_df_part1$date)) - window)]
  names(path_df_part1)[2] <- pathogen

  path_df_part1 <- path_df_part1 %>%
    dplyr::select(date, pathogen, all_of(metrics)) %>%
    na.omit()

  df <- expand.grid(X = metrics, Y = days) %>% arrange(X)
  y <- function(X, Y) {
    df <- path_df_part1[, c("date", as.character(X), pathogen)] %>%
      dplyr::filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window))
    start <- min(df$date)
    end <- max(df$date)

    t1 <- df[, 2] %>% pull(1)
    t2 <- df[, 3] %>% pull(1)

    n <- length(t1)
    exp_w <- exp_decay_func(n, l)

    # exp_plot <- ggplot() +
    #   geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
    # exp_plot

    res <- manual_ccf(a = t1, b = t2, lag.max = 21, w = exp_w)
    max_ccf_lag <- res[which.max(abs(res$cor)), ]$lag
    max_ccf <- res[which.max(abs(res$cor)), ]$cor
    actual_data <- data.frame(
      mobility_metric = as.character(X),
      start_week = as.character(Y),
      obs_max_ccf_lag = as.numeric(max_ccf_lag),
      obs_max_ccf = as.numeric(max_ccf),
      lambda_par = l
    )
    return(actual_data)
  }
  # y("amusement_and_recreation",days[1])
  tmp <- as.data.frame(t(df))
  result <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], SIMPLIFY = F, mc.cores = parallel::detectCores() - 1)

  actual_data <- do.call(rbind.data.frame, result)
  actual_data$start_week <- as.Date(actual_data$start_week)
  actual_data$month <- zoo::as.yearmon(actual_data$start_week, "%Y %m")

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
    df2 <- x %>% dplyr::filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window))
    return(df2)
  }

  df <- expand.grid(X = metrics, Y = days) %>% arrange(X)

  y <- function(X = metrics, Y = days, perms = 1000, l = 0.1) {
    df_rt <- path_df_part1[, c("date", pathogen)] %>%
      dplyr::filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window)) %>%
      as.data.frame()

    mob_df <- path_df_part1[, c("date", as.character(X))] %>%
      filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window)) %>%
      as.data.frame()

    t1 <- mob_df[, 2]
    t2 <- df_rt[, 2]

    n <- length(t1)
    exp_w <- exp_decay_func(n, l)

    set.seed(123)
    ts_boot <- tseries::tsbootstrap(x = mob_df[, 2], nb = perms, statistic = NULL, type = "block", b = 14)
    l <- lapply(seq_len(ncol(ts_boot)), function(i) ts_boot[, i])
    list_df <- lapply(l, function(x) data.frame(date = mob_df[, 1], mob = x))

    filtered_df <- lapply(list_df, function(x) x %>% filter_epi(Y = Y, window = window))
    models <- lapply(filtered_df, function(x) manual_ccf(a = x[, 2], b = df_rt[, 2], lag.max = 28, w = exp_w))
    output_list <- lapply(models, FUN = extract_ccf)
    output_df <- do.call(rbind.data.frame, output_list)
    return_df <- data.frame(mobility_metric = X, start_week = as.character(Y), output_df)
    return(return_df)
  }
  tmp <- as.data.frame(t(df))
  mc.cores <- parallelly::availableCores()
  system.time(result_list <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, perms = 1000, SIMPLIFY = F, mc.cores = mc.cores)) # analysis takes 2.5 hours
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
  actual_data_and_perm$pathogen <- pathogen
  names(actual_data_and_perm)[2] <- "start_date"
  return(actual_data_and_perm)
}
h1n1_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "h1n1", l = 0.001)
save(h1n1_block_bs, file = "//data//perofskyamc//h1n1_snowstorm_sliding_window_actual_and_null_output.RData")

h3n2_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "h3n2", l = 0.001)
save(h3n2_block_bs, file = "//data//perofskyamc//h3n2_snowstorm_sliding_window_actual_and_null_output.RData")

rsvb_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rsv_b", l = 0.001)
save(rsvb_block_bs, file = "//data//perofskyamc//rsv_b_snowstorm_sliding_window_actual_and_null_output.RData")

rsva_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rsv_a", l = 0.001)
save(rsva_block_bs, file = "//data//perofskyamc//rsv_a_snowstorm_sliding_window_actual_and_null_output.RData")

hpiv_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "hpiv_3_4", l = 0.001)
save(hpiv_block_bs, file = "//data//perofskyamc//hpiv_snowstorm_sliding_window_actual_and_null_output.RData")

hmpv_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "hmpv", l = 0.001)
save(hmpv_block_bs, file = "//data//perofskyamc//hmpv_snowstorm_sliding_window_actual_and_null_output.RData")

rhino_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rhino", l = 0.001)
save(rhino_block_bs, file = "//data//perofskyamc//rhino_snowstorm_sliding_window_actual_and_null_output.RData")

adeno_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "adeno", l = 0.001)
save(adeno_block_bs, file = "//data//perofskyamc//adeno_snowstorm_sliding_window_actual_and_null_output.RData")

scov_229E_OC43_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "scov_229E_OC43", l = 0.001)
save(scov_229E_OC43_block_bs, file = "//data//perofskyamc//scov_229E_OC43_snowstorm_sliding_window_actual_and_null_output.RData")

scov_HKU1_NL63_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "scov_HKU1_NL63", l = 0.001)
save(scov_HKU1_NL63_block_bs, file = "//data//perofskyamc//scov_HKU1_NL63_snowstorm_sliding_window_actual_and_null_output.RData")
