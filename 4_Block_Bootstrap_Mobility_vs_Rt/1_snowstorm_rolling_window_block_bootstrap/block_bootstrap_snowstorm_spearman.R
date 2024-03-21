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
source("utils.R")

########################################################################################
# import data
########################################################################################

combined <- read_rds("combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()

########################################################################################
### block bootstrap
########################################################################################

## pick lambda value so that edges of time windows have ~50% weight of midpoint
path_ccf_mv_window_snowstorm <- function(pathogen = "h1n1", l = 0.05) {
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
    dplyr::select(date, all_of(pathogen)) %>%
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

    t1_rank <- rank(t1)
    t2_rank <- rank(t2)

    res <- manual_ccf(a = t1_rank, b = t2_rank, lag.max = 21, w = exp_w)
    res <- res %>% filter(lag < 1)
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
  ## constrain to lags <1
  extract_ccf <- function(cv) {
    cor <- cv$cor
    lag <- cv$lag
    res <- data.frame(cor, lag)
    res <- res %>% filter(lag < 1)
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
  y <- function(X = metrics, Y = days, perms = 1000, l = 0.05) {
    df_rt <- path_df_part1[, c("date", pathogen)] %>%
      dplyr::filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window)) %>%
      as.data.frame()

    mob_df <- path_df_part1[, c("date", as.character(X))] %>%
      filter(date >= lubridate::as_date(Y) - (window) & date <= lubridate::as_date(Y) + (window)) %>%
      as.data.frame()

    t1 <- mob_df[, 2] # mobility time series
    t2 <- df_rt[, 2] # rt time series
    n <- length(t1)
    exp_w <- exp_decay_func(n, l)

    t1_rank <- rank(t1)
    t2_rank <- rank(t2)

    set.seed(123)
    ts_boot <- tseries::tsbootstrap(x = t1_rank, nb = perms, statistic = NULL, type = "block", b = 14) # shuffle mobility time series
    lst <- lapply(seq_len(ncol(ts_boot)), function(i) ts_boot[, i])
    list_df <- lapply(lst, function(x) data.frame(date = mob_df[, 1], mob = x))

    filtered_df <- lapply(list_df, function(x) x %>% filter_epi(Y = Y, window = window))
    models <- lapply(filtered_df, function(x) manual_ccf(a = x[, 2], b = t2_rank, lag.max = 21, w = exp_w))
    output_list <- lapply(models, FUN = extract_ccf)
    output_df <- do.call(rbind.data.frame, output_list)
    return_df <- data.frame(mobility_metric = X, start_week = as.character(Y), output_df)
    return(return_df)
  }
  tmp <- as.data.frame(t(df))
  mc.cores <- parallelly::availableCores()
  system.time(result_list <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, perms = 1000, SIMPLIFY = F, mc.cores = mc.cores))
  result_df <- do.call(rbind.data.frame, result_list)

  result_df$start_week <- as.Date(result_df$start_week)
  null_and_obs_ccf <- left_join(result_df, actual_data, by = c("mobility_metric", "start_week"))
  null_and_obs_ccf$max_ccf <- as.numeric(null_and_obs_ccf$max_ccf)
  null_and_obs_ccf$obs_max_ccf <- as.numeric(null_and_obs_ccf$obs_max_ccf)

  # 95% CIs
  # d2 <- null_and_obs_ccf %>%
  #   group_by(mobility_metric, start_week) %>%
  #   summarize(lower = quantile(max_ccf, probs = .025),
  #             upper = quantile(max_ccf, probs = .975))%>%
  #   ungroup

  # 90% CIs
  d2 <- null_and_obs_ccf %>%
    group_by(mobility_metric, start_week) %>%
    summarize(
      lower = quantile(max_ccf, probs = .05),
      upper = quantile(max_ccf, probs = .95)
    ) %>%
    ungroup()

  output_df <- left_join(d2, actual_data, by = c("mobility_metric", "start_week")) %>%
    mutate(sig = ifelse(obs_max_ccf > upper | obs_max_ccf < lower, "yes", "no"))

  null_and_obs_ccf$mobility_metric <- as.factor(null_and_obs_ccf$mobility_metric)

  output_df$mobility_metric <- as.factor(output_df$mobility_metric)
  d2$mobility_metric <- as.factor(d2$mobility_metric)
  actual_data$mobility_metric <- as.factor(actual_data$mobility_metric)

  actual_data_and_perm <- left_join(actual_data, output_df %>% dplyr::select(mobility_metric, start_week, sig), by = c("mobility_metric", "start_week"))
  actual_data_and_perm$pathogen <- pathogen
  names(actual_data_and_perm)[2] <- "start_date"
  return(actual_data_and_perm)
}
h1n1_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "h1n1", l = 0.05)
save(h1n1_block_bs, file = "//data//perofskyamc//h1n1_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

h3n2_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "h3n2", l = 0.05)
save(h3n2_block_bs, file = "//data//perofskyamc//h3n2_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

rsvb_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rsv_b", l = 0.05)
save(rsvb_block_bs, file = "//data//perofskyamc//rsv_b_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

rsva_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rsv_a", l = 0.05)
save(rsva_block_bs, file = "//data//perofskyamc//rsv_a_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

hpiv_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "hpiv_3_4", l = 0.05)
save(hpiv_block_bs, file = "//data//perofskyamc//hpiv_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

hmpv_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "hmpv", l = 0.05)
save(hmpv_block_bs, file = "//data//perofskyamc//hmpv_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

rhino_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "rhino", l = 0.05)
save(rhino_block_bs, file = "//data//perofskyamc//rhino_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

entero_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "entero", l = 0.05)
save(entero_block_bs, file = "//data//perofskyamc//entero_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

adeno_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "adeno", l = 0.05)
save(adeno_block_bs, file = "//data//perofskyamc//adeno_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

scov_229E_OC43_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "scov_229E_OC43", l = 0.05)
save(scov_229E_OC43_block_bs, file = "//data//perofskyamc//scov_229E_OC43_snowstorm_sliding_window_actual_and_null_output_spearman.RData")

scov_HKU1_NL63_block_bs <- path_ccf_mv_window_snowstorm(pathogen = "scov_HKU1_NL63", l = 0.05)
save(scov_HKU1_NL63_block_bs, file = "//data//perofskyamc//scov_HKU1_NL63_snowstorm_sliding_window_actual_and_null_output_spearman.RData")
