########################################################################################
## Moving window cross correlations between Rt and mobility, averaged by calendar month
## Figure S6: February 2019 snowstorm
########################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(gmodels) # for estimating CIs
library(ggplot2)
library(cowplot)

source("utils.R")

## pick lambda value so that edges of time windows have ~50% weight of midpoint
# exp_w <- exp_decay_func(30, l=0.05)
# exp_plot <- ggplot() +
#   geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
# exp_plot

########################################################################################
## Import data
########################################################################################

combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()

###########################################
## 2019 Snowstorm Moving Window
###########################################

## function for estimating monthly average of correlation coefficients and optimal lags
path_ccf_mv_window_snowstorm <- function(pathogen = "h1n1", l = 0.05) {
  metrics <- c(
    "within_neighborhood_movement",
    "within_city_movement",
    "within_state_movement",
    "out_of_state_movement",
    "performing_arts_or_sports_events",
    # "fb_leaving_home_custom", #fb % leaving home not available
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

  min_date <- combined %>%
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

  window <- 15 # window midpoint

  days <- unique(path_df_part1$date)[window:(length(unique(path_df_part1$date)) - window)]
  names(path_df_part1)[2] <- pathogen

  path_df_part1 <- path_df_part1 %>%
    dplyr::select(date, all_of(pathogen), all_of(metrics)) %>%
    na.omit()

  days <- unique(path_df_part1$date)[window:(length(unique(path_df_part1$date)) - window)]
  # length(days)
  names(path_df_part1)[2] <- pathogen

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

    # pearson correlations
    # res <- manual_ccf(a = t1, b = t2, lag.max = 21, w = exp_w)
    # # res <- res %>% filter(lag <1)
    # max_ccf_lag <- res[which.max(abs(res$cor)), ]$lag
    # max_ccf <- res[which.max(abs(res$cor)), ]$cor

    # spearman correlations
    t1_rank <- rank(t1)
    t2_rank <- rank(t2)
    res_sp <- manual_ccf(a = t1_rank, b = t2_rank, lag.max = 21, w = exp_w)
    # res_sp <- res_sp %>% filter(lag < 1)
    max_ccf_lag_sp <- res_sp[which.max(abs(res_sp$cor)), ]$lag
    max_ccf_sp <- res_sp[which.max(abs(res_sp$cor)), ]$cor

    actual_data <- data.frame(
      mobility_metric = as.character(X),
      start_week = as.character(Y),
      # obs_max_ccf_lag = as.numeric(max_ccf_lag),
      # obs_max_ccf = as.numeric(max_ccf),
      obs_max_ccf_lag_sp = as.numeric(max_ccf_lag_sp), # optimal lag
      obs_max_ccf_sp = as.numeric(max_ccf_sp), # correlation coefficient at optimal lag
      lambda_par = l
    )
    return(actual_data)
  }
  # y("amusement_and_recreation",days[1])
  tmp <- as.data.frame(t(df))
  mc.cores <- parallelly::availableCores() - 1
  result <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], SIMPLIFY = F, mc.cores = mc.cores)

  actual_data <- do.call(rbind.data.frame, result)
  actual_data$start_week <- as.Date(actual_data$start_week)
  # head(actual_data)
  # unique(actual_data$mobility_metric)

  actual_data$start_week <- as.Date(actual_data$start_week)
  actual_data$month <- zoo::as.yearmon(actual_data$start_week, "%Y %m")
  actual_data2 <- actual_data

  df <- actual_data2 %>%
    group_by(mobility_metric, month, lambda_par) %>%
    summarize(
      # mean_cor = mean(obs_max_ccf),
      # mean_lag = mean(obs_max_ccf_lag),
      mean_cor_sp = mean(obs_max_ccf_sp),
      mean_lag_sp = mean(obs_max_ccf_lag_sp),
    ) %>%
    arrange(month, -abs(mean_cor_sp)) %>%
    ungroup()

  actual_data3 <- left_join(actual_data2, df[, c("mobility_metric", "month", "mean_cor_sp", "lambda_par")],
    by = c("mobility_metric", "month", "lambda_par")
  )
  actual_data3$mobility_metric <- as.factor(actual_data3$mobility_metric)
  # levels(actual_data3$mobility_metric)

  actual_data3$mobility_metric <- factor(actual_data3$mobility_metric,
    levels = c(
      "within_neighborhood_movement",
      "within_city_movement",
      "within_state_movement",
      "out_of_state_movement",
      "full_service_restaurants",
      "groceries_and_pharmacies",
      "transit",
      "performing_arts_or_sports_events",
      "religious_orgs",
      "child_day_care",
      "elementary_and_secondary_schools",
      "colleges"
    )
  )

  levels(actual_data3$mobility_metric)
  moblabs <- c(
    "within-neighborhood movement",
    "between-neighborhood movement",
    "influx visitors other WA counties",
    "influx out-of-state visitors",
    "restaurants",
    "groceries and pharmacies",
    "transit",
    "performing arts or sports events",
    "religious organizations",
    "child daycare",
    "elementary and high schools",
    "colleges"
  )

  title <- actual_data3 %>%
    mutate(full_name = case_when(
      pathogen == "h1n1" ~ "Influenza A/H1N1",
      pathogen == "h3n2" ~ "Influenza A/H3N2",
      pathogen == "ivb" ~ "Influenza B",
      pathogen == "rsv_b" ~ "RSV B",
      pathogen == "rsv_a" ~ "RSV A",
      pathogen == "covid" ~ "COVID-19",
      pathogen == "hpiv_1_2" ~ "HPIV 1 + HPIV 2",
      pathogen == "hpiv_3_4" ~ "HPIV 3 + HPIV 4",
      pathogen == "hmpv" ~ "Human Metapneumovirus",
      pathogen == "scov_229E_OC43" ~ "HCoV 229E + HCoV OC43",
      pathogen == "scov_HKU1_NL63" ~ "HCoV HKU1 + HCoV NL63",
      pathogen == "adeno" ~ "Adenovirus",
      pathogen == "rhino" ~ "Rhinovirus",
      pathogen == "entero" ~ "Enterovirus"
    )) %>%
    pull(full_name) %>%
    unique()


  actual_data3$month <- as.factor(actual_data3$month)
  # unique(actual_data3$month)
  actual_data4 <- actual_data3 %>% filter(!(month %in% c("Apr 2019", "May 2019")))

  # mobility metric colors
  color_vec <- c(
    "blue",
    "#843C39",
    "#DE9ED6",
    "darkgreen",
    # "#c7eae5",
    "#f46d43",
    "#9C9EDE",
    "#E7969C",
    "#A55194",
    "#6DDE88",
    "#5254A3",
    "#BD9E39",
    "#637939",
    "skyblue"
  )


  sum_df <- actual_data3 %>%
    group_by(month, mobility_metric, lambda_par) %>%
    summarise(
      ## pearson correlation lag
      # obs_max_ccf_lag.mean = gmodels::ci(obs_max_ccf_lag,confidence = 0.95,alpha = 0.05, na.rm = T)[1],
      # obs_max_ccf_lag.lowCI = gmodels::ci(obs_max_ccf_lag,confidence = 0.95,alpha = 0.05, na.rm = T)[2],
      # obs_max_ccf_lag.hiCI = gmodels::ci(obs_max_ccf_lag,confidence = 0.95, alpha = 0.05, na.rm = T)[3],

      ## spearman correlation lag
      obs_max_ccf_lag_sp.mean = gmodels::ci(obs_max_ccf_lag_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[1],
      obs_max_ccf_lag_sp.lowCI = gmodels::ci(obs_max_ccf_lag_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[2],
      obs_max_ccf_lag_sp.hiCI = gmodels::ci(obs_max_ccf_lag_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[3],

      # pearson correlation coefficient
      # mean_cor = gmodels::ci(obs_max_ccf, confidence = 0.95, alpha = 0.05, na.rm = T)[1],
      # mean_cor.lowCI = gmodels::ci(obs_max_ccf, confidence = 0.95, alpha = 0.05, na.rm = T)[2],
      # mean_cor.hiCI = gmodels::ci(obs_max_ccf, confidence = 0.95, alpha = 0.05, na.rm = T)[3],

      # spearman correlation coefficient
      mean_cor_sp = gmodels::ci(obs_max_ccf_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[1],
      mean_cor_sp.lowCI = gmodels::ci(obs_max_ccf_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[2],
      mean_cor_sp.hiCI = gmodels::ci(obs_max_ccf_sp, confidence = 0.95, alpha = 0.05, na.rm = T)[3]
    ) %>%
    ungroup()

  # range(sum_df$obs_max_ccf_lag.mean)
  # range(sum_df$obs_max_ccf_lag_sp.mean)

  ## pearson correlations
  # q <- ggplot() +
  #   geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
  #   geom_errorbar(
  #     data = sum_df,
  #     aes(
  #       y = mean_cor,
  #       xmin = obs_max_ccf_lag.lowCI, xmax = obs_max_ccf_lag.hiCI,
  #       color = mobility_metric
  #     ), width = 0.08, lwd = 1
  #   ) +
  #   geom_point(
  #     data = sum_df,
  #     aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21
  #   ) +
  #   geom_hline(yintercept = 0, lty = "dashed") +
  #   geom_vline(xintercept = 0, lty = "dashed") +
  #   facet_wrap(~ as.factor(month), nrow = 1) +
  #   theme_bw(base_size = 16) +
  #   ylab("Cross-Correlation Coefficient") +
  #   xlab("Temporal Lag (days)") +
  #   scale_color_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #   scale_fill_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #   theme(
  #     legend.position = "bottom", legend.text = element_text(size = 14),
  #     axis.text.y = element_text(size = 14), strip.text = element_text(size = 16),
  #     title = element_text(size = 20), strip.background = element_blank()
  #   ) +
  #   ggtitle(title) +
  #   # scale_x_continuous(breaks = c(-10,-5, 0), limits = c(-12, 3))
  #   scale_x_continuous(breaks = c(-14, -7, 0, 7, 14),labels=c(-14, -7, 0, 7, 14), limits = c(-19, 19))
  # print(q)

  # spearman correlations
  r <- ggplot() +
    geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
    geom_errorbar(
      data = sum_df,
      aes(
        y = mean_cor_sp,
        xmin = obs_max_ccf_lag_sp.lowCI, xmax = obs_max_ccf_lag_sp.hiCI,
        color = mobility_metric
      ), width = 0.08, lwd = 1
    ) +
    geom_point(
      data = sum_df,
      aes(x = obs_max_ccf_lag_sp.mean, y = mean_cor_sp, fill = mobility_metric), size = 5, pch = 21
    ) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_vline(xintercept = 0, lty = "dashed") +
    facet_wrap(~ as.factor(month), nrow = 1) +
    theme_bw(base_size = 16) +
    ylab("Cross-Correlation Coefficient") +
    xlab("Temporal Lag (days)") +
    scale_color_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
    scale_fill_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 14),
      axis.text.y = element_text(size = 14), strip.text = element_text(size = 16),
      title = element_text(size = 20), strip.background = element_blank()
    ) +
    ggtitle(title) +
    scale_x_continuous(breaks = c(-14, -7, 0, 7, 14), labels = c(-14, -7, 0, 7, 14), limits = c(-19, 19))

  print(r)

  return(list(r, actual_data4, sum_df))
}
h1n1_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "h1n1", l = 0.05)
h1n1_mv_win_plot[[1]]

h3n2_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "h3n2", l = 0.05)
h3n2_mv_win_plot[[1]]

rsvb_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rsv_b", l = 0.05)
rsvb_mv_win_plot[[1]]

rsva_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rsv_a", l = 0.05)
rsva_mv_win_plot[[1]]

hpiv_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "hpiv_3_4", l = 0.05)
hpiv_mv_win_plot[[1]]

hmpv_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "hmpv", l = 0.05)
hmpv_mv_win_plot[[1]]
## no positive correlations

rhino_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rhino", l = 0.05)
rhino_mv_win_plot[[1]]
## weak correlations

entero_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "entero", l = 0.05)
entero_mv_win_plot[[1]]

adeno_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "adeno", l = 0.05)
adeno_mv_win_plot[[1]]

scov_229E_OC43_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "scov_229E_OC43", l = 0.05)
scov_229E_OC43_mv_win_plot[[1]]

scov_HKU1_NL63_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "scov_HKU1_NL63", l = 0.05)
scov_HKU1_NL63_mv_win_plot[[1]]
## no positive correlations

com <- plot_grid(
  adeno_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rsva_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rsvb_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  entero_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  h1n1_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  h3n2_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  scov_229E_OC43_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  scov_HKU1_NL63_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rhino_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  hmpv_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  nrow = 5, ncol = 2, byrow = F
)
com
leg <- get_legend(h1n1_mv_win_plot[[1]] + theme(legend.direction = "horizontal"))
com2 <- plot_grid(com, leg, nrow = 2, rel_heights = c(5, 0.3))
com2
save_plot(com2, filename = "figures/fig_s6_CCF_endemic_rt_and_mobility_monthly_snowstorm_spearman_exp_decay.png", base_width = 24, base_height = 16)
