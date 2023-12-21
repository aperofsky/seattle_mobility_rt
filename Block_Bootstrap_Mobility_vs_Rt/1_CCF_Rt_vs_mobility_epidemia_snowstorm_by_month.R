########################################################################################
## Moving window cross correlations between Rt and mobility, averaged by calendar month
## Figure S6: February 2019 snowstorm
########################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(tseries)
library(gmodels)

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
## Import data
########################################################################################

combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

rt_df <- read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

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

###########################################
## 2019 Snowstorm Moving Window
###########################################
path_ccf_mv_window_snowstorm <- function(pathogen = "h1n1", l = 0.01) {
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
    dplyr::select(date, all_of(pathogen), all_of(metrics)) %>%
    na.omit()

  days <- unique(path_df_part1$date)[window:(length(unique(path_df_part1$date)) - window)]
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

    exp_plot <- ggplot() +
      geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
    exp_plot

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
  head(actual_data)
  unique(actual_data$mobility_metric)

  actual_data$start_week <- as.Date(actual_data$start_week)
  actual_data$month <- zoo::as.yearmon(actual_data$start_week, "%Y %m")
  actual_data2 <- actual_data

  df <- actual_data2 %>%
    group_by(mobility_metric, month, lambda_par) %>%
    summarize(
      mean_cor = mean(obs_max_ccf),
      mean_lag = mean(obs_max_ccf_lag)
    ) %>%
    arrange(month, -abs(mean_cor)) %>%
    ungroup()

  actual_data3 <- left_join(actual_data2, df[, c("mobility_metric", "month", "mean_cor", "lambda_par")],
    by = c("mobility_metric", "month", "lambda_par")
  )
  actual_data3$mobility_metric <- as.factor(actual_data3$mobility_metric)
  levels(actual_data3$mobility_metric)

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
      pathogen == "covid_wt" ~ "COVID-19",
      pathogen == "hpiv_1_2" ~ "HPIV 1 + HPIV 2",
      pathogen == "hpiv_3_4" ~ "HPIV 3 + HPIV 4",
      pathogen == "hmpv" ~ "Human Metapneumovirus",
      pathogen == "scov_229E_OC43" ~ "HCoV 229E + HCoV OC43",
      pathogen == "scov_HKU1_NL63" ~ "HCoV HKU1 + HCoV NL63",
      pathogen == "adeno" ~ "Adenovirus",
      pathogen == "rhino" ~ "Rhinovirus"
    )) %>%
    pull(full_name) %>%
    unique()


  actual_data3$month <- as.factor(actual_data3$month)
  unique(actual_data3$month)
  actual_data4 <- actual_data3 %>% filter(!(month %in% c("Apr 2019", "May 2019")))

  color_vec <- c(
    "blue", "#843C39", "#DE9ED6", "darkgreen",
    # "#c7eae5",
    "#f46d43", "#9C9EDE", "#E7969C", "#A55194", "#6DDE88", "#5254A3", "#BD9E39", "#637939", "skyblue"
  )


  sum_df <- actual_data3 %>%
    group_by(month, mobility_metric, lambda_par) %>%
    summarise(
      obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
      obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
      obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3],
      mean_cor = ci(obs_max_ccf, na.rm = T)[1],
      mean_cor.lowCI = ci(obs_max_ccf, na.rm = T)[2],
      mean_cor.hiCI = ci(obs_max_ccf, na.rm = T)[3]
    ) %>%
    ungroup()

  q <- ggplot() +
    geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
    geom_errorbar(
      data = sum_df,
      aes(
        y = mean_cor,
        xmin = obs_max_ccf_lag.lowCI, xmax = obs_max_ccf_lag.hiCI,
        color = mobility_metric
      ), width = 0.08, lwd = 1
    ) +
    geom_point(
      data = sum_df,
      aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21
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
    scale_x_continuous(breaks = c(-10, -5, 0, 5, 10), limits = c(-14, 14))

  print(q)


  return(list(q, actual_data4, sum_df))
}
h1n1_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "h1n1", l = 0.001)
h1n1_mv_win_plot[[1]]
h1n1_mv_win_plot[[3]] %>% filter(month == "Feb 2019")
h1n1_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

h3n2_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "h3n2", l = 0.001)
h3n2_mv_win_plot[[1]]
h3n2_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

rsvb_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rsv_b", l = 0.001)
rsvb_mv_win_plot[[1]]
rsvb_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

rsva_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rsv_a", l = 0.001)
rsva_mv_win_plot[[1]]
rsva_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

hpiv_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "hpiv_3_4", l = 0.001)
hpiv_mv_win_plot[[1]]
hpiv_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

hmpv_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "hmpv", l = 0.001)
hmpv_mv_win_plot[[1]]
hmpv_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)
## no positive correlations

rhino_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "rhino", l = 0.001)
rhino_mv_win_plot[[1]]
rhino_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)
## weak correlations

adeno_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "adeno", l = 0.001)
adeno_mv_win_plot[[1]]
adeno_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019") %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

scov_229E_OC43_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "scov_229E_OC43", l = 0.001)
scov_229E_OC43_mv_win_plot[[1]]
scov_229E_OC43_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)

scov_HKU1_NL63_mv_win_plot <- path_ccf_mv_window_snowstorm(pathogen = "scov_HKU1_NL63", l = 0.001)
scov_HKU1_NL63_mv_win_plot[[1]]
scov_HKU1_NL63_mv_win_plot[[3]] %>%
  filter(month == "Feb 2019" & obs_max_ccf_lag.mean < 1 & mean_cor > 0) %>%
  arrange(-mean_cor)
## no positive correlations

com <- plot_grid(
  adeno_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rsva_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rsvb_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  h1n1_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  h3n2_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  scov_229E_OC43_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  rhino_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  hmpv_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation"),
  nrow = 4, ncol = 2, byrow = F
)
com
leg <- get_legend(h1n1_mv_win_plot[[1]] + theme(legend.direction = "horizontal"))
com2 <- plot_grid(com, leg, nrow = 2, rel_heights = c(5, 0.3))
com2
save_plot(com2, filename = "figures/fig_s6_CCF_endemic_rt_and_mobility_monthly_snowstorm.png", base_width = 22, base_height = 14)
