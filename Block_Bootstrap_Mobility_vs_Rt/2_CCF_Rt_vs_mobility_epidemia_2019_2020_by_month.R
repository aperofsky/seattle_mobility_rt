########################################################################################
## Moving window cross correlations between Rt and mobility, averaged by calendar month
## Figure 4. Sept 2019 - May 2020
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
## Import data
########################################################################################

combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

# 2 week moving average
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

## RSV B 2021
rt_df_wide %>%
  filter(date > as.Date("2021-01-01")) %>%
  dplyr::select(date, rsv_b) %>%
  na.omit() %>%
  slice_min(date)
## first date: 2021-05-08

## HPIV 2021
rt_df_wide %>%
  filter(date > as.Date("2021-01-01")) %>%
  dplyr::select(date, hpiv_3_4) %>%
  na.omit() %>%
  slice_min(date)
## first date: 2021-02-07

rt_df_wide %>%
  filter(date < as.Date("2021-05-03")) %>%
  dplyr::select(date, rsv_b) %>%
  na.omit() %>%
  slice_max(date)
## last date rsv b 2019-2020: 2020-05-07

rt_df_wide %>%
  filter(date < as.Date("2021-06-20")) %>%
  dplyr::select(date, rsv_a) %>%
  na.omit() %>%
  slice_max(date)
# last date rsv a 2019-2020: 2020-04-22

########################################################################################
## Moving Window CCF
########################################################################################
path_ccf_mv_window <- function(pathogen = "covid", l = 0.07) {
  combined2 <- combined
  metrics <- c(
    "within_neighborhood_movement",
    "within_city_movement",
    "within_state_movement",
    "out_of_state_movement",
    "fb_leaving_home_custom",
    "full_service_restaurants",
    "groceries_and_pharmacies",
    "transit",
    "performing_arts_or_sports_events",
    "religious_orgs",
    "child_day_care",
    "elementary_and_secondary_schools",
    "colleges"
  )

  combined2$fb_leaving_home_custom <- 100 - combined2$fb_stay_put_custom

  combined2[all_of(pathogen)][is.na(combined2[all_of(pathogen)])] <- 0

  min_date <- as.Date("2019-06-01")
  min_date <- if_else(pathogen == "covid", as.Date("2019-11-01"), min_date)

  last_date <- combined2 %>%
    filter(eval(parse(text = pathogen)) > 0) %>%
    slice_max(epi_date) %>%
    pull(epi_date) %>%
    max()
  last_date

  max_date <- as.Date("2020-08-01")
  max_date <- if_else(pathogen %in% c("rhino", "adeno"), as.Date("2022-05-15"), max_date)
  max_date <- if_else(pathogen %in% c("covid"), as.Date("2021-05-05"), max_date)

  comb_weekly <- combined2 %>%
    filter(epi_date >= min_date & epi_date < max_date) %>%
    group_by(epi_date) %>%
    summarize_at(c(all_of(pathogen), metrics), ~ mean(.x, na.rm = T))
  comb_weekly[all_of(pathogen)][is.na(comb_weekly[all_of(pathogen)])] <- 0

  df1 <- comb_weekly
  window <- 10
  weeks <- unique(df1$epi_date)[window:(length(unique(df1$epi_date)) - window)]

  df <- expand.grid(X = metrics, Y = weeks) %>% arrange(X)
  y <- function(X, Y, l = l) {
    df <- df1[, c("epi_date", as.character(X), pathogen)] %>%
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
  # y("amusement_and_recreation",days[1])
  tmp <- as.data.frame(t(df))
  result <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, SIMPLIFY = F, mc.cores = parallel::detectCores() - 1)

  actual_data <- do.call(rbind.data.frame, result)
  actual_data$start_week <- as.Date(actual_data$start_week)
  range(actual_data$start_week)
  actual_data$start_week <- as.Date(actual_data$start_week)
  actual_data$month <- zoo::as.yearmon(actual_data$start_week, "%Y %m")
  actual_data2 <- actual_data

  df <- actual_data2 %>%
    group_by(mobility_metric, month) %>%
    summarize(
      mean_cor = mean(obs_max_ccf),
      mean_lag = mean(obs_max_ccf_lag)
    ) %>%
    arrange(month, -abs(mean_cor)) %>%
    ungroup()

  actual_data3 <- left_join(actual_data2, df[, c("mobility_metric", "month", "mean_cor")],
    by = c("mobility_metric", "month")
  )
  actual_data3$mobility_metric <- as.factor(actual_data3$mobility_metric)

  actual_data3$mobility_metric <- factor(actual_data3$mobility_metric,
    levels = c(
      "within_neighborhood_movement",
      "within_city_movement",
      "within_state_movement",
      "out_of_state_movement",
      "fb_leaving_home_custom",
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
      pathogen == "rhino" ~ "Rhinovirus"
    )) %>%
    pull(full_name) %>%
    unique()

  actual_data3$month <- as.factor(actual_data3$month)

  color_vec <- c(
    "blue", "#843C39", "#DE9ED6", "darkgreen","#c7eae5",
    "#f46d43", "#9C9EDE", "#E7969C", "#A55194", "#6DDE88", "#5254A3", "#BD9E39", "#637939", "skyblue"
  )

    moblabs <- c(
    "within-neighborhood movement",
    "between-neighborhood movement",
    "influx visitors other WA counties",
    "influx out-of-state visitors",
    "% devices leaving home",
    "restaurants",
    "groceries and pharmacies",
    "transit",
    "performing arts or sports events",
    "religious organizations",
    "child daycare",
    "elementary and high schools",
    "colleges"
  )

  output <- {
    if (pathogen %in% c("rhino", "adeno")) {
      actual_data4 <- actual_data3 %>%
        filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
      sum_df <- actual_data4 %>%
        filter(start_week < as.Date("2021-06-01")) %>%
        group_by(month, mobility_metric, mean_cor) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
          axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
          title = element_text(size = 16), strip.background = element_blank()
        ) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4)) +
        ggtitle(title) +
        ylim(-1, 1) +
        xlim(-4, 4)
      p
    }

    if (!(pathogen %in% c("rhino", "covid", "adeno"))) {
      actual_data4 <- actual_data3 %>%
        filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))

      sum_df <- actual_data4 %>%
        group_by(month, mobility_metric, mean_cor) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
          axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
          title = element_text(size = 16), strip.background = element_blank()
        ) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4)) +
        ggtitle(title) +
        ylim(-1, 1) +
        xlim(-4, 4)
    }

    ## SC2 has shorter time window of circulation, compared to endemic pathogens
    if (pathogen %in% c("covid")) {
      actual_data4 <- actual_data3 %>% filter(month %in% c("Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
      sum_df <- actual_data4 %>%
        group_by(month, mobility_metric, mean_cor) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
          axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
          title = element_text(size = 16), strip.background = element_blank()
        ) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4)) +
        ggtitle(title) +
        ylim(-1, 1) +
        xlim(-4, 4)
    }

    p
  }

  return(list(output, actual_data3, actual_data4, sum_df))
}
##########################################
## Influenza B
ivb_mv_win_plot <- path_ccf_mv_window(pathogen = "ivb", l = 0.07)
ivb_mv_win_plot[[1]]

##########################################
## Influenza A/H1N1
h1n1_mv_win_plot <- path_ccf_mv_window(pathogen = "h1n1", l = 0.07)
h1n1_mv_win_plot[[1]]

##########################################
## RSV B
rsvb_mv_win_plot <- path_ccf_mv_window(pathogen = "rsv_b", l = 0.07)
rsvb_mv_win_plot[[1]]

##########################################
## RSV A
rsva_mv_win_plot <- path_ccf_mv_window(pathogen = "rsv_a", l = 0.07)
rsva_mv_win_plot[[1]]

##########################################
## HPIV 3 +4
hpiv_mv_win_plot <- path_ccf_mv_window(pathogen = "hpiv_3_4", l = 0.07)
hpiv_mv_win_plot[[1]]

##########################################
## HPIV 1 +2: outbreak ended in mid-feb 2020 before start of covid npis (don't include in figure)
# hpiv_1_2_mv_win_plot = path_ccf_mv_window(pathogen = "hpiv_1_2",l=0.07)

##########################################
## hMPV
hmpv_mv_win_plot <- path_ccf_mv_window(pathogen = "hmpv", l = 0.07)
hmpv_mv_win_plot[[1]]

##########################################
## hCoV 229E + OC43
scov_229E_OC43_mv_win_plot <- path_ccf_mv_window(pathogen = "scov_229E_OC43", l = 0.07)
scov_229E_OC43_mv_win_plot[[1]]

##########################################
## hCoV HKU1 + NL63
scov_HKU1_NL63_mv_win_plot <- path_ccf_mv_window(pathogen = "scov_HKU1_NL63")
scov_HKU1_NL63_mv_win_plot[[1]]

##########################################
## SC2
covid_mv_win_plot <- path_ccf_mv_window(pathogen = "covid")
covid_mv_win_plot[[1]]

##########################################
## rhinovirus
rhino_mv_win_plot <- path_ccf_mv_window(pathogen = "rhino")
rhino_mv_win_plot[[1]]

##########################################
## adenovirus
adeno_mv_win_plot <- path_ccf_mv_window(pathogen = "adeno")
adeno_mv_win_plot[[1]]

##########################################
## Figure 4: Combine plots for all pathogens
##########################################
com <- plot_grid(
  h1n1_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  ivb_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  rsva_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  rsvb_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  hmpv_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  scov_229E_OC43_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  scov_HKU1_NL63_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  hpiv_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  rhino_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  adeno_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  nrow = 5, ncol = 2, byrow = F
)
com
leg <- get_legend(h1n1_mv_win_plot[[1]] + theme(legend.direction = "horizontal"))
com2 <- plot_grid(com, leg, nrow = 2, rel_heights = c(5, 0.3))
com2

bottom <- plot_grid(leg, covid_mv_win_plot[[1]] + ylab("Cross-Correlation") + theme(legend.position = "none") + ggtitle("SARS-CoV-2") + xlab("Temporal Lag (weeks)"),
  nrow = 1, rel_heights = c(1, 1),
  rel_widths = c(0.71, 0.29)
)
com2 <- plot_grid(com, bottom, nrow = 2, rel_heights = c(5, 1))
com2
save_plot(com2, filename = "figures/fig_4_CCF_pathogen_rt_and_mobility_by_month_2019_2020.png", base_width = 24, base_height = 16)
