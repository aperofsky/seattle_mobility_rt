########################################################################################
## Moving window cross correlations between Rt and mobility, averaged by calendar month
## Figure 4. Sept 2019 - May 2020 (select mobility indicators)
########################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(parallelly)
library(zoo)
library(lubridate)
library(tseries)
library(gmodels)
library(ggplot2)
library(cowplot)
library(patchwork)

source("utils.R")

## pick lambda value so that edges of time windows have ~50% weight of midpoint
# pr = exp_decay_func(n=20,l=0.077) #20 weeks
# ggplot()+
#   geom_line(aes(x=seq(1:length(pr)),y=pr))

########################################################################################
## Import data
########################################################################################
combined <- read_rds("3_Combine_Mobility_and_Rt_data/combined_rt_mobility_15day_mv_avg.rds") %>% as_tibble()

## RSV B 2021
combined %>%
  filter(date > as.Date("2021-01-01")) %>%
  dplyr::select(date, rsv_b) %>%
  na.omit() %>%
  slice_min(date)
## first date: 2021-05-08

## HPIV 2021
combined %>%
  filter(date > as.Date("2021-01-01")) %>%
  dplyr::select(date, hpiv_3_4) %>%
  na.omit() %>%
  slice_min(date)
## first date: 2021-02-07

combined %>%
  filter(date < as.Date("2021-05-03")) %>%
  dplyr::select(date, rsv_b) %>%
  na.omit() %>%
  slice_max(date)
## last date rsv b 2019-2020: 2020-05-07

combined %>%
  filter(date < as.Date("2021-06-20")) %>%
  dplyr::select(date, rsv_a) %>%
  na.omit() %>%
  slice_max(date)
# last date rsv a 2019-2020: 2020-04-22

combined %>% filter(covid > 0) # 2020-02-25

# limit covid Rt estimates to after 3/1/2020
combined <- combined %>% mutate(covid = if_else(date < as.Date("2020-03-01"), NA, covid))
combined %>% filter(covid > 0) # 2020-03-01

########################################################################################
## Moving Window CCF
########################################################################################
path_ccf_mv_window <- function(pathogen = "covid", l = 0.077) {
  combined2 <- combined

  # limit to metrics that correlate with endemic virus Rt in Fall 2019
  metrics <- c(
    # "within_neighborhood_movement",
    "within_city_movement",
    # "within_state_movement",
    "out_of_state_movement",
    "fb_leaving_home_custom",
    # "full_service_restaurants",
    # "groceries_and_pharmacies",
    # "transit",
    # "performing_arts_or_sports_events",
    "religious_orgs",
    "child_day_care",
    "elementary_and_secondary_schools",
    "colleges"
  )

  combined2$fb_leaving_home_custom <- 100 - combined2$fb_stay_put_custom

  combined2[pathogen][is.na(combined2[pathogen])] <- 0

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
  # length(unique(df1$epi_date))
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

    # exp_plot <- ggplot() +
    #   geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
    # exp_plot

    res <- manual_ccf(a = t1, b = t2, lag.max = 4, w = exp_w)
    res
    # res <- res %>% filter(lag < 1)
    max_ccf_lag <- res[which.max(abs(res$cor)), ]$lag
    max_ccf <- res[which.max(abs(res$cor)), ]$cor

    t1_rank <- rank(t1)
    t2_rank <- rank(t2)
    res_sp <- manual_ccf(a = t1_rank, b = t2_rank, lag.max = 4, w = exp_w)
    # res_sp <- res_sp %>% filter(lag < 1)
    max_ccf_lag_sp <- res_sp[which.max(abs(res_sp$cor)), ]$lag
    max_ccf_sp <- res_sp[which.max(abs(res_sp$cor)), ]$cor


    actual_data <- data.frame(
      mobility_metric = as.character(X),
      start_week = as.character(Y),
      # obs_max_ccf_lag = as.numeric(max_ccf_lag),
      # obs_max_ccf = as.numeric(max_ccf),
      obs_max_ccf_lag_sp = as.numeric(max_ccf_lag_sp),
      obs_max_ccf_sp = as.numeric(max_ccf_sp)
    )
    return(actual_data)
  }
  # y("amusement_and_recreation",days[1])
  tmp <- as.data.frame(t(df))
  mc.cores <- parallelly::availableCores() - 1
  result <- parallel::mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, SIMPLIFY = F, mc.cores = mc.cores)
  actual_data <- do.call(rbind.data.frame, result)
  actual_data$start_week <- as.Date(actual_data$start_week)
  # range(actual_data$start_week)
  actual_data$start_week <- as.Date(actual_data$start_week)
  actual_data$month <- zoo::as.yearmon(actual_data$start_week, "%Y %m")
  actual_data2 <- actual_data
  # names(actual_data2)

  df <- actual_data2 %>%
    group_by(mobility_metric, month) %>%
    summarize(
      # mean_cor = mean(obs_max_ccf),
      # mean_lag = mean(obs_max_ccf_lag),
      mean_cor_sp = mean(obs_max_ccf_sp),
      mean_lag_sp = mean(obs_max_ccf_lag_sp)
    ) %>%
    arrange(month, -abs(mean_cor_sp)) %>%
    ungroup()
  # head(df)

  actual_data3 <- left_join(actual_data2, df[, c("mobility_metric", "month", "mean_cor_sp")],
    by = c("mobility_metric", "month")
  )
  actual_data3$mobility_metric <- as.factor(actual_data3$mobility_metric)

  actual_data3$mobility_metric <- factor(actual_data3$mobility_metric,
    levels = c(
      # "within_neighborhood_movement",
      "within_city_movement",
      # "within_state_movement",
      "out_of_state_movement",
      "fb_leaving_home_custom",
      # "full_service_restaurants",
      # "groceries_and_pharmacies",
      # "transit",
      # "performing_arts_or_sports_events",
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
      pathogen == "rhino" ~ "Rhinovirus",
      pathogen == "entero" ~ "Enterovirus"
    )) %>%
    pull(full_name) %>%
    unique()

  actual_data3$month <- as.factor(actual_data3$month)

  # "within_neighborhood_movement",
  # "within_city_movement",
  # "within_state_movement",
  # "out_of_state_movement",
  # "fb_leaving_home_custom",
  # "full_service_restaurants",
  # "groceries_and_pharmacies",
  # "transit",
  # "performing_arts_or_sports_events",
  # "religious_orgs",
  # "child_day_care",
  # "elementary_and_secondary_schools",
  # "colleges"

  # mobility metric colors
  # color_vec <- c(
  #   "blue",
  #   "#843C39",
  #   "#DE9ED6",
  #   "darkgreen",
  #   "#c7eae5",
  #   "#f46d43",
  #   "#9C9EDE",
  #   "#E7969C",
  #   "#A55194",
  #   "#6DDE88",
  #   "#5254A3",
  #   "#BD9E39",
  #   "#637939",
  #   "skyblue"
  # )

  color_vec <- c(
    # "blue",
    "#843C39",
    # "#DE9ED6",
    "darkgreen",
    "#c7eae5",
    # "#f46d43",
    # "#9C9EDE",
    # "#E7969C",
    # "#A55194",
    "#6DDE88",
    "#5254A3",
    "#BD9E39",
    "#637939",
    "skyblue"
  )

  moblabs <- c(
    # "within-neighborhood movement",
    "between-neighborhood movement",
    # "influx visitors other WA counties",
    "influx out-of-state visitors",
    "% devices leaving home",
    # "restaurants",
    # "groceries and pharmacies",
    # "transit",
    # "performing arts or sports events",
    "religious organizations",
    "child daycare",
    "elementary and high schools",
    "colleges"
  )

  # pearson correlations
  # output <- {
  #   if (pathogen %in% c("rhino", "adeno")) {
  #     actual_data4 <- actual_data3 %>%
  #       filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
  #
  #     sum_df <- actual_data4 %>%
  #       filter(start_week < as.Date("2021-06-01")) %>%
  #       group_by(month, mobility_metric, mean_cor) %>%
  #       summarise(
  #         obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
  #         obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
  #         obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
  #       ) %>%
  #       ungroup()
  #     # filter(obs_max_ccf_lag.mean < 0.5)
  #
  #     p <- ggplot() +
  #       geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
  #       geom_point(
  #         data = sum_df,
  #         aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
  #       ) +
  #       geom_hline(yintercept = 0, lty = "dashed") +
  #       geom_vline(xintercept = 0, lty = "dashed") +
  #       facet_wrap(~ as.factor(month), nrow = 1) +
  #       theme_bw(base_size = 16) +
  #       ylab("Cross-Correlation Coefficient") +
  #       xlab("Temporal Lag (days)") +
  #       scale_color_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       scale_fill_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       theme(
  #         legend.position = "bottom", legend.text = element_text(size = 14),
  #         axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
  #         title = element_text(size = 16), strip.background = element_blank()
  #       ) +
  #       ylim(-1, 1) +
  #       # xlim(-5, 2) +
  #       # scale_x_continuous(breaks = c(-4,-3,-2,-1, 0), labels = c(-4,-3,-2,-1, 0),limits=(c(-4.5,0.5))) +
  #       scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
  #       ggtitle(title)
  #     p
  #   }
  #
  #   if (!(pathogen %in% c("rhino", "covid", "adeno"))) {
  #     actual_data4 <- actual_data3 %>%
  #       filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
  #
  #     sum_df <- actual_data4 %>%
  #       group_by(month, mobility_metric, mean_cor) %>%
  #       summarise(
  #         obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
  #         obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
  #         obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
  #       ) %>%
  #       ungroup()
  #     # filter(obs_max_ccf_lag.mean < 0.5)
  #
  #     p <- ggplot() +
  #       geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
  #       geom_point(
  #         data = sum_df,
  #         aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
  #       ) +
  #       geom_hline(yintercept = 0, lty = "dashed") +
  #       geom_vline(xintercept = 0, lty = "dashed") +
  #       facet_wrap(~ as.factor(month), nrow = 1) +
  #       theme_bw(base_size = 16) +
  #       ylab("Cross-Correlation Coefficient") +
  #       xlab("Temporal Lag (days)") +
  #       scale_color_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       scale_fill_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       theme(
  #         legend.position = "bottom", legend.text = element_text(size = 14),
  #         axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
  #         title = element_text(size = 16), strip.background = element_blank()
  #       ) +
  #       ylim(-1, 1) +
  #       # xlim(-5, 2) +
  #       # scale_x_continuous(breaks = c(-4,-3,-2,-1, 0), labels = c(-4,-3,-2,-1, 0),limits=(c(-4.5,0.5))) +
  #       scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
  #       ggtitle(title)
  #   }
  #
  #   ## SC2 has shorter time window of circulation, compared to endemic pathogens
  #   if (pathogen %in% c("covid")) {
  #     actual_data4 <- actual_data3 %>% filter(month %in% c("Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
  #     sum_df <- actual_data4 %>%
  #       group_by(month, mobility_metric, mean_cor) %>%
  #       summarise(
  #         obs_max_ccf_lag.mean = ci(obs_max_ccf_lag, na.rm = T)[1],
  #         obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag, na.rm = T)[2],
  #         obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag, na.rm = T)[3]
  #       ) %>%
  #       ungroup()
  #     # filter(obs_max_ccf_lag.mean < 0.5)
  #
  #     p <- ggplot() +
  #       geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
  #       geom_point(
  #         data = sum_df,
  #         aes(x = obs_max_ccf_lag.mean, y = mean_cor, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
  #       ) +
  #       geom_hline(yintercept = 0, lty = "dashed") +
  #       geom_vline(xintercept = 0, lty = "dashed") +
  #       facet_wrap(~ as.factor(month), nrow = 1) +
  #       theme_bw(base_size = 16) +
  #       ylab("Cross-Correlation Coefficient") +
  #       xlab("Temporal Lag (days)") +
  #       scale_color_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       scale_fill_manual(values = color_vec, labels = moblabs, breaks = unique(sum_df$mobility_metric), name = "mobility metric") +
  #       theme(
  #         legend.position = "bottom", legend.text = element_text(size = 14),
  #         axis.text.y = element_text(size = 12), strip.text = element_text(size = 16),
  #         title = element_text(size = 16), strip.background = element_blank()
  #       ) +
  #       ylim(-1, 1) +
  #       # xlim(-5, 2) +
  #       # scale_x_continuous(breaks = c(-4,-3,-2,-1, 0), labels = c(-4,-3,-2,-1, 0),limits=(c(-4.5,0.5))) +
  #       scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
  #       ggtitle(title)
  #   }
  #
  #   p
  # }

  output_s <- {
    if (pathogen %in% c("rhino", "adeno")) {
      actual_data4 <- actual_data3 %>%
        filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))

      sum_df <- actual_data4 %>%
        filter(start_week < as.Date("2021-06-01")) %>%
        group_by(month, mobility_metric, mean_cor_sp) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag_sp, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag_sp, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag_sp, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor_sp, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
        ylim(-1, 1) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        ggtitle(title)
      p
    }

    if (!(pathogen %in% c("rhino", "covid", "adeno"))) {
      actual_data4 <- actual_data3 %>%
        filter(month %in% c("Sep 2019", "Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))

      sum_df <- actual_data4 %>%
        group_by(month, mobility_metric, mean_cor_sp) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag_sp, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag_sp, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag_sp, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor_sp, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
        ylim(-1, 1) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        ggtitle(title)
    }

    ## SC2 has shorter time window of circulation, compared to endemic pathogens
    if (pathogen %in% c("covid")) {
      actual_data4 <- actual_data3 %>% filter(month %in% c("Jan 2020", "Feb 2020", "Mar 2020", "Apr 2020", "May 2020"))
      sum_df <- actual_data4 %>%
        group_by(month, mobility_metric, mean_cor_sp) %>%
        summarise(
          obs_max_ccf_lag.mean = ci(obs_max_ccf_lag_sp, na.rm = T)[1],
          obs_max_ccf_lag.lowCI = ci(obs_max_ccf_lag_sp, na.rm = T)[2],
          obs_max_ccf_lag.hiCI = ci(obs_max_ccf_lag_sp, na.rm = T)[3]
        ) %>%
        ungroup()

      p <- ggplot() +
        geom_rect(data = sum_df, aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), alpha = 0.02, fill = "yellow") +
        geom_point(
          data = sum_df,
          aes(x = obs_max_ccf_lag.mean, y = mean_cor_sp, fill = mobility_metric), size = 5, pch = 21, alpha = 0.8
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
        ylim(-1, 1) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        ggtitle(title)
    }

    p
  }

  return(list(output_s, actual_data3, actual_data4, sum_df))
}
##########################################
## Influenza B
ivb_mv_win_plot <- path_ccf_mv_window(pathogen = "ivb", l = 0.077)
ivb_mv_win_plot[[1]]

##########################################
## Influenza A/H1N1
h1n1_mv_win_plot <- path_ccf_mv_window(pathogen = "h1n1", l = 0.077)
h1n1_mv_win_plot[[1]]

##########################################
## RSV B
rsvb_mv_win_plot <- path_ccf_mv_window(pathogen = "rsv_b", l = 0.077)
rsvb_mv_win_plot[[1]]

##########################################
## RSV A
rsva_mv_win_plot <- path_ccf_mv_window(pathogen = "rsv_a", l = 0.077)
rsva_mv_win_plot[[1]]
##########################################
## HPIV 3 +4
hpiv_mv_win_plot <- path_ccf_mv_window(pathogen = "hpiv_3_4", l = 0.077)
hpiv_mv_win_plot[[1]]
##########################################
## HPIV 1 +2: outbreak ended in mid-feb 2020 before start of covid npis (don't include in figure)
# hpiv_1_2_mv_win_plot = path_ccf_mv_window(pathogen = "hpiv_1_2",l=0.077)
##########################################
## hMPV
hmpv_mv_win_plot <- path_ccf_mv_window(pathogen = "hmpv", l = 0.077)
hmpv_mv_win_plot[[1]]
##########################################
## hCoV 229E + OC43
scov_229E_OC43_mv_win_plot <- path_ccf_mv_window(pathogen = "scov_229E_OC43", l = 0.077)
scov_229E_OC43_mv_win_plot[[1]]
##########################################
## hCoV HKU1 + NL63
scov_HKU1_NL63_mv_win_plot <- path_ccf_mv_window(pathogen = "scov_HKU1_NL63", l = 0.077)
scov_HKU1_NL63_mv_win_plot[[1]]
##########################################
## SC2
covid_mv_win_plot <- path_ccf_mv_window(pathogen = "covid", l = 0.077)
covid_mv_win_plot[[1]]
##########################################
## rhinovirus
rhino_mv_win_plot <- path_ccf_mv_window(pathogen = "rhino", l = 0.077)
rhino_mv_win_plot[[1]]
##########################################
## adenovirus
adeno_mv_win_plot <- path_ccf_mv_window(pathogen = "adeno", l = 0.077)
adeno_mv_win_plot[[1]]
##########################################
## enterovirus
entero_mv_win_plot <- path_ccf_mv_window(pathogen = "entero", l = 0.077)
entero_mv_win_plot[[1]]
##########################################
## Figure 4: Combine plots for all pathogens
##########################################

##########################################
## Spearman CCF Exp Decay
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
  entero_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  adeno_mv_win_plot[[1]] + theme(legend.position = "none") + ylab("Cross-Correlation") + xlab("Temporal Lag (weeks)"),
  nrow = 6, ncol = 2, byrow = F
)
com
# leg <- get_legend(h1n1_mv_win_plot[[1]] + theme(legend.direction = "horizontal"))

leg <- cowplot::get_plot_component(h1n1_mv_win_plot[[1]] +
                                     guides(color = "none") +
                                     theme(
                                       legend.position = "bottom",
                                       legend.direction = "horizontal",
                                       legend.justification = "center",
                                       legend.box.just = "bottom",
                                       legend.text = element_text(size = 14),
                                       legend.title = element_text(size = 16)),
                                   'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(leg)
plot_component_names(h1n1_mv_win_plot[[1]])

covid <- covid_mv_win_plot[[1]] + ylab("Cross-Correlation") + theme(legend.position = "none") + ggtitle("SARS-CoV-2") + xlab("Temporal Lag (weeks)")
comb <- com + inset_element(covid, right = 1, bottom = 0, top = 0.165, left = 0.705)
comb

com2 <- plot_grid(comb, leg, nrow = 2, rel_heights = c(7, 0.3))
com2

save_plot(com2, filename = "figures/fig_4_CCF_pathogen_rt_and_mobility_by_month_2019_2020_spearman_exp_decay_select_indicators.png", base_width = 24, base_height = 16)
