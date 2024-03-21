ts_block_bootstrap_sp <- function(df1 = comb_weekly, pathogen = "hpiv_3_4", window = 10, l = 0.077, perms = 1000) {
  
  weeks <- unique(df1$epi_date)[window:(length(unique(df1$epi_date)) - window)]

  df <- expand.grid(X = metrics, Y = weeks) %>% arrange(X)

  y <- function(X, Y, l = 0.077) {
    df <- df1[, c("epi_date", as.character(X), pathogen)] %>%
      filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window))
    start <- min(df$epi_date)
    end <- max(df$epi_date)

    t1 <- df[, 2] %>% pull(1)
    t2 <- df[, 3] %>% pull(1)

    t1 <- rank(t1)
    t2 <- rank(t2)

    n <- length(t1)
    exp_w <- exp_decay_func(n, l)

    # exp_plot <- ggplot() +
    #   geom_line(aes(x = seq(1:length(exp_w)), y = exp_w))
    # exp_plot

    res <- manual_ccf(a = t1, b = t2, lag.max = 4, w = exp_w)
    res <- res %>% filter(lag < 1)

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
  result <- mapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, SIMPLIFY = F)

  actual_data <- do.call(rbind.data.frame, result)
  actual_data$start_week <- as.Date(actual_data$start_week)

  ###########################################
  ## Null Distribution CCFs
  ###########################################

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
    df2 <- x %>% dplyr::filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window))
    return(df2)
  }

  df <- expand.grid(X = metrics, Y = weeks) %>% arrange(X)
  y <- function(X = metrics, Y = weeks, perms = perms, l = l) {
    actual_inc <- df1[, c("epi_date", pathogen)] %>%
      filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window)) %>%
      as.data.frame()

    mob_df <- df1[, c("epi_date", X)] %>%
      filter(epi_date >= lubridate::as_date(Y) - (7 * window) & epi_date <= lubridate::as_date(Y) + (7 * window)) %>%
      as.data.frame()

    t1 <- mob_df[, 2]
    t2 <- actual_inc[, 2]

    t1 <- rank(t1) # mobility
    t2 <- rank(t2) # incidence

    n <- length(t1)
    exp_w <- exp_decay_func(n, l)

    set.seed(123)
    ts_boot <- tseries::tsbootstrap(x = t1, nb = perms, statistic = NULL, type = "block", b = 2) # shuffle mobility
    l <- lapply(seq_len(ncol(ts_boot)), function(i) ts_boot[, i])
    list_df <- lapply(l, function(x) data.frame(epi_date = mob_df[, 1], mob = x))

    filtered_df <- lapply(list_df, function(x) x %>% filter_epi(Y = Y, window = window))
    models <- lapply(filtered_df, function(x) manual_ccf(a = x[, 2], b = t2, lag.max = 4, w = exp_w))
    output_list <- lapply(models, FUN = extract_ccf)
    output_df <- do.call(rbind.data.frame, output_list)
    return_df <- data.frame(mobility_metric = X, start_week = as.character(Y), output_df)
    return(return_df)
  }
  tmp <- as.data.frame(t(df))

  mc.cores <- parallelly::availableCores()
  system.time(result_list <- mcmapply(FUN = y, X = tmp[1, ], Y = tmp[2, ], l = l, perms = perms, SIMPLIFY = F, mc.cores = mc.cores))
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
  return(actual_data_and_perm)
}
