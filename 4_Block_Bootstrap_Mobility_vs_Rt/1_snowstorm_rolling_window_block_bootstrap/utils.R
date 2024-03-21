
extract_Rt <- function(x) {
  post_rt <- posterior_rt(x)
  post_rt_quantiles <- get_quantiles(post_rt, levels = c(90))
  
  post_rt_median <- data.frame(
    date = post_rt$time,
    median = apply(post_rt$draws, 2, function(x) quantile(x, 0.5)),
    group = post_rt$group
  )
  
  df <- left_join(post_rt_quantiles, post_rt_median, by = c("date", "group"))
  return(df)
}

getLOESSCases <- function(df, days_incl = 21, degree = 1) {
  
  count_data <- df %>%
    mutate(confirm = round(1000* all_sites_per_pos_ili_scaled_sum))%>%
    mutate(confirm = zoo::na.approx(confirm, maxgap = 3, na.rm = F)) %>%
    mutate(confirm = tidyr::replace_na(confirm, 0)) %>%
    dplyr::select(date, confirm)

  n_points <- length(unique(count_data$date))
  dates <- unique(count_data$date)
  
  sel_span <- days_incl / n_points
  n_pad <- round(length(count_data$confirm) * sel_span * 0.5)
  c_data <- data.frame(value = c(rep(0, n_pad), count_data$confirm),
                       date_num = c(seq(as.numeric(dates[1]) - n_pad, as.numeric(dates[1]) - 1), as.numeric(dates)))
  
  c_data.lo <- loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  
  smoothed <- predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <- raw_smoothed_counts * sum(count_data$confirm, na.rm = T) / sum(raw_smoothed_counts, na.rm = T)
  
  smoothed_df = data.frame(date=dates,confirm=round(normalized_smoothed_counts))
  
  return(smoothed_df)
}

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
